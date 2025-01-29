workflow imputation {

    Int chr
    Array[String] beds
    Array[String] batch_qc_exclude_variants
    Array[String] panel_exclude_variants
    Map[Int, String] force_impute_variants
    String exclude_samples
    Map[Int, String] ref_panel
    Boolean add_batch_suffix = true
    String docker

    scatter (i in range(length(beds))) {
        call plink_to_vcf {
            input: chr=chr, bed=beds[i],
            batch_qc_exclude_variants=batch_qc_exclude_variants[i],
            panel_exclude_variants=panel_exclude_variants[i],
            exclude_samples=exclude_samples,
            force_impute_variants=force_impute_variants[chr],
            docker=docker
        }
        call phase_impute {
            input: chr=chr, vcf=plink_to_vcf.vcf,
            ref_panel=ref_panel, docker=docker
        }
        call post_imputation {
            input: chr=chr, vcf=phase_impute.out_imputed,
            add_batch_suffix = add_batch_suffix,
            docker=docker
        }
    }

    output {
        Array[File] vcfs = post_imputation.tags_edited
        Array[File] vcf_idxs = post_imputation.tags_edited_tbi
    }
}

task plink_to_vcf {

    Int chr
    File bed
    File bim = sub(bed, ".bed$", ".bim")
    File fam = sub(bed, ".bed$", ".fam")
    String base = basename(bed, ".bed")
    File batch_qc_exclude_variants
    File panel_exclude_variants
    File exclude_samples
    File force_impute_variants
    String docker

    command <<<

        set -euxo pipefail

        mem=$((`free -m | grep -oP '\d+' | head -n 1`-500))
        plink2_cmd="plink2 --allow-extra-chr --memory $mem --chr ${chr}"

        # get list of variants to exclude
        cat \
        <(cut -f1 ${batch_qc_exclude_variants}) \
        <(cut -f1 ${panel_exclude_variants}) ${force_impute_variants} | \
        sort -u > exclude_variants.txt

        ## Removing postfix
        awk -v base=${base} ' BEGIN{ batch=base; sub("_chr..?$","",batch)} $1==batch && !seen[$2]++ {print 0,$2}' ${exclude_samples} > exclude_samples.txt

        $plink2_cmd --bfile ${sub(bed, ".bed$", "")} --remove exclude_samples.txt --exclude exclude_variants.txt --make-bed --out temp
        # remove _dup suffix from fam
        sed -i 's/_dup[0-9]\+//' temp.fam
        # underscores cause trouble with vcf conversions
        sed -i 's/_/-/g' temp.fam
        $plink2_cmd --bfile temp --recode vcf-iid bgz --output-chr chrM --remove exclude_samples.txt --exclude exclude_variants.txt --out temp
        $plink2_cmd --vcf temp.vcf.gz --recode vcf-iid bgz --vcf-half-call haploid --output-chr chrM --out ${base}

        bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' ${base}.vcf.gz -Ou | \
        bcftools +fill-tags -Oz -o ${base}_chr${chr}_qcd.vcf.gz -- -t AC,AN,AF

    >>>

    output {
        File vcf = base + "_chr" + chr + "_qcd.vcf.gz"
        Array[File] png = glob("*.png")
    }
    runtime {
        docker: "${docker}"
        memory: "7 GB"
        cpu: 4
        disks: "local-disk 200 HDD"
        preemptible: 2
        zones: "europe-west1-b europe-west1-c europe-west1-d"
    }
}

task phase_impute {

    Int chr
    String chrname = if chr == 23 then "chrX" else chr
    File vcf
    String base = basename(vcf, ".vcf.gz")
    Map[Int, String] ref_panel
    Map[Int, String] genetic_maps_eagle
    Map[Int, String] genetic_maps_beagle
    File panel = ref_panel[chr]
    File genetic_map_eagle = genetic_maps_eagle[chr]
    File genetic_map_beagle = genetic_maps_beagle[chr]
    String docker
    String dollar = "$"

    command <<<

        n_cpu=`grep -c ^processor /proc/cpuinfo`
        g_mem=$((`free -g | grep -oP '\d+' | head -n 1`-1))

        eagle \
            --vcf            ${vcf} \
            --chrom          ${chrname} \
            --geneticMapFile ${genetic_map_eagle} \
            --outPrefix      ${base}_for_imputation \
            --Kpbwt          20000 \
            --numThreads     $n_cpu

        java -Xss5m -Xmx${dollar}{g_mem}g -jar /tools/beagle/beagle.27Jan18.7e1.jar \
            gt=${base}_for_imputation.vcf.gz \
            ref=${panel} \
            map=${genetic_map_beagle} \
            out=${base}_imputed \
            niterations=10 \
            ne=20000 \
            impute=true \
            gprobs=true \
            seed=-99999 \
            nthreads=$n_cpu

    >>>

    output {
        File out_phased = base + "_for_imputation.vcf.gz"
        File out_imputed = base + "_imputed.vcf.gz"
    }

    runtime {
        docker: "${docker}"
        memory: "24 GB"
        cpu: 32
        disks: "local-disk 100 HDD"
        preemptible: 2
        zones: "europe-west1-b europe-west1-c europe-west1-d"
    }
}

task post_imputation {

    String chr
    File vcf
    String base = basename(vcf, ".vcf.gz")
    File ref_panel_freq
    File annot_hdr
	File annot_tab
	File annot_tab_index
	String annot_col_incl
    Boolean add_batch_suffix
    String docker
    String dollar = "$"

    command <<<

        set -euxo pipefail

        # add info tags
        bcftools index -t -f ${vcf}
        bcftools +fill-tags ${vcf} -Ou -- -t AF,AC_Hom,AC_Het,NS,HWE | \
        bcftools +impute-info -Oz -o ${base}_imputed_infotags.vcf.gz

        # create info/af plots
        echo -e 'CHR\tSNP\tREF\tALT\tAF\tINFO\tAF_GROUP' > ${base}_varID_AF_INFO_GROUP.txt
        time bcftools query -f '%CHROM\t%CHROM\_%POS\_%REF\_%ALT\t%REF\t%ALT\t%INFO/AF\t%INFO/INFO\t-\n' ${base}_imputed_infotags.vcf.gz | \
        awk '{if ($5>=0.05 && $5<=0.95) $7=1; else if(($5>=0.005 && $5<0.05) || ($5<=0.995 && $5>0.95)) $7=2; else $7=3} { print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7 }' \
        >> ${base}_varID_AF_INFO_GROUP.txt
        time Rscript --no-save /tools/r_scripts/plot_INFO_and_AF_for_imputed_chrs.R ${base} ${ref_panel_freq} ${chr} ${base}_varID_AF_INFO_GROUP.txt

        # annotate, edit tags and index
        edits="AR2: DR2: IMP:"
        if ${add_batch_suffix}; then
            [[ ${base} =~ (.*)_qcd_imputed ]]
            batch=${dollar}{BASH_REMATCH[1]}
            edits="AF:AF_$batch INFO:INFO_$batch CHIP:CHIP_$batch AC_Het:AC_Het_$batch AC_Hom:AC_Hom_$batch HWE:HWE_$batch NS:NS_$batch AR2: DR2: IMP:"
        fi
        time bcftools index -t -f ${base}_imputed_infotags.vcf.gz
		time bcftools annotate -i 'INFO/IMP=0' -k -a ${annot_tab} -h ${annot_hdr} -c ${annot_col_incl} ${base}_imputed_infotags.vcf.gz -Ou | \
		bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' -Ov | \
        java -jar /tools/tagEditing_v1.1.jar $edits | \
        bgzip > ${base}_tags_edited.vcf.gz
        time bcftools index -t -f ${base}_tags_edited.vcf.gz

    >>>

    runtime {
        docker: "${docker}"
        memory: "16 GB"
        cpu: 2
        disks: "local-disk 200 HDD"
        preemptible: 2
        zones: "europe-west1-b europe-west1-c europe-west1-d"
    }

    output {
        File out_varid_af_info_group = "${base}_varID_AF_INFO_GROUP.txt"
        File tags_edited = base + "_tags_edited.vcf.gz"
        File tags_edited_tbi = base + "_tags_edited.vcf.gz.tbi"
        Array[File] png = glob("*.png")
    }
}
