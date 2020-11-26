workflow chr_qc {

    String name
    Int chr
    File vcf_loc
    Array[String] vcfs = read_lines(vcf_loc)
    File exclude_samples_loc
    Array[String] exclude_samples = read_lines(exclude_samples_loc)
    String panel_vcf
    String include_regex

    scatter (i in range(length(vcfs))) {
        call vcf_to_bed_chr {
            input: chr=chr, vcf=vcfs[i], exclude_samples=exclude_samples[i], include_regex=include_regex
        }
        call compare_panel {
            input: chr=chr, bed=vcf_to_bed_chr.out_bed, panel_vcf=panel_vcf
        }
    }

    call joint_qc {
        input: chr=chr, outname=name+"_chr"+chr, beds=vcf_to_bed_chr.out_bed, bims=vcf_to_bed_chr.out_bim, fams=vcf_to_bed_chr.out_fam
    }

    output {
        Array[File] out_beds = vcf_to_bed_chr.out_bed
        Array[File] out_bims = vcf_to_bed_chr.out_bim
        Array[File] out_fams = vcf_to_bed_chr.out_fam
        File exclude_variants = joint_qc.exclude_variants
        File snpstats = joint_qc.snpstats
        Array[File] glm = compare_panel.glm
        Array[File] freq = compare_panel.freq
        Array[File] freq_panel = compare_panel.freq_panel
        Array[File] exclude_variants_panel = compare_panel.exclude_variants
    }
}

task vcf_to_bed_chr {

    Int chr
    String chrstr = if chr == 23 then "X" else chr
    File vcf
    String base = sub(sub(basename(vcf, ".vcf"), "_original", ""), ".calls", "")
    String basechr = base + "_chr" + chr
    File exclude_samples
    String include_regex
    File ref_fasta

    command <<<

        set -euxo pipefail

        catcmd="cat"
        if [[ ${vcf} == *.gz ]] || [[ ${vcf} == *.bgz ]]; then catcmd="zcat"; fi
        mem=$((`free -m | grep -oP '\d+' | head -n 1`-500))

        # filter by chr, PASS variants and non-CNVs, rename variant id to chr1_1234_A_T
        echo -e "`date`\tfilter to chr"
        $catcmd ${vcf} | awk 'BEGIN{OFS="\t"} $1 ~ "^#" {print $0} $1 ~ "^(chr)?${chrstr}$" && $5 !~ "CNV" && $7 == "PASS" {sub("chr", "", $1); $1="chr"$1; $3=$1"_"$2"_"$4"_"$5; print $0}' | bgzip -@2 > vcf.gz
        tabix -p vcf vcf.gz

        # filter by given samples and prefix, align to reference, remove multiallelics, rename variant id to chr1_1234_A_T
        echo -e "`date`\talign"
        comm -23 <($catcmd ${vcf} | grep -E "^#CHROM"  | head -1 | tr '\t' '\n' | tail -n+10 | grep -E "${include_regex}" | sort) <(cut -f1 ${exclude_samples} | sort) > include_samples.txt
        bcftools view -S include_samples.txt vcf.gz -Ou | \
        bcftools norm -f ${ref_fasta} -c ws -Ou | \
        bcftools view -m 2 -M 2 -Ou | \
        bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' -Oz -o vcf_anno.vcf.gz

        # convert to plink
        echo -e "`date`\tconvert to plink"
        plink2 --allow-extra-chr --memory $mem --vcf vcf_anno.vcf.gz --make-bed --out ${basechr}
        echo -e "`date`\tdone"

    >>>

    output {
        File out_bed = basechr + ".bed"
        File out_bim = basechr + ".bim"
        File out_fam = basechr + ".fam"
    }

    runtime {
        docker: "eu.gcr.io/finngen-refinery-dev/bioinformatics:0.6"
        memory: "3 GB"
        cpu: 4
        disks: "local-disk 200 HDD"
        preemptible: 2
    }
}

task compare_panel {

    Boolean compare
    Int chr
    String chrstr = if chr == 23 then "X" else chr
    File bed
    File bim = sub(bed, ".bed$", ".bim")
    File fam = sub(bed, ".bed$", ".fam")
    String base = basename(bed, ".bed")
    File panel_vcf
    String base_panel = basename(panel_vcf, ".vcf.gz")
    File panel_freq
    String pca_ld
    Float pca_maf
    Float af_panel
    Float af_diff
    Float af_fc

    command <<<

        set -euxo pipefail

        touch ${base}.exclude_variants.txt
        touch ${base}.PHENO1.glm.logistic.hybrid
        touch ${base}_chip.frq
        touch ${base_panel}.frq
        touch ${base}_chip_AFs.png
        touch ${base}.eigenvec
        if ! ${compare}; then
            exit 0
        fi

        mem=$((`free -m | grep -oP '\d+' | head -n 1`-500))
        plink_cmd="plink --allow-extra-chr --keep-allele-order --memory $mem --chr ${chr}"
        plink2_cmd="plink2 --allow-extra-chr --memory $mem --chr ${chr}"

        # compare to panel and create af plots
        echo -e "`date`\tcompare to panel r"
        $plink2_cmd --bfile ${sub(bed, ".bed$", "")} --freq --out ${base}
        awk 'BEGIN{OFS="\t"} NR==1{print "CHR","SNP","REF","ALT","AF"} NR>1{print "chr"$1,$2,$3,$4,$5}' ${base}.afreq > ${base}_chip.frq
        awk 'NR==1 || $1 ~ "^(chr)?${chrstr}$"' ${panel_freq} > ${base_panel}.frq
        Rscript /tools/r_scripts/plot_AF_FC_v1.3.3.R ${base} chip ${base_panel}.frq ${af_diff} ${af_fc}

        # get variants to exclude from panel comparison
        awk 'BEGIN{OFS="\t";} {print($1, "not_in_panel", "NA")}' ${base}_chip_nonpanel_exclude.txt >> ${base}.exclude_variants.txt
        awk 'BEGIN{OFS="\t";} {print($1, "af_diff_panel", "NA")}' ${base}_chip_exclude.txt >> ${base}.exclude_variants.txt
        awk 'BEGIN{OFS="\t";} FNR==NR&&NR>1{a[$2]=1} FNR<NR&&NR>1&&$2 in a &&$5<${af_panel} {print($2, "af_panel", $5)}' ${base}_chip.frq ${base_panel}.frq >> ${base}.exclude_variants.txt

        # create plink files for dataset and panel with common variants between them
        echo -e "`date`\tcompare to panel glm"
        $plink2_cmd --bfile ${sub(bed, ".bed$", "")} --exclude <(cut -f1 ${base}.exclude_variants.txt | sort -u) --make-bed --out panel_isec
        # rename dataset ids in case there are panel samples
        awk '{OFS="\t"; $2="DATA-"$2; print $0}' panel_isec.fam > temp && mv temp panel_isec.fam
        awk 'FNR==NR{a[$2]=1} FNR<NR && ($0 ~ "^#" || $3 in a)' panel_isec.bim <(zcat ${panel_vcf}) | bgzip > panel.vcf.gz
        $plink2_cmd --vcf panel.vcf.gz --make-bed --out panel
        awk 'BEGIN{OFS="\t"} {print "0", $2, "1"}' panel.fam > pheno
        awk 'BEGIN{OFS="\t"} {print "0", $2, "2"}' panel_isec.fam >> pheno

        # merge panel and data and run logreg between them
        $plink_cmd --bfile panel_isec --bmerge panel --make-bed --out merged
        if [[ ${chr} == 23 ]]; then
            $plink2_cmd --bfile merged --const-fid --glm firth-fallback allow-no-covars --pheno pheno --out ${base}
        else
            $plink2_cmd --bfile merged --maf ${pca_maf} --indep-pairwise ${pca_ld} --out pruned
            $plink2_cmd --bfile merged --pca --extract pruned.prune.in --out ${base}
            $plink2_cmd --bfile merged --const-fid --glm firth-fallback --covar ${base}.eigenvec --pheno pheno --out ${base}
            grep -E "TEST|ADD" ${base}.PHENO1.glm.logistic.hybrid > temp && mv temp ${base}.PHENO1.glm.logistic.hybrid
        fi
        $plink2_cmd --bfile panel_isec --freq --out ${base}
        $plink2_cmd --bfile panel --freq --out ${base_panel}
        awk 'BEGIN{OFS="\t"} NR==1{print "CHR","SNP","REF","ALT","AF"} NR>1{print "chr"$1,$2,$3,$4,$5}' ${base}.afreq > ${base}.frq
        awk 'BEGIN{OFS="\t"} NR==1{print "CHR","SNP","REF","ALT","AF"} NR>1{print "chr"$1,$2,$3,$4,$5}' ${base_panel}.afreq > ${base_panel}.frq
        echo -e "`date`\tdone"

    >>>

    output {
        File exclude_variants = base + ".exclude_variants.txt"
        File glm  = base + ".PHENO1.glm.logistic.hybrid"
        File freq = base + ".frq"
        File freq_panel = base_panel + ".frq"
        File eigenvec = base + ".eigenvec"
        File png = base + "_chip_AFs.png"
    }

    runtime {
        docker: "eu.gcr.io/finngen-refinery-dev/pftools:0.1.2"
        memory: "7 GB"
        cpu: 1
        disks: "local-disk 200 HDD"
        preemptible: 2
    }
}

task joint_qc {

    Boolean run_joint_qc
    Int chr
    String outname
    Array[File] beds
    Array[File] bims
    Array[File] fams
    Float af
    Float hw
    Float variant_missing
    String dollar = "$"

    command <<<

        set -euxo pipefail

        touch ${outname}.snpstats.txt
        touch ${outname}.exclude_variants.txt
        if ! ${run_joint_qc}; then
            exit 0
        fi

        # move plink files to current directory
        for file in ${sep=" " beds}; do
            mv $file .; mv ${dollar}{file/\.bed/\.bim} .; mv ${dollar}{file/\.bed/\.fam} .
        done

        # merge plink datasets
        echo -e "`date`\tmerge plink data"
        ls -1 *.bed | sed 's/\.bed//' > mergelist.txt
        mem=$((`free -m | grep -oP '\d+' | head -n 1`-500))
        plink --allow-extra-chr --keep-allele-order --memory $mem --merge-list mergelist.txt --make-bed --out ${outname}

        # calculate snpstats and get variants to exclude by af, hwe, missingness
        echo -e "`date`\tcalculate snpstats"
        qctool -g ${outname}.bed -snp-stats -osnp ${outname}.snpstats.txt
        grep -Ev "^#" ${outname}.snpstats.txt | awk '
        BEGIN {OFS="\t"} NR==1 {for (i=1;i<=NF;i++) h[$i]=i}
        NR>1 {
            id=$h["rsid"]; af1=$h["alleleA_frequency"]; af2=$h["alleleB_frequency"]; hw=$h["HW_exact_p_value"]; miss=$h["missing_proportion"];
            if (af1<${af}) print id,"alleleA_frequency",af1;
            if (af2<${af}) print id,"alleleB_frequency",af2;
            if ($h["chromosome"] !~ "(chr)?X|23" && hw<${hw}) print id,"HW_exact_p_value",hw;
            if (miss>${variant_missing}) print id,"missing_proportion",miss;
        }' > ${outname}.exclude_variants.txt
        echo -e "`date`\tdone"

    >>>

    output {
        File snpstats = outname + ".snpstats.txt"
        File exclude_variants = outname + ".exclude_variants.txt"
    }

    runtime {
        docker: "eu.gcr.io/finngen-refinery-dev/bioinformatics:0.6"
        memory: "24 GB"
        cpu: 6
        disks: "local-disk 500 HDD"
        preemptible: 2
    }
}
