import "qc_sub.wdl" as qc_sub
import "imp_sub.wdl" as imp_sub

workflow qc_imputation {

    String name
    Boolean run_imputation
    File vcf_loc
    Array[String] vcfs = read_lines(vcf_loc)
    File fam_loc
    Array[String] fams_sexcheck = read_lines(fam_loc)
    File exclude_samples_loc
    Array[String] exclude_samples = read_lines(exclude_samples_loc)
    Map[Int, String] ref_panel
    Array[Int] chrs
    String include_regex
    Int genome_build
    Float variant_missing
    Float sample_missing
    Array[Float] f
    Float het_sd
    Int pi_hat_min_n
    Int pi_hat_min_n_excess

    # run joint qc and panel comparison per chromosome
    scatter (chr in chrs) {
        call qc_sub.chr_qc as chr_qc {
            input: name=name, chr=chr, panel_vcf=ref_panel[chr], include_regex=include_regex
        }
    }
    # run qc per batch across all chromosomes
    scatter (i in range(length(vcfs))) {
        String base1 = sub(sub(basename(vcfs[i], ".vcf"), "_original", ""), ".calls", "")
        # combine panel comparison across chrs and create plot
        call panel_comparison {
            input: base=base1, freq=transpose(chr_qc.freq)[i], freq_panel=transpose(chr_qc.freq_panel)[i],
            files_glm=transpose(chr_qc.glm)[i], exclude_variants_panel=transpose(chr_qc.exclude_variants_panel)[i]
        }
        # run qc
        call batch_qc {
            input: base=base1, beds=transpose(chr_qc.out_beds)[i], bims=transpose(chr_qc.out_bims)[i], fams=transpose(chr_qc.out_fams)[i],
            fam_sexcheck=fams_sexcheck[i], exclude_variants_joint=chr_qc.exclude_variants, exclude_variants_panel=panel_comparison.exclude_variants,
            genome_build=genome_build, exclude_samples=exclude_samples[i], variant_missing=variant_missing, sample_missing=sample_missing,
            f=f, het_sd=het_sd, pi_hat_min_n=pi_hat_min_n, pi_hat_min_n_excess=pi_hat_min_n_excess, include_regex=include_regex
        }
    }
    # get duplicates to remove across batches
    call duplicates {
        input: samples_include=batch_qc.samples_include, samples_exclude=batch_qc.samples_exclude, sample_missingness=batch_qc.sample_missingness
    }
    # create plots per batch
    scatter (i in range(length(vcfs))) {
        String base2 = sub(sub(basename(vcfs[i], ".vcf"), "_original", ""), ".calls", "")
        call plots {
            input: name=base2, joint=false, sexcheck=[batch_qc.sexcheck[i]], heterozygosity=[batch_qc.heterozygosity[i]], exclude_samples=duplicates.allbatches_samples_exclude,
            exclude_variants=[batch_qc.variants_exclude[i]], variant_missingness=[batch_qc.variant_missingness[i]], sample_missingness=[batch_qc.sample_missingness[i]], pihat_n=[batch_qc.pihat_n[i]],
            pihat_n_raw=[batch_qc.pihat_n_raw[i]], variant_missing=variant_missing, sample_missing=sample_missing, f=f, het_sd=het_sd, pi_hat_min_n=pi_hat_min_n
        }
    }
    # create plots across all batches
    call plots as joint_plots {
        input: name=name, joint=true, sexcheck=batch_qc.sexcheck, heterozygosity=batch_qc.heterozygosity, exclude_samples=duplicates.allbatches_samples_exclude,
        exclude_variants=batch_qc.variants_exclude, variant_missingness=batch_qc.variant_missingness, sample_missingness=batch_qc.sample_missingness, pihat_n=batch_qc.pihat_n,
        pihat_n_raw=batch_qc.pihat_n_raw, variant_missing=variant_missing, sample_missing=sample_missing, f=f, het_sd=het_sd, pi_hat_min_n=pi_hat_min_n
    }
    # run imputation per chromosome
    if (run_imputation) {
        scatter (i in range(length(chrs))) {
            call imp_sub.imputation {
                input: chr=chrs[i], beds=chr_qc.out_beds[i], joint_qc_exclude_variants=chr_qc.exclude_variants[i],
                batch_qc_exclude_variants=batch_qc.variants_exclude, exclude_samples=duplicates.allbatches_samples_exclude
            }
        }
    }
}

task panel_comparison {

    String base
    Array[File] freq
    Array[File] freq_panel
    Array[File] files_glm
    Array[File] exclude_variants_panel
    Float p

    command <<<

        awk 'FNR==NR||FNR>1' ${sep=' ' freq} > data.frq
        awk 'FNR==NR||FNR>1' ${sep=' ' freq_panel} > panel.frq
        awk 'FNR==NR||FNR>1' ${sep=' ' files_glm} | sed 's/#CHR/CHR/' > ${base}.glm

        # plot_AF_comparison.R from factory
        Rscript - <<EOF

        library(qqman)
        options(bitmapType='cairo')

        # Input variables
        output_dir <- "./"
        indataset <- "${base}" # dataset name tag
        pval <- "${base}.glm" # path to plink .hybrid file
        af <- "data.frq" # .frq file of the dataset
        af_panel <- "panel.frq" # .frq file of panel
        pval_limit <- ${p} # threshold for p-value

        # Read input files
        pval <- read.delim(pval, header = T, stringsAsFactors = F)
        af <- read.delim(af, header = T, stringsAsFactors = F)
        af_panel <- read.delim(af_panel, header = T, stringsAsFactors = F)
        print(head(pval))

        ## Dataframe modifications and objects for plotting
        # For manhattan plot
        #colnames(pval)[c(1,2)] <- c("CHR", "BP") # required for manhattan
        colnames(pval)[1] <- "CHR" # required for manhattan
        colnames(pval)[2] <- "BP"
        print(head(pval))
        pval[["CHR"]] <- sub("X", "23", pval[["CHR"]])
        pval[["CHR"]] <- as.integer(pval[["CHR"]]) # required for manhattan

        print(head(pval))
        print(head(af_panel))
        print(head(af))

        # For output writing and AF plots
        colnames(af) <- colnames(af_panel) <- c("CHROM", "ID", "REF", "ALT", "AF")
        af_panel <- af_panel[,c(2,5)]
        output <- merge(merge(af, af_panel, by ="ID"), pval, by = c("ID", "REF", "ALT"))
        colnames(output)[c(5,6,9)] <- c("AF_CHIP", "AF_PANEL", "FIRTH")

        # For qqplot
        logobspval <- -log10(sort(output[["P"]]) + 1.0E-100)
        logexppval <- -log10((1:length(logobspval)-0.5)/length(1:length(logobspval)))
        lambda <- round(median(logobspval/logexppval), 3)

        # For AF comparison to panel plot
        pval_limit_dyn <- pval_limit/(median(logobspval/logexppval)^3)
        sig <- output[output[["P"]] < pval_limit_dyn,]
        sig[["AF_CHIP"]] <- as.numeric(sig[["AF_CHIP"]])

        # Plot p-value distribution with a red line at p-value threshold, and
        # the same histogram only for significant AFs,
        # AF dataset vs. AF panel with significant as ,
        # qqplot,
        # another plot for dataset AF vs p-values with a red line at 5e-08,
        # and a manhattan plot
        png(paste0(output_dir, "/", indataset, "panel_AF_glmfirthfallback_plots.png"), width = 1200, height = 1800)
        par(mfrow = c(3,2))
        par(cex.axis = 1.6, cex.lab = 1.5, cex.main = 1.6)

        hist(-log10(output[["P"]]), breaks = 500, main = paste0(indataset, " vs. Reference panel\np-value"), xlab = "-log10(p)")
        abline(v=-log10(pval_limit_dyn), col = "red", lwd = 2)

        hist(sig[["AF_CHIP"]], breaks = 100, main = paste0(indataset, " vs. Reference panel,\nAF (p < ", pval_limit_dyn, ")"), xlab = paste0(indataset, " AF"))

        plot(output[["AF_PANEL"]], output[["AF_CHIP"]], col = 1, pch = 20, cex = 1.5,  main = paste0(indataset, " vs. Reference panel\nAF"), xlab = "Reference panel AF", ylab = paste0(indataset, " AF"))
        points(sig[["AF_PANEL"]], sig[["AF_CHIP"]], col = "red", pch = 20)
        legend("topleft", legend = c(paste0("ns., n = ", nrow(output[output[["P"]] >= pval_limit_dyn,])), paste0("p < ", round(pval_limit_dyn, 14), ", n = ", nrow(sig))), col = c("black", "red"), pch=20, cex=2)

        qq(output[["P"]], lwd = 2, main = paste0(indataset, ", lambda: ", lambda), las = 1, cex = 1.5)

        plot(output[["AF_CHIP"]], -log10(output[["P"]]), col = 1, pch = 20, cex = 1.5, main = paste0(indataset, " vs. Reference panel\nAF vs. p-value"), ylab = "-log10(p)", xlab = paste0(indataset, " AF"))
        abline(h=-log10(pval_limit_dyn), col = "red", lwd = 2)

        manhattan(output, main = indataset, cex = 1.5, suggestiveline = -log10(pval_limit_dyn))

        dev.off()

        # Store table with AF and p-value included
        write.table(output, paste0(output_dir, "/", indataset, "_panel_AF_glmfirthfallback.txt"), row.names = F, quote = F, sep = "\t")
        write.table(sig, paste0(output_dir, "/", indataset, "_panel_AF_glmfirthfallback_exclusion.txt"), row.names = F, quote = F, sep = "\t")
        write.table(sig[["ID"]], paste0(output_dir, "/", indataset, "_panel_AF_glmfirthfallback_exclusion_SNPID.txt"), row.names = F, col.names = F, quote = F, sep = "\t")

        sink(paste0(output_dir, "/", indataset, "panel_AF_glmfirthfallback_pvalue_report.txt"))
        cat(paste0("pvalue_dynamic:", pval_limit_dyn))
        cat("\n")
        cat(paste0("pvalue_defined:", pval_limit))
        cat("\n")
        cat(paste0("lambda value:", lambda))
        sink()
        EOF

        cat ${sep=' ' exclude_variants_panel} > ${base}.exclude_variants.txt
        awk 'BEGIN{OFS="\t"}NR==1{for(i=1;i<=NF;i++) a[$i]=i} NR>1 && $a["P"] < ${p} {print $1,"glm_panel",$a["P"]}' ${base}_panel_AF_glmfirthfallback_exclusion.txt >> ${base}.exclude_variants.txt

    >>>

    output {
        File glm = base + ".glm"
        File plots = base + "panel_AF_glmfirthfallback_plots.png"
        File exclude_variants = base + ".exclude_variants.txt"
    }

    runtime {
        docker: "eu.gcr.io/finngen-refinery-dev/bioinformatics:0.6"
        memory: "7 GB"
        cpu: 1
        disks: "local-disk 200 HDD"
        preemptible: 2
    }
}

task batch_qc {

    String base
    Array[File] beds
    Array[File] bims
    Array[File] fams
    Int genome_build
    Array[File] exclude_variants_joint
    File exclude_variants_panel
    File exclude_samples
    String include_regex
    File fam_sexcheck
    Boolean check_ssn_sex
    File ssn_sex
    File high_ld_regions
    Float hw
    Float variant_missing
    Float sample_missing
    Array[Float] f
    Float het_sd
    Float pi_hat
    String pihat_ld
    Int pi_hat_min_n_excess
    Int pi_hat_min_n
    String dollar = "$"

    command <<<

        set -euxo pipefail

        mem=$((`free -m | grep -oP '\d+' | head -n 1`-500))
        plink_cmd="plink --allow-extra-chr --keep-allele-order --memory $mem"
        plink2_cmd="plink2 --allow-extra-chr --memory $mem"

        # get variants to exclude based on joint qc and panel comparison
        for file in ${sep=" " exclude_variants_joint}; do cat $file >> upfront_exclude_variants.txt; done
        cat ${exclude_variants_panel} >> upfront_exclude_variants.txt

        # move per-chr plink files to current directory
        for file in ${sep=" " beds}; do
            mv $file .; mv ${dollar}{file/\.bed/\.bim} .; mv ${dollar}{file/\.bed/\.fam} .
        done

        # get samples in data to exclude by given list not counting the same sample more than once
        # TODO using sex check fam here because samples may have been excluded earlier and we want to list all excluded samples here
        awk 'NR==FNR {a[$2]=1} NR>FNR && $1 in a && !seen[$1]++' ${fam_sexcheck} ${exclude_samples} | grep -E "${include_regex}" | \
        awk 'BEGIN{OFS="\t"} {print $1,$2,"NA"}' | sort -u >> ${base}.samples_exclude.txt

        # merge chrs excluding variants that did not pass joint qc or panel comparison and samples in given list
        echo -e "`date`\tmerge"
        ls -1 *.bed | sed 's/\.bed//' > mergelist.txt
        $plink_cmd --merge-list mergelist.txt --remove <(awk '{print 0,$1}' ${base}.samples_exclude.txt) --exclude <(cut -f1 upfront_exclude_variants.txt | sort -u) --make-bed --out ${base}

        # get initial missingness for all remaining samples
        $plink2_cmd --bfile ${base} --missing --out missing
        cp missing.smiss ${base}.raw.smiss

        # check that given sex check fam file matches data
        echo -e "`date`\tsex_check"
        join -t $'\t' -1 3 -2 2 <(awk '{OFS="\t"; $1=$1; print $0}' ${base}.fam | nl -nln | sort -k3,3) <(sort -k2,2 ${fam_sexcheck}) | \
        sort -b -k2,2g | awk '{OFS="\t"; print $8,$1,$9,$10,$11,$12}' > sexcheck.fam
        if [[ `diff <(awk '{print $2}' sexcheck.fam) <(awk '{print $2}' ${base}.fam) | wc -l | awk '{print $1}'` != 0 ]]; then
            >2& echo "samples differ between ${fam_sexcheck} and ${base}.fam"
            exit 1
        fi
        # check sex excluding PAR region
        mv sexcheck.fam ${base}.fam
        $plink_cmd --bfile ${base} --split-x b${genome_build} no-fail --make-bed --out plink_data
        $plink_cmd --bfile plink_data --check-sex ${sep=" " f} --out plink_data
        awk 'BEGIN{OFS="\t"} NR>1&&$5!="OK" {print $2,"sex_check",$6}' plink_data.sexcheck >> ${base}.samples_exclude.txt
        cp plink_data.sexcheck ${base}.sexcheck
        awk 'NR==1||$5=="OK"' ${base}.sexcheck > sexcheck.ok
        # check against ssn sex for samples that passed sex check
        if ${check_ssn_sex}; then
            join -t $'\t' -1 2 -2 1 \
            <(awk 'BEGIN{OFS="\t"} NR==1 {print "IID","SNPSEX"} NR>1 {sub("_dup(.)*", "", $2); print $2,$4}' sexcheck.ok | nl -nln | sort -b -k2,2) \
            <(awk 'BEGIN{OFS="\t"} NR==1 {print "IID","SSN_SEX"} NR>1 {print $1,$2}' ${ssn_sex} | sort -k1,1) | \
            sort -k2,2g | cut -f1,3- | awk 'NR>1&&$2!=$3 {OFS="\t"; print $1,"sex_check_ssn","NA"}' >> ${base}.samples_exclude.txt
        fi

        # remove samples not passing sex check and given samples and filter by hw
        echo -e "`date`\thwe"
        $plink2_cmd --bfile ${base} --remove <(awk '{OFS="\t"; print "0",$1}' ${base}.samples_exclude.txt | sort -u) \
              --hardy --hwe ${hw} --make-bed --out plink_data_after_hwe
        awk 'NR>1&&$10<${hw} {OFS="\t"; print $2,"hwe",$10}' plink_data_after_hwe.hardy >> ${base}.variants_exclude.txt
        cp plink_data_after_hwe.hardy ${base}.hardy

        # remove variants by missingness
        echo -e "`date`\tmissingness"
        $plink2_cmd --bfile plink_data_after_hwe --geno ${variant_missing} --make-bed --out plink_data
        $plink2_cmd --bfile plink_data --missing --out missing
        # get samples to exclude by missingness
        awk 'BEGIN{OFS="\t"} NR>1 && $5>${sample_missing} {print $2,"missing_proportion",$5}' missing.smiss >> ${base}.samples_exclude.txt
        cp missing.smiss ${base}.smiss

        # remove samples by missingness and sex check
        # going back to when variants were not yet excluded by missingness
        echo -e "`date`\tfilter_samples"
        $plink2_cmd --bfile plink_data_after_hwe --remove <(awk '{OFS="\t"; print "0",$1}' ${base}.samples_exclude.txt | sort -u) --missing --make-bed --out plink_data

        # remove variants by missingness once more after removing samples by missingness
        echo -e "`date`\tmissingness"
        awk 'BEGIN{OFS="\t"} NR>1 && $5!="nan" && $5>${variant_missing} {print $2,"missing_proportion",$5}' plink_data.vmiss >> ${base}.variants_exclude.txt
        cp plink_data.vmiss ${base}.vmiss
        # NOTE variant exclusion now done, filter with high AF threshold for the rest of QC
        $plink2_cmd --bfile plink_data --maf 0.05 --exclude <(cut -f1 ${base}.variants_exclude.txt) --make-bed --out plink_data

        # get samples to exclude due to high heterozygosity (contamination)
        echo -e "`date`\theterozygosity"
        $plink_cmd --bfile plink_data --het --out plink_data
        awk 'BEGIN{OFS="\t"} NR==1{$1=$1; print $0,"het_rate"} NR>1{$1=$1; print $0,($5-$3)/$3}' plink_data.het > ${base}.heterozygosity.txt
        mean=`tail -n+2 ${base}.heterozygosity.txt | cut -f7 | datamash mean 1`
        sdev=`tail -n+2 ${base}.heterozygosity.txt | cut -f7 | datamash sstdev 1`
        awk -v mean=$mean -v sdev=$sdev 'BEGIN{OFS="\t"} NR>1&&$7!="NA"&&($7-mean>${het_sd}*sdev) {print $2,"heterozygosity",$7}' ${base}.heterozygosity.txt >> ${base}.samples_exclude.txt

        # get samples to exclude due to extreme pi-hat if any are left after heterozygosity filtering
        echo -e "`date`\tpi_hat"
        $plink2_cmd --bfile plink_data --exclude range ${high_ld_regions} --remove <(awk '{print 0,$1}' ${base}.samples_exclude.txt | sort -u) --indep-pairwise ${pihat_ld} --out pruned
        $plink_cmd --bfile plink_data --extract pruned.prune.in --remove <(awk '{print 0,$1}' ${base}.samples_exclude.txt | sort -u) --genome gz --out genome
        zcat genome.genome.gz | awk 'NR>1&&$10>${pi_hat} {n[$2]++;n[$4]++} END{for(id in n) print id,n[id]}' > ${base}.pihat_n_raw.txt
        awk 'BEGIN{OFS="\t"} $2>=${pi_hat_min_n_excess} {print $1,"pi_hat_excess",$2}' ${base}.pihat_n_raw.txt >> ${base}.samples_exclude.txt
        mv genome.genome.gz ${base}.genome_raw.gz
        # exclude samples with an extreme number of relatives based on pi-hat (contamination) and run again
        $plink_cmd --bfile plink_data --remove <(awk '{print 0,$1}' ${base}.samples_exclude.txt | sort -u) --extract pruned.prune.in --genome gz --out genome
        zcat genome.genome.gz | awk 'NR>1&&$10>${pi_hat} {n[$2]++;n[$4]++} END{for(id in n) print id,n[id]}' > ${base}.pihat_n.txt
        awk 'BEGIN{OFS="\t"} $2>=${pi_hat_min_n} {print $1,"pi_hat",$2}' ${base}.pihat_n.txt >> ${base}.samples_exclude.txt
        mv genome.genome.gz ${base}.genome.gz

        # get included samples
        comm -23 <(cut -f2 ${base}.fam | sort) <(cut -f1 ${base}.samples_exclude.txt | sort) > ${base}.samples_include.txt
        echo -e "`date`\tdone"

        # add variants from panel comparison to exclusion
        cat ${exclude_variants_panel} >> ${base}.variants_exclude.txt

    >>>

    runtime {
        docker: "eu.gcr.io/finngen-refinery-dev/bioinformatics:0.6"
        memory: "10 GB"
        cpu: 4
        disks: "local-disk 200 HDD"
        preemptible: 2
    }

    output {
        File variants_exclude = base + ".variants_exclude.txt"
        File samples_exclude = base + ".samples_exclude.txt"
        File samples_include = base + ".samples_include.txt"
        File heterozygosity = base + ".heterozygosity.txt"
        File variant_missingness = base + ".vmiss"
        File sample_missingness = base + ".smiss"
        File sample_missingness_raw = base + ".raw.smiss"
        File sexcheck = base + ".sexcheck"
        File genome_raw = base + ".genome_raw.gz"
        File genome = base + ".genome.gz"
        File pihat_n_raw = base + ".pihat_n_raw.txt"
        File pihat_n = base + ".pihat_n.txt"
        File hardy = base + ".hardy"
    }
}

task duplicates {

    Array[File] samples_include
    Array[File] samples_exclude
    Array[File] sample_missingness
    String dollar = "$"

    command <<<

        set -euxo pipefail

        for file in ${sep=" " samples_include}; do mv $file .; done
        for file in ${sep=" " samples_exclude}; do mv $file .; done
        for file in ${sep=" " sample_missingness}; do mv $file .; done

        # concat inclusion and missingness files including batch id
        awk 'BEGIN{OFS="\t"} {sub(".samples_include.txt", "", FILENAME); print FILENAME,$0}' *.samples_include.txt > inc
        awk 'BEGIN{OFS="\t"} FNR==NR&&NR==1{print "batch",$0} FNR>1{sub(".smiss", "", FILENAME); print FILENAME,$0}' *.smiss > smiss

        # get missingness for included samples
        awk 'BEGIN{OFS="\t"} NR==FNR{a[$0]=1} NR>FNR && (FNR==1 || $1"\t"$3 in a)' inc smiss > inc_smiss
        if [[ `tail -n+2 inc_smiss | wc -l | awk '{print $1}'` != `wc -l inc | awk '{print $1}'` ]]; then
            >2& echo "samples do not match between sample inclusion files and smiss files"
            exit 1
        fi

        # get duplicates (2 or more _dup or same id) to remove keeping the sample with most variants called
        awk 'NR>1 { split($3, d, "_dup"); if ($5-$4>n[d[1]]) { n[d[1]]=$5-$4; i[d[1]]=$1"\t"$3 }} END {for (s in i) print i[s]}' inc_smiss > dupremoved
        comm -23 <(sort inc) <(sort dupremoved) > excdup

        # remove duplicates from each batch and create joint files
        for file in ${sep=" " samples_include}; do
            bname=`basename $file`
            batch=${dollar}{bname/".samples_include.txt"/}
            comm -12 <(sort $batch.samples_include.txt) <(awk -vbatch=$batch '$1==batch{print $2}' dupremoved | sort) > include && mv include $batch.samples_include.txt
            awk -v batch=$batch '$1==batch {OFS="\t"; print $2,"duplicate","NA"}' excdup >> $batch.samples_exclude.txt
            awk -v batch=$batch '{OFS="\t"; print batch,$0}' $batch.samples_include.txt >> all_batches.samples_include.txt
            awk -v batch=$batch '{OFS="\t"; print batch,$0}' $batch.samples_exclude.txt >> all_batches.samples_exclude.txt
        done

    >>>

    output {
        File allbatches_samples_include = "all_batches.samples_include.txt"
        File allbatches_samples_exclude = "all_batches.samples_exclude.txt"
    }

    runtime {
        docker: "eu.gcr.io/finngen-refinery-dev/bioinformatics:0.6"
        memory: "1 GB"
        cpu: 1
        disks: "local-disk 100 HDD"
        preemptible: 2
    }
}

task plots {

    String name
    Boolean joint
    String JOINT = if joint then "TRUE" else "FALSE"
    Array[File] sexcheck
    Array[File] variant_missingness
    Array[File] sample_missingness
    Array[File] heterozygosity
    Array[File] pihat_n
    Array[File] pihat_n_raw
    File exclude_samples
    Array[File] exclude_variants
    Array[Float] f
    Float variant_missing
    Float sample_missing
    Float het_sd
    Int pi_hat_min_n

    command <<<

        set -euxo pipefail

        # concat variants if not creating joint plots
        if ! ${joint}; then
            head -1 ${variant_missingness[0]} > vmiss; for file in ${sep=" " variant_missingness}; do tail -n+2 $file >> vmiss; done
        fi
        # concat files
        head -1 ${sexcheck[0]} > sexcheck; for file in ${sep=" " sexcheck}; do tail -n+2 $file >> sexcheck; done
        head -1 ${sample_missingness[0]} > smiss; for file in ${sep=" " sample_missingness}; do tail -n+2 $file >> smiss; done
        head -1 ${heterozygosity[0]} > heterozygosity; for file in ${sep=" " heterozygosity}; do tail -n+2 $file >> heterozygosity; done
        cat ${sep=" " pihat_n} > pi_hat
        cat ${sep=" " pihat_n_raw} > pi_hat_raw
        # get excluded samples only for this batch if not creating joint plots
        if ! ${joint}; then
            awk '$1=="${name}"' ${exclude_samples} > excluded_samples
        else
            mv ${exclude_samples} excluded_samples
        fi
        for file in ${sep=" " exclude_variants}; do mv $file .; done
        awk '{OFS="\t"; sub(".variants_exclude.txt", "", FILENAME); print FILENAME,$0}' *variants_exclude.txt | grep -Ev "not_in_panel|af_panel" > excluded_variants

        # create plots
        Rscript - <<EOF

        require(data.table)
        require(tidyverse)
        require(ggplot2)
        options(bitmapType='cairo')

        name <- "${name}"
        f1 <- ${f[0]}
        f2 <- ${f[1]}
        miss_thres <- c(${sample_missing}, ${variant_missing})
        miss_lim <- list() # sample/variant missingness x limits
        miss_lim[[1]] <- c(0.01, 0.05)
        miss_lim[[2]] <- c(0, 0.02)
        het_sd <- ${het_sd}
        pihat_n_lim <- ${pi_hat_min_n}
        pihat_x_lim <- 30

        file <- "sexcheck"
        data <- fread(file)
        png(paste0(name, "_", file, ".png"), width=1000, height=1000, units="px")
        p <- ggplot(data, aes(F)) +
          geom_density() +
          xlim(-6,1) +
          geom_vline(xintercept=f2, colour="black", linetype = "longdash") +
          geom_vline(xintercept=f1, colour="black", linetype = "longdash") +
          ggtitle(paste(name, file, f1, f2)) +
          theme_minimal(base_size=24) +
          scale_colour_manual(values=c("red", "blue"))
          print(p)
          dev.off()

        if (${JOINT}) {
            files <- c("smiss")
        } else {
            files <- c("smiss", "vmiss")
        }
        for (i in 1:length(files)) {
          file=files[i]
          data <- fread(file)
          for (j in 1:2) {
            png(paste0(name, "_", file, "_", j, "_", miss_lim[[j]][1], "_", miss_lim[[j]][2], ".png"), width=1000, height=1000, units="px")
            p <- ggplot(data, aes(F_MISS)) +
              geom_density() +
              xlim(miss_lim[[j]][1],miss_lim[[j]][2]) +
              ggtitle(paste(name, file, miss_lim[[j]][1], miss_lim[[j]][2])) +
              geom_vline(xintercept=miss_thres[i], colour="black", linetype = "longdash") +
              theme_minimal(base_size=24) +
              scale_colour_manual(values=c("red", "blue"))
            print(p)
            dev.off()
          }
        }

        file <- "heterozygosity"
        data <- fread(file)
        png(paste0(name, "_", file, ".png"), width=1000, height=1000, units="px")
        p <- ggplot(data, aes(het_rate)) +
          geom_density() +
          geom_vline(xintercept=mean(data[["het_rate"]])+het_sd*sd(data[["het_rate"]]), colour="black", linetype = "longdash") +
          ggtitle(paste(name, file)) +
          theme_minimal(base_size=24) +
          scale_colour_manual(values=c("red", "blue"))
        print(p)
        dev.off()

        files <- c("pi_hat_raw", "pi_hat")
        for (i in 1:length(files)) {
          file=files[i]
          data <- fread(file, header = F) %>% rename(IID=V1, pihat_n=V2)
          for (j in 1:2) {
            png(paste0(name, "_", file, "_",j,".png"), width=1000, height=1000, units="px")
            p <- ggplot(data, aes(pihat_n)) +
              geom_histogram(binwidth=1) +
              geom_vline(xintercept=pihat_n_lim, colour="black", linetype = "longdash") +
              ggtitle(paste(name, file)) +
              theme_minimal(base_size=24) +
              scale_colour_manual(values=c("red", "blue"))
            if (j==2) {
              p <- p + xlim(0, pihat_x_lim)
            }
            print(p)
            dev.off()
          }
        }

        files <- c("excluded_samples", "excluded_variants")
        for (i in 1:length(files)) {
          file=files[i]
          data <- fread(file, header=F) %>%
            rename(batch=V1,id=V2,reason=V3,value=V4) %>%
            group_by(batch, reason) %>% count()
          png(paste0(name, "_", file, ".png"), width=1000, height=1000, units="px")
          p <- ggplot(data, aes(x=batch, y=n, fill=reason)) +
            geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
            theme_minimal(base_size=16) +
            xlab("") + ylab("") +
            ggtitle(paste(name, file)) +
            coord_flip()
          print(p)
          dev.off()
        }
        EOF
    >>>

    output {
        Array[File] png = glob("*.png")
    }

    runtime {
        docker: "eu.gcr.io/finngen-refinery-dev/bioinformatics:0.6"
        memory: "3 GB"
        cpu: 1
        disks: "local-disk 100 HDD"
        preemptible: 2
    }
}
