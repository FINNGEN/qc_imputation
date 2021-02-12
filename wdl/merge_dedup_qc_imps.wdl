### after separate runs of e.g. legacy and affy and with list of same individual ids
## we need to pick single sample to keep across those same individuals
import "dedup.wdl" as subset_samples

workflow submerge {

    ## use scripts/generate_dedup_merge_conf to generate input from multiple cromwell qc_import runs for this
    String remove_samples
    ## batches by chr order

    String per_chr_batch_files

    Array[Array[String]] all_batch_per_chrom_shards=read_tsv(per_chr_batch_files)

    ## per chr output names of merged files
    String basename
    Array[Int] chrs

    scatter (chr in range(length(all_batch_per_chrom_shards))) {
        call subset_samples.dedup {
            input: batch_files=all_batch_per_chrom_shards[chr], remove_samples=remove_samples
        }

        call merge {
            input: vcfs=dedup.subvcfs, outfile=basename+"_chr"+chrs[chr]+".vcf.gz"
        }
    }

}

task merge {
    Array[File] vcfs
    String outfile
    Int cpus = 32

    command <<<
        set -euxo pipefail
        echo "Execute vcf-fusion&paste&bgzip command"
        time vcf-fusion ${sep=" " vcfs}| bgzip -@${cpus} > ${outfile}
        tabix -p vcf ${outfile}
    >>>

    output {
        File out = outfile
        File idx = outfile + ".tbi"
    }

    runtime {
        docker: "gcr.io/finngen-refinery-dev/qc_imputation:test_paste_0.1.0"
        memory: "12 GB"
        cpu: cpus
        disks: "local-disk 600 HDD"
        preemptible: 0
    }



}
