### after separate runs of e.g. legacy and affy and with list of same individual ids
## we need to pick single sample to keep across those same individuals
import "dedup.wdl" as subset_samples
import "merge_chunks.wdl" as merge_sub

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

        call merge_sub.merge_in_chunks {
            input: batches=dedup.subvcfs, outfile=basename+"_chr"+chrs[chr]+".vcf.gz",chunksize=20000,docker="gcr.io/finngen-refinery-dev/qc_imputation:test_paste_0.1.0"
        }
    }

}

