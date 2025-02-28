version 1.0
### after separate runs of e.g. legacy and affy and with list of same individual ids
## we need to pick single sample to keep across those same individuals
import "dedup2.wdl" as subset_samples
import "merge_chunks_X.wdl" as merge_sub

workflow submerge {

    input {
        ## use scripts/generate_dedup_merge_conf to generate input from multiple cromwell qc_import runs for this
        File remove_samples
        File chr_vcfs_loc

        Int chunk_window = 5000000
        String basename
        String chr = "23"
        String docker
    }

    Array[File] chr_vcfs = read_lines(chr_vcfs_loc)

    call subset_samples.dedup {
        input:
            batch_files = chr_vcfs,
            remove_samples = remove_samples
    }

    call merge_sub.merge_in_chunks {
        input:
            vcfs = dedup.subvcfs,
            outfile = basename + "_chr" + chr + ".vcf.gz",
            chunksize = chunk_window,
            docker = docker
    }

    output {
        File out = merge_in_chunks.out
        File idx = merge_in_chunks.idx
    }

}