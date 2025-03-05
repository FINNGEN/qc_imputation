import "imp_sub.wdl" as imp_sub
import "post_subset_sub.wdl" as post_subset_sub
import "merge_chunks.wdl" as  merge_sub

workflow impute {

    String docker
    String name
    Array[Int] imputation_chrs
    File bed_loc
    Array[String] beds = read_lines(bed_loc)
    File batch_qc_exclude_variants_loc
    Array[String] batch_qc_exclude_variants = read_lines(batch_qc_exclude_variants_loc)
    File panel_exclude_variants_loc
    Array[String] panel_exclude_variants = read_lines(panel_exclude_variants_loc)
    File allbatches_samples_exclude
    File exclude_denials
    File duplicate_samples
    Map[Int, String] ref_panel
    File sample_summaries
    Boolean merge_imputed_batches

    # run imputation per chromosome
    scatter (i in range(length(imputation_chrs))) {
        call imp_sub.imputation as imputation {
            input: chr=imputation_chrs[i], beds=beds,
            batch_qc_exclude_variants=batch_qc_exclude_variants,
            panel_exclude_variants=panel_exclude_variants,
            exclude_samples=allbatches_samples_exclude,
            ref_panel=ref_panel,
            add_batch_suffix=false
        }
        # exclude duplicates (same sample with two different ids) and denials
        call post_subset_sub.subset_samples as subset_samples {
            input: vcfs=imputation.vcfs, vcf_idxs=imputation.vcf_idxs,
            already_excluded_samples=allbatches_samples_exclude,
            exclude_denials=exclude_denials, duplicate_samples=duplicate_samples,
            sample_summaries=sample_summaries,
            add_batch_suffix=true,
            docker=docker
        }
    }
    # paste imputed batches per chr
    if ( length(beds)>1 && merge_imputed_batches ) {
        scatter (i in range(length(imputation_chrs))) {
            call merge_sub.merge_in_chunks {
                  input: vcfs=subset_samples.subset_vcfs[i],
                  outfile=name+"_all_chr"+imputation_chrs[i]+".vcf.gz",
                  chunksize=20000
            }
        }
    }
}