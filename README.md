# qc_imputation
Genotype QC and imputation pipeline

## Configuration

Note that all input text files should be UTF-8 encoded

QC/imputation pipeline imputation.wdl inputs

```
    "qc_imputation.docker" docker image to use
    "qc_imputation.name": name of the run, e.g. r7_legacy
    "qc_imputation.run_imputation": true/false, whether to run also imputation or only qc
    "qc_imputation.chr_qc.joint_qc.run_joint_qc": true/false, whether to run joint qc across all batches or not
    "qc_imputation.chr_qc.compare_panel.compare": true/false, whether to run comparison against panel
    "qc_imputation.batch_qc.check_ssn_sex": true/false, whether to check sex against provided list of social security number based sex
    "qc_imputation.vcf_loc": location of genotype data VCF files to run
    "qc_imputation.fam_loc": location of genotype data FAM files to run
    "qc_imputation.exclude_samples_loc": location of samples to exclude from the run
    "qc_imputation.duplicate_samples": location of tab-delimited list of duplicate ids (same individual genotyped many times)
    "qc_imputation.include_regex": regex of samples to include in the run, e.g. "^FG" to only include FinnGen ids
    "qc_imputation.chrs": list of chromosomes to run
    "qc_imputation.genome_build": genome build version (38)
    "qc_imputation.panel_comparison.p": p-value threshold to use in excluding variants based on GWAS against panel, e.g. 5e-8
    "qc_imputation.high_ld_regions": location of list of high-LD regions of the genome
    "qc_imputation.chr_qc.vcf_to_bed_chr.ref_fasta": location of reference genome FASTA file
    "qc_imputation.chr_qc.compare_panel.panel_freq": location of allele frequencies in the imputation panel
    "qc_imputation.chr_qc.compare_panel.pca_ld": PLINK LD parameters to use in PCA in comparison against imputation panel
    "qc_imputation.chr_qc.compare_panel.pca_maf": minimum allele frequency to use in PCA in comparison against imputation panel
    "qc_imputation.chr_qc.compare_panel.af_panel": minimum allele frequency in imputation panel - variants below this frequency will be excluded from imputation
    "qc_imputation.chr_qc.compare_panel.af_diff": maximum absolute frequency between data and imputation panel - variants above this difference will be excluded from imputation
    "qc_imputation.chr_qc.compare_panel.af_fc": log of maximum allele frequency fold change against imputation panel - variants above this fold change will be excluded from imputation
    "qc_imputation.chr_qc.joint_qc.af": minimum allele frequency in the genotype data - variants below this frequency will be excluded from imputation - this can be 0 as we're using the imputation panel to determine minimum frequency to use
    "qc_imputation.chr_qc.joint_qc.hw": HWE p-value threshold - variants below this threshold across all batches will be excluded
    "qc_imputation.chr_qc.joint_qc.variant_missing_overall": maximum variant missingness across all batches - variants above this missingness will be excluded - this can be 1 because different batches may contain different variants which can lead to high overall missingness
    "qc_imputation.chr_qc.joint_qc.variant_missing_single_batch": maximum variant missingness in a single batch - variants above this missingness in any batch will be excluded from all batches
    "qc_imputation.f": PLINK F thresholds for determining genotype sex - individuals between this range will be excluded
    "qc_imputation.batch_qc.ssn_sex": location of list of social security number based sexes
    "qc_imputation.batch_qc.pi_hat": pi-hat threshold to use in detecting sample contamination
    "qc_imputation.batch_qc.variant_missing": maximum variant missingess in batch QC - variants above this missingess will be excluded from imputation of the batch
    "qc_imputation.batch_qc.sample_missing": maximum sample missingess in batch QC - samples above this missingess will be excluded from imputation
    "qc_imputation.batch_qc.pihat_ld": PLINK LD parameters for pruning for pi-hat calculation
    "qc_imputation.batch_qc.hw": HWE p-value threshold for the batch - variants below this will be excluded from imputation of the batch
    "qc_imputation.batch_qc.het_sd": heterozygosity standard deviation threshold X - samples with heterozygosity lower than mean - X*sdev or higher than mean + X*sdev will be excluded
    "qc_imputation.batch_qc.pi_hat_min_n_excess": maximum excessive number of relatives based on pi-hat - in the first round of pi-hat calculation samples with more than this number of relatives will be excluded
    "qc_imputation.batch_qc.pi_hat_min_n": maximum number of relatives based on pi-hat - in the second round of pi-hat calculation samples with more than this number of relatives will be excluded
    "qc_imputation.ref_panel": chromosome-to-imputation_panel_vcf_location_with_SNPID dictionary 
    "qc_imputation.imputation.phase_impute.ref_panel": chromosome-to-imputation_panel_vcf_location dictionary
    "qc_imputation.imputation.phase_impute.genetic_maps_eagle": chromosome-to-genetic_map_for_eagle dictionary
    "qc_imputation.imputation.phase_impute.genetic_maps_beagle": chromosome-to-genetic_map_for_beagle dictionary
    "qc_imputation.imputation.post_imputation.ref_panel_freq": location of allele frequencies in the imputation panel
    "qc_imputation.imputation.post_imputation.annot_hdr": location of header file to use in post-imputation annotation
    "qc_imputation.imputation.post_imputation.annot_tab": location of gzipped .tab file to use in post-imputation annotation
    "qc_imputation.imputation.post_imputation.annot_tab_index": location of tabix index of the above
    "qc_imputation.imputation.post_imputation.annot_col_incl": list of inclusion columns in post-imputation annotation
```

## Running
    *TODO*

## Copying output to targe buckets

### Copy imputed file results

These examples from R7 final imputation runs.

Files before duplicate removal created in `qc_imputation.imputation.vcfs` task:
```
CromwellInteract.py outfiles 7c9804b0-1850-4c76-8643-c628332db399 qc_imputation.imputation.vcfs | awk '{ print $1; print $1".tbi"}' | gsutil -m cp -I  gs://from_analysis_team/r7imputation/release_run/before_duplicate_removal_for_BB/
```

```
CromwellInteract.py outfiles a371ed62-cf5c-43ee-b755-8790eb3efe0c qc_imputation.imputation.vcfs | awk '{ print $1; print $1".tbi"}' | gsutil -m cp -I  gs://from_analysis_team/r7imputation/release_run/before_duplicate_removal_for_BB/
```

Files after merging in `qc_imputation.paste.out` output tag:
```
CromwellInteract.py outfiles a371ed62-cf5c-43ee-b755-8790eb3efe0c qc_imputation.paste.out | awk '{ print $1; print $1".tbi"}' | gsutil -m cp -I  gs://output_bucket
```

## Merge multiple runs and remove duplicates across runs

After successful run of different qc/imputation pipelines( e.g. separately for legacy and affy), we need to do duplicate removal across all batches in all runs.

qc/imputation pipeline creates sample qc summary stats and pick sample with the most variants after qc as the genotype ID to keep.

*NOTE: you need to have socks ssh tunnel running in port 5000 connected to cromwell machine in all below commands*.
Create tunnel if necessary:
```
CromwellInteract.py connect cromwell-machine-name
```


Create configuration file using:
```
scripts/generate_dedup_merge_conf.py hash1,hash2 duplicateids_file outprefix
```

Parameters:
1. comma separated list of cromwell job hashes
2. file with duplicate ids for a single individual. Each row should have tab separated list of alternate genotyping IDs.
3. outputprefix of created files

Example from R7 production runs:
 ```
 scripts/generate_dedup_merge_conf.py "a371ed62-cf5c-43ee-b755-8790eb3efe0c,7c9804b0-1850-4c76-8643-c628332db399" finngen_R7_duplicate_list.txt finngen_R7
 ```

**Upload the created configuration files to bucket**

Files:
- {outprefix}_final_sample_removals . Add any sample IDs to be removed here in case additional removals are needed (e.g. samples in that are not in exclusions list).
- {outprefix}_all_batches_per_chr



In our example:
- finngen_R7_final_sample_removals
- finngen_R7_all_batches_per_chr.

Modify json in `wdl/merge_dedub_qc_imps.json` to point to those those config files.

**Submit duplicate removal and final merge by chromosome**

zip dependencies
```
zip -j wdl/merge_dedup_qc_imps_deps.zip wdl/dedup.wdl
```

**Run pipeline:**
```
CromwellInteract.py submit --wdl wdl/merge_dedup_qc_imps.wdl --inputs wdl/merge_dedup_qc_imps_remove_non_inclusion_r7.json --deps wdl/merge_dedup_qc_imps_deps.zip
```

**copy output**

Per batch duplicates removed
Final merged files


## Production runs
**Release 7 runs **

Data was ran in cromwell-fg-1
- Legacy data batches hash 7c9804b0-1850-4c76-8643-c628332db399  
- Affy data a371ed62-cf5c-43ee-b755-8790eb3efe0c
- Deduplication across batches f188e4b7-f004-483b-b250-319b087766ac
- Removal of final non-inclusions samples and deduplication combined 8934d7cc-97d1-47fe-8893-c8799943d2ff

##
List of final non included samples in above run
comm -23 <(gsutil cat gs://thl-incoming-data/from-data-team/R7/fgfactory_pass_samples_R7.txt | awk 'BEGIN{ FS=":";} NR>1{ print $6}' | sort -b) <(gsutil cat gs://r7_data/R7_lists_from_Tero/finngen_R7_finngenid_inclusion_list.txt | sort -b)  > r7_samples_not_in_inclusion
