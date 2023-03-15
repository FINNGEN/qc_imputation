# qc_imputation

Genotype QC and imputation pipeline

Takes a list of VCF files and input parameters and performs variant and sample wise QC for both chip QC and imputation purposes

This pipeline has been developed based on characteristics of FinnGen chip genotype data

## QC steps ([wdl/qc_imputation.wdl](wdl/qc_imputation.wdl))

1. `vcf_to_bed`: Convert each VCF file to PLINK bed

- done per genotyping batch
- .vcf or .vcf.gz accepted
- excludes samples in a given text file (one sample per row)
- possibility to filter samples based on prefix (in FinnGen we use prefix ^FG here to rid non-FinnGen samples up front)
- aligns to given human reference genome
- removes non-PASS variants and non-ATCG variants (CNVs)

2. `glm_vs_panel`: Run PLINK glm firth-fallback against imputation panel to detect frequency deviations

- done per genotyping batch
- merges each bed with the panel on shared variants and runs glm dataset vs. panel
- covariates for glm: sex, PCs

3. `missingness`: Compute variant missingness and list variants with high missingness

- done per genotyping batch
- lists variants with missingness above the given threshold for exclusion

4. `snpstats`: Compute qctool snpstats (HWE p, missingness etc.)

- done per chromosome across all genotyping batches
- HWE p-value for X chromosome is based on genetic females only

5. `gather_snpstats_joint_qc`: Combines snpstats across chromosomes and lists variants to be excluded from all batches based on stats

- done across all genotyping batches
- creates a list of variants to exclude from all batches
- excludes variants that were non-PASS (output of `vcf_to_bed`) in at least the given number of batches
- excludes variants with high missingness (output of `missingness`) in at least the given number of batches
- excludes variants with overall MAF below the given threshold (usually not used in practice i.e. the threshold is 0)
- excludes variants with overall HWE exact p-value less than the given threshold (for X chromosome HWE is computed with genetic females only) - exception: variants below a given frequency that have a deficiency in homozygotes escape this exclusion, e.g. rs201754030 in TSFM is severe in recessive state and would be excluded based on HWE without this exception
- excludes variants with overall missingness above the given threshold (usually not used in practice i.e. the threshold is 1)

6. `panel_comparison`: Creates batch-panel comparison plots and lists variants to exclude based on glm

- done per genotyping batch
- takes glm stats computed in `glm_vs_panel`
- lists variants with glm p-value less than the given threshold for exclusion
- allows variants in given regions pass the glm check (such regions should be defined for disease-enriched batches with known genetic signal for the disease)
- for imputation lists variants that are not in the panel or whose AF in the panel is less than the given threshold
- creates AF comparison plots batch vs. panel

7. `batch_qc`: Performs batch-wise QC

- done per genotyping batch
- creates a list of variants to exclude from each batch
- variants listed for exclusion by `vcf_to_bed`, `gather_snpstats_joint_qc` and `panel_comparison` are excluded up front before batch-wise QC
- excludes variants with HWE exact p-value less than the given threshold
- excludes variants with missingness above the given threshold (can be a different threshold for Y chromosome)
- excludes samples that fail PLINK sex check with the given f thresholds
- excludes samples that don't have the same sex as indicated by SSN
- excludes samples with missingness above the given threshold
- excludes samples with heterozygosity more than the given number of standard deviations above the mean (contamination)
- excludes samples with an excessive number of relatives based on two round of pi-hat calculation with the given thresholds

8. `duplicates`: Lists samples to include/exclude in each batch removing duplicates (same sample ID or \_dupX suffix)

- in a duplicate pair, the sample with more variants called is included in the data

9. `plots`: Creates QC plots and creates variant and sample exclusion summaries

- done both for each genotyping batch and across batches

10. `filter_batch`: Filters the PLINK files of each batch file excluding samples and variants based on the above QC

- done per genotyping batch
- done for chip qc purpose only, not for imputation
- excludes variants based on QC across batches (output of `gather_snpstats_joint_qc`)
- excludes variants based on batch-wise QC (output of `batch_qc`)
- excludes samples based on batch-wise QC with duplicates removed (output of `duplicates`)

11. `merge_batches`: Merges genotyping batches

- done per chromosome
- done for chip qc purpose only, not for imputation
- merges genotyping batch PLINK datasets (output of `filter_batch`) to a VCF

Steps 12-14 are done in the `imp_sub.wdl` subworkflow for imputation if the `run_imputation` input flag is set

12. `plink_to_vcf`: Converts the PLINK files of each batch to VCFs excluding samples and variants based on the above QC for imputation

- done per genotyping batch and per chrosomome
- excludes variants based on QC across batches (output of `gather_snpstats_joint_qc`)
- excludes variants based on batch-wise QC (output of `batch_qc`)
- excludes variants not in the imputation panel or low AF in the panel (output of `panel_comparison`)
- excludes samples based on batch-wise QC with duplicates removed (output of `duplicates`)

13. `phase_impute`: Phases (eagle) and imputes (beagle) each chromosome of each batch

- done per genotyping batch and per chrosomome
- uses the filtered genotype files (output of `plink_to_vcf`)
- outputs phased genotype data and imputed genotype data

14. `post_imputation`: Annotates imputed genotype files with INFO tags

- done per genotyping batch and per chrosomome
- adds AF, INFO, CHIP, AC_Het, AC_Hom, HWE tags to the VCF INFO field
- creates AF/INFO plots
- creates a tabix index

15. `post_subset_sub.subset_samples`: Removes SSN duplicates and denials across batches

- done across all genotyping batches per chromosome
- done separately for chip qc and imputation purposes
- removes remaining SSN duplicates and denials based on given lists
- biobank returns need to contain denials and SSN duplicates so this is done as the last thing

16. `merge_chip_data`: Merges chip data into one dataset

- done across all genotyping batches and chromosomes
- merges chip data after `post_subset_sub.subset_samples` to one VCF and PLINK dataset

## Configuration

Note that all input text files should be UTF-8 encoded

QC/imputation pipeline qc_imputation.wdl inputs

```
    "qc_imputation.docker" docker image to use in QC tasks
    "qc_imputation.imputation.docker" docker image to use in imputation tasks
    "qc_imputation.merge_in_chunks.docker": docker image to use in imputed data merging
    "qc_imputation.name": name of the run, e.g. r12_affy
    "qc_imputation.create_chip_dataset": true/false, whether to create a chip dataset across batches and chromosomes
    "qc_imputation.run_imputation": true/false, whether to run also imputation or only QC - best to first set this to false and check that the QC is fine before running imputation
    "qc_imputation.merge_imputed_batches": true/false, whether to merge imputed batches (usually not necessary as merging will be later done across all legacy and affy batches)
    "qc_imputation.imputation.force_impute_variants": chromosome-to-variant_list_location dictionary, per-chr lists of variants to force imputation of even if they otherwise pass QC (can be an empty file)
    "qc_imputation.snpstats.run_joint_qc": true/false, whether to run joint qc across all batches or not
    "qc_imputation.batch_qc.check_ssn_sex": true/false, whether to check sex against provided list of social security number based sex
    "qc_imputation.vcf_loc": location of the list of genotype data VCF files to run
    "qc_imputation.fam_loc": [OPTIONAL] location of genotype data FAM files corresponding to the VCF files - used for sex check if given
    "qc_imputation.ignore_panel_comparison_regions_loc": location of the list of files per batch that include genetic regions that are allowed to fail glm comparison to imputation panel (useful if there's disease enrichment in a batch and disease-associated regions would otherwise be QC'd out), files in this list can be empty files
    "qc_imputation.exclude_samples_loc": location of the list of files per batch that include samples to exclude from the run
    "qc_imputation.duplicate_samples": location of tab-delimited list of duplicate ids (same individual genotyped many times)
    "qc_imputation.exclude_denials": location of denial list
    "qc_imputation.include_regex": regex of samples to include in the run, e.g. "^FG" to only include FinnGen ids
    "qc_imputation.vcf_to_bed.ref_fasta": location of reference genome FASTA file
    "qc_imputation.chrs": list of chromosomes to run QC for (can include 24 (Y) and 26 (MT) that won't be imputed)
    "qc_imputation.imputation_chrs": list of chromosomes to impute - it's possible to e.g. run imputation of a specific chromosome after doing QC for all chromosomes
    "qc_imputation.genome_build": genome build version (38)
    "qc_imputation.panel_comparison.p": p-value threshold to use in excluding variants based on GWAS against panel, e.g. 5e-8
    "qc_imputation.panel_comparison.af_panel": variants with AF smaller than this threshold will not be given to imputation
    "qc_imputation.gather_panel_comparison.glm_panel_n_batches_prop": variants failing panel comparison in at least this proportion of batches will be excluded from all batches
    "qc_imputation.high_ld_regions": location of list of high-LD regions of the genome
    "qc_imputation.panel_bed": location of imputation panel PLINK bed file
    "qc_imputation.glm_vs_panel.pca_ld": PLINK LD parameters to use in PCA in comparison against imputation panel
    "qc_imputation.glm_vs_panel.pca_maf": minimum allele frequency to use in PCA in comparison against imputation panel
    "qc_imputation.missingness.variant_missing_n_batches": if there are at least `qc_imputation.gather_snpstats_joint_qc.missingness_n_batches_prop` times n_batches batches with missingness higher than this threshold, the variant will be excluded from all batches
    "qc_imputation.gather_snpstats_joint_qc.missingness_n_batches_prop": variants with missingness higher than `qc_imputation.missingness.variant_missing_n_batches` in at least this proportion of batches will be excluded from all batches
    "qc_imputation.gather_snpstats_joint_qc.af": minimum allele frequency in the genotype data - variants below this frequency will be excluded from imputation - this can be 0 as we're using the imputation panel to determine minimum frequency to use
    "qc_imputation.gather_snpstats_joint_qc.hw": HWE p-value threshold - variants below this threshold across all batches will be excluded
    "qc_imputation.gather_snpstats_joint_qc.escape_hwe_maf": variants with frequency smaller than this threshold and with homozygote deficiency escape HWE check across batches
    "qc_imputation.gather_snpstats_joint_qc.non_pass_n_batches_prop": variants that didn't pass Affy QC pipeline in at least this proportion of batches will be removed from all batches
    "qc_imputation.gather_snpstats_joint_qc.variant_missing_overall": maximum variant missingness computed across all batches - variants above this missingness will be excluded - this can be 1 because different batches may contain different variants which can lead to high overall missingness
    "qc_imputation.f": PLINK F thresholds for determining genotype sex - individuals between this range will be excluded
    "qc_imputation.batch_qc.ignore_sexcheck_samples": list of sample ids one per line that are allowed to pass sex check (e.g. XXY individuals), can be an empty file
    "qc_imputation.batch_qc.ssn_sex": location of list of social security number based sexes
    "qc_imputation.batch_qc.pi_hat": pi-hat threshold to use in detecting sample contamination
    "qc_imputation.batch_qc.variant_missing": maximum variant missingess in batch QC - variants above this missingness will be excluded from imputation of the batch
    "qc_imputation.batch_qc.variant_missing_y": maximum variant missingess for Y chromosome in batch QC - variants above this missingness will be excluded from imputation of the batch
    "qc_imputation.batch_qc.sample_missing": maximum sample missingess in batch QC - samples above this missingess will be excluded from imputation
    "qc_imputation.batch_qc.pihat_ld": PLINK LD parameters for pruning for pi-hat calculation
    "qc_imputation.batch_qc.hw": HWE p-value threshold for the batch - variants below this will be excluded from imputation of the batch
    "qc_imputation.batch_qc.het_sd": heterozygosity standard deviation threshold X - samples with heterozygosity lower than mean - X*sdev or higher than mean + X*sdev will be excluded
    "qc_imputation.batch_qc.pi_hat_min_n_excess": maximum excessive number of relatives based on pi-hat - in the first round of pi-hat calculation samples with more than this number of relatives will be excluded
    "qc_imputation.batch_qc.pi_hat_min_n": maximum number of relatives based on pi-hat - in the second round of pi-hat calculation samples with more than this number of relatives will be excluded
    "qc_imputation.ref_panel": chromosome-to-imputation_panel_vcf_location dictionary
    "qc_imputation.imputation.phase_impute.genetic_maps_eagle": chromosome-to-genetic_map_for_eagle dictionary
    "qc_imputation.imputation.phase_impute.genetic_maps_beagle": chromosome-to-genetic_map_for_beagle dictionary
    "qc_imputation.imputation.post_imputation.ref_panel_freq": location of allele frequencies in the imputation panel
    "qc_imputation.imputation.post_imputation.annot_hdr": location of header file to use in post-imputation annotation
    "qc_imputation.imputation.post_imputation.annot_tab": location of gzipped .tab file to use in post-imputation annotation
    "qc_imputation.imputation.post_imputation.annot_tab_index": location of tabix index of the above
    "qc_imputation.imputation.post_imputation.annot_col_incl": list of inclusion columns in post-imputation annotation
```

## Running

The QC/imputation pipeline can be run using Cromwell. [CromwellInteract](https://github.com/FINNGEN/CromwellInteract) command-line tool can be used

Example:

zip dependencies

```
cd wdl
zip sub.zip imp_sub.wdl post_subset_sub.wdl merge_chunks.wdl
cd ..
```

_NOTE: you need to have socks ssh tunnel running in port 5000 connected to cromwell machine_

Create tunnel if necessary:

```
CromwellInteract.py connect cromwell-machine-name
```

Then run pipeline, e.g. for legacy batches

```
CromwellInteract.py submit --wdl wdl/qc_imputation.wdl --inputs wdl/imputation.r7.legacy.json --deps wdl/sub.zip
```

## Copying output to targe buckets

### Copy imputed file results

These examples from R7 final imputation runs.

Files before duplicate removal created in `qc_imputation.imputation.vcfs` task:

```
cromwell_interact.py outfiles 7c9804b0-1850-4c76-8643-c628332db399 qc_imputation.imputation.vcfs | awk '{ print $1; print $1".tbi"}' | gsutil -m cp -I  gs://from_analysis_team/r7imputation/release_run/before_duplicate_removal_for_BB/
```

```
cromwell_interact.py outfiles a371ed62-cf5c-43ee-b755-8790eb3efe0c qc_imputation.imputation.vcfs | awk '{ print $1; print $1".tbi"}' | gsutil -m cp -I  gs://from_analysis_team/r7imputation/release_run/before_duplicate_removal_for_BB/
```

Files after merging in `qc_imputation.paste.out` output tag:

```
cromwell_interact.py outfiles a371ed62-cf5c-43ee-b755-8790eb3efe0c qc_imputation.paste.out | awk '{ print $1; print $1".tbi"}' | gsutil -m cp -I  gs://output_bucket
```

## Merge multiple runs and remove duplicates across runs

After successful run of different qc/imputation pipelines( e.g. separately for legacy and affy), we need to do duplicate removal across all batches in all runs.

qc/imputation pipeline creates sample qc summary stats and we pick the sample with the most variants after qc as the genotype ID to keep.

_NOTE: you need to have socks ssh tunnel running in port 5000 connected to cromwell machine in all below commands_.
Create tunnel if necessary:

```
cromwell_interact.py connect cromwell-machine-name
```

Create configuration file using:

```
scripts/generate_dedup_merge_conf.py hash1,hash2 duplicateids_file outprefix
```

Parameters:

1. comma separated list of cromwell job hashes (e.g. affy run, legacy run)
2. file with duplicate ids for a single individual. Each row should have tab separated list of alternate genotyping IDs.
3. outputprefix of created files

Example from R7 production runs:

```
scripts/generate_dedup_merge_conf.py "a371ed62-cf5c-43ee-b755-8790eb3efe0c,7c9804b0-1850-4c76-8643-c628332db399" finngen_R7_duplicate_list.txt finngen_R7
```

**Upload the created configuration files to bucket**

Files:

- {outprefix}\_final_sample_removals . Add any sample IDs to be removed here in case additional removals are needed (e.g. samples in that are not in exclusions list).
- {outprefix}\_all_batches_per_chr

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
cromwell_interact.py submit --wdl wdl/merge_dedup_qc_imps.wdl --inputs wdl/merge_dedup_qc_imps_remove_non_inclusion_r7.json --deps wdl/merge_dedup_qc_imps_deps.zip
```

**copy output**

Per batch duplicates removed
Final merged files
