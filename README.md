# qc_imputation
Genotype QC and imputation pipeline

## Configuration
    *TODO*
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
CromwellInteract connect cromwell-machine-name
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
- {outprefix}_final_sample_removals`
- {outprefix}_all_batches_per_chr`

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
CromwellInteract submit --wdl wdl/merge_dedub_qc_imps.json --inputs wd/merge_dedup_qc_imps_modified.json --deps wdl/merge_dedup_qc_imps_deps.zip
```

**copy output**

Per batch duplicates removed
Final merged files
