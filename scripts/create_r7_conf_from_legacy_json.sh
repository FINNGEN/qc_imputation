

### get samples to exclude from
# get plink data from preimputation_QC.IN_FILE_BED
#* vcf_conversion_v1_1.ADDITIONAL_SAMPLE_EXCLUSION_LIST
#
# there is this inclusion list everywhere listing all FG ids... is this really necessary?
#gs://fgf-ref-files/R7_lists_from_Tero/finngen_R7_finngenid_inclusion_list.tx
# DUplicates per release are in paired file: gs://fgf-ref-files/R7_lists_from_Tero/finngen_R7_duplicate_list.txt
# MAybe we should add duplicate chooser at relevant point. We can and should impute everything anyway and just get phenos for one or the
# other.


scripts/get_legacy_data.sh "gs://r7_data/R7_pilot_runs/legacy/factory_run_configs/*.json" "gs://r7_data/R7_pilot_runs/legacy/batches/"

gsutil -m cp gs://fg-cromwell_fresh/convert_plink_to_fg_vcf/f349308c-1df7-497d-b066-f5badea25257/call-convert/**/*.vcf.gz \
    gs://r7_data/R7_pilot_runs/legacy/batches/converted_vcf/

##
## After plink conversion, generate fam and vcf file lists for imputation pipe in the same order.
#gs://r7_data/R7_pilot_runs/legacy/batches/converted_vcf/

fam_file=$(mktemp)
vcf_file=$(mktemp)
gsutil cat gs://r7_data/R7_pilot_runs/legacy/batches/legacy_plink | \
    while read f; do bf=$(basename $f);  \
                    echo "gs://r7_data/R7_pilot_runs/legacy/batches/converted_vcf/"${bf%.bed}".vcf.gz" >> $vcf_file; \
                    echo ${f%.bed}".fam" >> $fam_file; \
                done
gsutil cp $vcf_file gs://r7_data/R7_pilot_runs/legacy/batches/legacy_vcf
gsutil cp $fam_file gs://r7_data/R7_pilot_runs/legacy/batches/legacy_fam

rm $vcf_file
rm $fam_file
