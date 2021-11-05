
json_glob=$1
outbucket=$2

legacy_plink=$(mktemp)
legacy_exclusion=$(mktemp)
legacy_fam=$(mktemp)
copy_legacy_files=$(mktemp)

gsutil cat $json_glob | \
    jq ".[\"preimputation_QC.IN_FILE_BED\"],.[\"vcf_conversion_v1_1.ADDITIONAL_SAMPLE_EXCLUSION_LIST\"] " | tr -d '"' | \
    tee -a $copy_legacy_files | while read f; do echo $outbucket$(basename $f); done | \
        awk -v pl=$legacy_plink -v excl=$legacy_exclusion '{ if(NR % 2==1) file=pl; else file=excl; print $0>>file }'

gsutil cat gs://r7_data/R7_pilot_runs/legacy/factory_run_configs/*.json | \
    jq ".[\"preimputation_QC.IN_FILE_BED\"],.[\"vcf_conversion_v1_1.ADDITIONAL_SAMPLE_EXCLUSION_LIST\"] " | tr -d '"' | \
    tee -a copy | while read f; do echo $outbucket$(basename $f); done | \
        awk -v pl=legacy_plink -v excl=legacy_exclusion '{ if(NR % 2==1) file=pl; else file=excl; print $0>>file }'

#cat $copy_legacy_files | awk ' BEGIN{ suf[".bim"]=0; suf[".fam"]=0} $1~/\.bed$/{ for ( s in suf) { f=$1; gsub(".bed",s,f); print f} } { print $1}' \
#  | gsutil -m cp -I $outbucket

gsutil cp $legacy_plink $outbucket"legacy_plink"
gsutil cp $legacy_exclusion $outbucket"legacy_exclusion"

cat $legacy_plink | sed 's/.bed/.fam/g' > $legacy_fam
gsutil cp $legacy_fam $outbucket/legacy_fam

gsutil cp $outbucket"legacy_exclusion" |while read f; do b=$(basename $f); gsutil cat $f | \
    awk  'BEGIN{OFS="\t"}{ print $1,"list1"}' >> temppi; gsutil cp temppi "gs://r7_data/R7_pilot_runs/legacy/batches/"$b; rm temppi; done

rm $legacy_plink $legacy_exclusion $copy_legacy_files

##4084f69a-4d61-4b59-a805-7fa65ad0b2bc
