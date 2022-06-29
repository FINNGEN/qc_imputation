

workflow subset_samples {

    Array[File] vcfs
    Array[File] vcf_idxs
    File sample_summaries

    File already_excluded_samples
    File? duplicate_samples
    File? exclude_denials

    String docker

    scatter( b in range(length(vcfs)) ) {
        ## add task that figures out which of same person with different IDs to remove and pass to next task.
        ## good to make all this optional for other imputation tasks, at least at some point.
        call subset { input: vcf=vcfs[b], vcf_idx=vcf_idxs[b], already_excluded_samples=already_excluded_samples,
            duplicate_samples=duplicate_samples, sample_summaries=sample_summaries,
            docker=docker, exclude_denials=exclude_denials
        }
    }

    output {
        Array[File] subset_vcfs=subset.exclude_vcf
        Array[File] subset_vcfs_idx=subset.exclude_vcf_idx
    }
}

task subset {
    File vcf
    File vcf_idx
    File sample_summaries

    ## this is used only to see if some duplicates are not duplicates after qc.
    File already_excluded_samples

    File? exclude_denials
    File? duplicate_samples
    String docker
    String bn=basename(sub(sub(sub(vcf,".vcf.gz",""),".vcf.bgz",""),".vcf",""))


    command <<<
        set -euxo pipefail

        excl_denials_file=${true='' false='' defined(exclude_denials)}${exclude_denials}
        dupl_sampl_file=${true='' false='' defined(duplicate_samples)}${duplicate_samples}

        touch removals

        cut -f 2 ${already_excluded_samples} | sort | uniq > removals

        # remove _dup suffixes from sample summaries as there are no suffixes in the duplicate list
        sed -i 's/_dup[0-9]\+//' ${sample_summaries}

        if [[ "$excl_denials_file" != "" ]]
        then
            cat $excl_denials_file >> removals
        fi

        if [[ "$dupl_sampl_file" != "" ]]
        then

            awk 'FNR==NR&&!match($0,"^#"){ rems[$1]=""; array[1]=""}
                NR>FNR&&!match($0,"^#"){
                    delete array
                    incl=0;
                    for (i=1;i<=NF;i++) {
                        if(! ($i in rems) ) {
                            #mark ids to be included in duplicate row.
                            array[incl++]=$i
                        }
                    }
                    if(incl>1){
                        line=array[0]
                        for(i=1;i<incl;i++) {
                            line=line"\t"array[i]
                        }
                        print(line)
                    }
            }'  removals $dupl_sampl_file > dups_remaining

            awk 'BEGIN {FS=OFS="\t"}
                NR==1{
                        print("reading header")
                        print(NF)
                        for(i=1;i<=NF;i++) { h[$i]=i; print("storing header " $i) };
                        if(!("IID" in h && "N_GENO_in_batch_qc" in h && "n_excluded_vars" in h && "n_vars_total" in h && "N_MISS_in_batch_qc" in h)) {
                            print "Required columns missing from sample summary. Needs IID,N_GENO, n_vars_total and n_excluded_vars.";
                            exit 1;
                        }
                }
                NR>1 && NR==FNR {
                    vars[$h["IID"]]=$h["n_vars_total"]-$h["n_excluded_vars"];
                    miss_vars[$h["IID"]]=$h["N_MISS_in_batch_qc"];
                }
                NR>FNR {
                        ## check all duplicate IDs which to remove
                        bestid=""
                        best_score=0

                        for(i=1;i<=NF;i++){
                            if(!( $i in vars )){
                                ## no need to remove as not part of sample
                                continue;
                            }
                            n_vars=vars[$i]-miss_vars[$i];

                            if(n_vars>best_score) {
                                if(bestid!="") {print bestid}
                                bestid=$i;
                                best_score=n_vars;
                            } else {
                                print $i;
                            }
                        }
                } ' ${sample_summaries} dups_remaining > remove_duplicates
            cat remove_duplicates >> removals
        fi

        n_remove=$(cat removals | wc -l)
        if [[ $n_remove -eq 0 ]]
        then
            echo "Subset samples file not specified and no duplicates to remove. Not subsetting."
            mv ${vcf} ${bn}_temp.vcf.gz
        else
            echo "Subsetting "$n_remove" samples"
            bcftools view -S ^removals --force-samples ${vcf} -Ov | bgzip > ${bn}_temp.vcf.gz
        fi

        # index and add info tags
        echo "`date`\ttags"
        bcftools +fill-tags ${bn}_temp.vcf.gz -Oz -o ${bn}_subset.vcf.gz -- -t AF,AN,AC,AC_Hom,AC_Het,NS
        echo "`date`\tindex vcf"
        tabix -p vcf ${bn}"_subset.vcf.gz"
        echo "`date`\tdone"
        
    >>>


    output {
        File exclude_vcf=bn+"_subset.vcf.gz"
        File exclude_vcf_idx=bn+"_subset.vcf.gz.tbi"
        File removed_samples="removals"
    }

    runtime {
        docker: "${docker}"
        cpu: 1
        memory: "7 GB"
        disks: "local-disk 200 HDD"
        preemptible: 2
        zones: "europe-west1-b europe-west1-c europe-west1-d"
    }
}
