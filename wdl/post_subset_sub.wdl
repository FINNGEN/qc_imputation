

workflow subset_samples {

    Array[File] vcfs
    Array[File] vcf_idxs
    File? exclude_samples

    scatter( b in range(length(vcfs)) ) {
        ## add task that figures out which of same person with different IDs to remove and pass to next task.
        ## good to make all this optional for other imputation tasks, at least at some point.
        call subset { input: vcf=vcfs[b], vcf_idx=vcf_idxs[b], exclude_samples=exclude_samples}
    }

    output {
        Array[File] subset_vcfs=subset.exclude_vcf
    }
}

task subset {
    File vcf
    File vcf_idx
    File? exclude_samples

    String bn=basename(sub(sub(sub(vcf,".vcf.gz",""),".vcf.bgz",""),".vcf",""))

    command <<<
        set -euxo pipefail

        excl_sampl_file=${true='' false='' defined(exclude_samples)}${exclude_samples}

        if [[ "$excl_sampl_file" == "" ]]
        then
            echo "Subset samples file not specified. Not subsetting."
            mv ${vcf} ${bn}"_subset.vcf.gz"
        else
            echo "Subsetting samples with file ${exclude_samples}"
            bcftools view -S ^${exclude_samples} --force-samples ${vcf} -Ov | bgzip > ${bn}"_subset.vcf.gz"
            tabix -p vcf ${bn}"_subset.vcf.gz"
        fi
    >>>

    output {
        File exclude_vcf=bn+"_subset.vcf.gz"
        File exclude_vcf_idx=bn+"_subset.vcf.gz.tbi"
    }

    runtime {
        docker: "eu.gcr.io/finngen-refinery-dev/bioinformatics:0.6"
        memory: "7 GB"
        cpu: 1
        disks: "local-disk 200 HDD"
        preemptible: 2
    }
}
