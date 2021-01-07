

workflow subset_samples {

    Array[File] vcfs
    File? exclude_samples

    scatter( b in range(vcfs) ) {
        ## add task that figures out which of same person with different IDs to remove and pass to next task.
        ## good to make all this optional for other imputation tasks, at least at some point.
        call subset { input: vcf=vcfs[b], exclude_samples=exclude_samples}
    }

    output {
        Array[File] subset_vcfs= subset.exclude_vcf
    }
}

task subset {
    File vcf
    File? exclude_samples

    command <<<

        excl_sampl_file=${true='' false='' defined(exclude_samples)}${exclude_samples}
        bn=`basename ${vcf} | sed -e 's|\.vcf\.gz$||g' -e 's|\.vcf\.bgz$||g' -e 's|\.vcf$||g'`

        if [[ "$excl_sampl_file" == "" ]]
        then
            echo "Subset samples file not specified. Not subsetting."
            mv ${vcf} $bn".vcf.gz"
        else
            echo "Subsetting samples with file ${exclude_samples}"
            bcftools view -S ^${exclude_samples} ${vcf} -Oz $bn".vcf.gz"
        fi
        echo $bn".vcf.gz" > off
    >>>

    output {
        File exclude_vcf=read_string("off")
    }

    runtime {
        docker: "eu.gcr.io/finngen-refinery-dev/bioinformatics:0.6"
        memory: "7 GB"
        cpu: 1
        disks: "local-disk 200 HDD"
        preemptible: 2
    }
}
