version 1.0

workflow dedup {

    input {
        String remove_samples
        ## batches by chr order
        Array[String] batch_files
        String docker
    }

    scatter (f in batch_files) {
        call subset {
            input: vcf=f, removals=remove_samples, docker=docker
        }
    }

    output {
        Array[File] subvcfs = subset.exclude_vcf
        Array[File] subsvcfs_idx = subset.exclude_vcf_idx
        Array[File] removed_samples = subset.removed_samples
    }

}

task subset {

    input {
        File vcf
        String bn = basename(sub(sub(sub(vcf,".vcf.gz",""),".vcf.bgz",""),".vcf",""))
        File tbi=vcf+".tbi"
        File removals
        String docker
    }

    Int disk_size = ceil(size(vcf, "GB")) * 3 + 5
    
    command <<<

        set -x
        
        tabix -H ~{vcf} | tail -1 | tr '\t' '\n' | tail -n+10 > included_samples
        comm -12 <(sort -b ~{removals} ) <(sort -b included_samples ) > included_to_be_removed
        
        bcftools view -S ^included_to_be_removed --force-samples -Ou ~{vcf} | bcftools +fill-tags -Ou -- -t NS,AC,AF,AN,AC_Het,AC_Hom | bcftools +impute-info -Ou | bcftools view -Oz --threads 5 -Wtbi -o ~{bn}"_subset.vcf.gz" 

    >>>

    output {
        File exclude_vcf=bn+"_subset.vcf.gz"
        File exclude_vcf_idx=bn+"_subset.vcf.gz.tbi"
        File removed_samples="included_to_be_removed"
    }

    runtime {
        docker: docker
        memory: "32 GB"
        cpu: 8
        disks: "local-disk " + disk_size + " HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 0
    }

}