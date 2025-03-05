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

    Int disk_size = ceil(size(vcf, "GB")) * 2 + 5
    
    command <<<

        set -x pipefail
        ## something silently "fails" in the next command although result is correct... had to drop -e from pipefail for now.
        zcat ~{vcf} | grep -m 1 '^#CHROM' | tr '\t' '\n' | tail -n+10 > included_samples
        comm -12 <(sort -b ~{removals} ) <(sort -b included_samples ) > included_to_be_removed
        n_rem=$(wc -l included_to_be_removed|awk '{ print $1}')
        if [[ $n_rem -eq 0 ]];
        then
            echo "No samples to be removed from this batch"
            cmd="zcat -f"
        else
            echo "Subsetting "$n_rem" samples"
            cmd="bcftools view -S ^included_to_be_removed --force-samples -Ov"
        fi
        $cmd ~{vcf} | java -jar /tools/tagEditing_v1.1.jar AF: AN: AC: AC_Hom: AC_Het: NS: INFO: | bgzip -@ 2 > ~{bn}"_subset.vcf.gz"
        tabix -p vcf ~{bn}"_subset.vcf.gz"
    >>>

    output {
        File exclude_vcf=bn+"_subset.vcf.gz"
        File exclude_vcf_idx=bn+"_subset.vcf.gz.tbi"
        File removed_samples="included_to_be_removed"
    }

    runtime {
        docker: docker
        memory: "7 GB"
        cpu: 4
        disks: "local-disk " + disk_size + " HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 2
    }

}