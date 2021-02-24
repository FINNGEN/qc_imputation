

workflow convert_plink_to_fg_vcf {

    File inputlist
    Array[String] plink_stems=read_lines(inputlist)

    scatter (p in plink_stems) {
        call convert {
            input: plinkstem=p
        }
    }
}

task convert {
    String plinkstem
    File bed=sub(plinkstem, ".bed","") +".bed"
    File bim=sub(plinkstem, ".bed","") +".bim"
    File fam=sub(plinkstem, ".bed","") +".fam"
    String base=basename(plinkstem,".bed")
    command <<<
        plink2 --bed ${bed} --bim ${bim} --fam ${fam} --export vcf-4.2 id-paste=iid --out temp
        cat temp.vcf | awk ' BEGIN {IFS="\t"; OFS="\t"} !/^#/{ $7="PASS";} { print $0}' | bgzip > ${base}.vcf.gz
        tabix -p vcf ${base}.vcf.gz
    >>>

    runtime {
        docker: "eu.gcr.io/finngen-refinery-dev/bioinformatics:0.6"
        memory: "10 GB"
        cpu: 1
        disks: "local-disk 200 HDD"
        preemptible: 2
    }

    output {
        File vcf=base + ".vcf.gz"
        File tbi=base + ".vcf.gz.tbi"
    }

}
