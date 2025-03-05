version 1.0

workflow merge_in_chunks {

    input {
        Array[File] vcfs
        String outfile
        Int chunksize
        String docker
    }

    scatter (i in range(length(vcfs))){
        call chunk {
            input:
                batch = vcfs[i],
                chunksize = chunksize,
                docker = docker
        }
    }

    call reorder_array{
        input:
            batch_chunk_array = chunk.globbed
    }

    Array[Array[File]] chunk_batch_array = transpose(reorder_array.batch_chunk)

    scatter (i in range( length( chunk_batch_array ) ) ){
        call merge {
            input:
                chunks = chunk_batch_array[i],
                number = i,
                docker = docker
        }
    }

    call concatenate_chunks {
        input:
            chunks = merge.out,
            docker = docker,
            filename = outfile
    }

    output {
        File out = concatenate_chunks.out
        File idx = concatenate_chunks.idx
    }
}

task reorder_array {

    input {
        Array[Array[String]] batch_chunk_array
        File tsv = write_tsv(batch_chunk_array)
    }

    command {
        set -euxo pipefail
        python3 <<EOF
        #read data
        with open("~{tsv}","r") as f:
            lines = f.readlines()
            chunks = [a.strip("\n").split("\t") for a in lines]
        #keep line ordering, as batches are on lines, but order them so that chunks are from 1,2,3..n
        with open("output","w") as OUT:
            for chunk in chunks:
                ch_tuples = [(a,int(a.split("_")[-1].split(".")[0]) ) for a in chunk]
                keyfunc = lambda x: x[1]
                sorted_chunks = [a[0] for a in sorted(ch_tuples,key=keyfunc)]
                outline = "\t".join(sorted_chunks)+"\n"
                OUT.write(outline)

        EOF
    }

    runtime {
        docker: "python:3.6-slim"
        memory: "2 GB"
        cpu: 1
        disks: "local-disk 10 HDD"
        preemptible: 2
        zones: "europe-west1-b europe-west1-c europe-west1-d"
    }

    output {
        Array[Array[String]] batch_chunk = read_tsv("output")
    }
}

task chunk {

    input {
        File batch
        Int chunksize
        String docker
        String prefix = basename(batch,".vcf.gz")
    }

    command <<<
        set -euxo pipefail
        cat << EOF > blocker.py
        """
        Separate a VCF file into smaller gzipped blocks. Each block has the original comments and headers, and is therefore a valid VCF file.
        The blocks contain variants from windows of positions so that batches with differing variants do not mess up the subsequent merging (e.g. legacy batches and chr X PAR).
        Files are bgzipped and can be indexed.
        """
        import gzip
        import bgzip
        import argparse
        import os.path as path

        def main(path:str, blocksize:int, prefix: str):
            
            block_names = []

            with gzip.open(path,"rt") as f:
                
                #read comments, header
                comment_holder = []
                header = None
                while True:
                    l = f.readline()
                    if l[0:2] == '##':
                        comment_holder.append(l)
                    elif l[0:2] == '#C':
                        header=l
                        break
                    else:
                        raise Exception("Header line not present before data! Invalid vcf")
                #after that, add comments & header to first block, then read data until block is full

                block_idx = 0
                last_line = None
                last_pos = -1
                file_read = False
                #write datalines
                while not file_read:
                    block_name = f"{prefix}_{block_idx}.gz"
                    block_names.append(block_name)
                    print(f"writing {block_name}")
                    #start writing block
                    with open(block_name, "wb") as blockfile:
                        with bgzip.BGZipWriter(blockfile) as fh:
                            for comment_line in comment_holder:
                                fh.write(bytes(comment_line, encoding='ascii'))
                            fh.write(bytes(header, encoding='ascii'))
                            #write last line from previous iteration if it was not written yet
                            if last_pos >= (blocksize * block_idx) and last_pos < (blocksize * (block_idx+1)):
                                fh.write(bytes(last_line, encoding='ascii'))
                            #write datalines
                            last_line = None
                            for l in f:
                                last_line = l
                                last_pos = int(last_line.strip().split('\t')[1])
                                if last_pos >= (blocksize * (block_idx+1)):
                                    break
                                fh.write(bytes(last_line, encoding='ascii'))
                            if last_line is None: #EOF
                                file_read=True
                    block_idx += 1
            print("Wrote blocks:")
            for b in block_names:
                print(b)
            

        def validate_path(p: str) -> str :
            if path.exists(p):
                return p
            raise Exception(f"path {p} invalid: File not found")

        if __name__=="__main__":
            parser = argparse.ArgumentParser("Separate VCF into chunks, writing comments & header to each of them.")
            parser.add_argument("file",type=validate_path,help="path to vcf")
            parser.add_argument("--block-size",type=int, help="interval size for blocks")
            parser.add_argument("--prefix",type=str,help="block prefix")
            args = parser.parse_args()
            main(args.file, args.block_size, args.prefix)
        EOF

        python3 blocker.py ~{batch} --block-size ~{chunksize} --prefix ~{prefix}

        ls -1 ~{prefix}_*.gz | sort -V > chunks_written
    >>>

    runtime {
        docker: docker
        memory: "2 GB"
        cpu: 1
        disks: "local-disk " + (ceil(size(batch, "GB")) * 2 + 3) + " HDD"
        preemptible: 2
        zones: "europe-west1-b europe-west1-c europe-west1-d"
    }

    output {
        Array[File] batch_chunked = read_lines("chunks_written")
        Array[File] globbed = glob(prefix+"_*.gz")
    }
}

task merge {

    input {
        Array[File] chunks
        Int number
        Int cpus = 10

        String docker
    }
    
    command <<<
        set -euxo pipefail

        echo "Indexing chunks"
        parallel tabix -p vcf ::: ~{sep=" " chunks}

        echo "Merge chunks"
        if [ ~{number} -eq 0 ]; then
            time bcftools merge -Ou ~{sep=" " chunks} | \
            bcftools +fill-tags -Ou -- -t AF,AN,AC,AC_Hom,AC_Het,NS | \
            bcftools +impute-info -Oz -o merged_chunk~{number}.gz
        else
            time bcftools merge -Ou ~{sep=" " chunks} | \
            bcftools +fill-tags -Ou -- -t AF,AN,AC,AC_Hom,AC_Het,NS | \
            bcftools +impute-info -Ou | \
            bcftools view --no-header -Oz -o merged_chunk~{number}.gz
        fi
    >>>

    output {
        File out = "merged_chunk" + number + ".gz"
    }

    runtime {
        docker: docker
        memory: "8 GB"
        cpu: cpus
        disks: "local-disk " + (ceil(size(chunks, "GB")) * 2 + 3) + " HDD"
        preemptible: 0
        zones: "europe-west1-b europe-west1-c europe-west1-d"
    }
}

task concatenate_chunks {

    input {
        Array[File] chunks
        String filename
        String docker
    }

    command {
        set -euxo pipefail
        #move chunks here so cromwell id etc do not mess up sorting
        mv -t ./ ~{sep=" " chunks}
        
        find ./ -name "merged_chunk*.gz" | sort -V | xargs cat > ${filename}

        tabix -p vcf ~{filename}
    }

    output {
        File out = filename
        File idx = filename + ".tbi"
    }

    runtime {
        docker: docker
        memory: "3 GB"
        cpu: 1
        disks: "local-disk " + (ceil(size(chunks, "GB")) * 2 + 3) + " SSD"
        preemptible: 0
        zones: "europe-west1-b europe-west1-c europe-west1-d"
    }
}
