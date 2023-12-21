workflow merge_in_chunks{
    Array[File] vcfs
    String outfile
    Int chunksize
    String docker

    scatter (i in range(length(vcfs))){
        call chunk {
            input: batch=vcfs[i],chunksize=chunksize
        }
    }
    call reorder_array{
        input: batch_chunk_array = chunk.globbed
    }
    Array[Array[File]] chunk_batch_array = transpose(reorder_array.batch_chunk)

    scatter (i in range( length( chunk_batch_array ) ) ){
        call merge {
            input: chunks = chunk_batch_array[i],number=i,docker=docker
        }
    }

    call concatenate_chunks {
        input: chunks = merge.out, docker = docker, filename = outfile
    }
    output {
        File out = concatenate_chunks.out
        File idx = concatenate_chunks.idx
    }
}

task reorder_array{

    Array[Array[String]] batch_chunk_array
    File tsv = write_tsv(batch_chunk_array)

    command {
        set -euxo pipefail
        python3 <<EOF
        #read data
        with open("${tsv}","r") as f:
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
        preemptible: 0
        zones: "europe-west1-b europe-west1-c europe-west1-d"
    }

    output {
        Array[Array[String]] batch_chunk = read_tsv("output")
    }
}

task chunk{
    File batch
    Int chunksize
    String prefix = basename(batch,".vcf.gz")

    command <<<
        set -euxo pipefail
        cat << EOF > blocker.py
        """
        Separate a VCF file into N-lines long gzipped blocks. Each block has the original comments and headers, and is therefore a valid VCF file. 
        NOTE: Since the files a gzipped, NOT bgzipped, they CAN NOT be indexed. This script should be used as only an intermediate step before vcf-fusion. 
        """
        import gzip
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
                file_read = False
                #write datalines
                while not file_read:
                    block_name = f"{prefix}_{block_idx}.gz"
                    block_names.append(block_name)
                    print(f"writing {block_name}")
                    #start writing block
                    with gzip.open(block_name,"wt") as blockfile:
                        
                        for comment_line in comment_holder:
                            blockfile.write(comment_line)
                        blockfile.write(header)
                        #write datalines
                        current_dataline=0
                        last_line=""
                        for l in f:
                            blockfile.write(l)
                            current_dataline +=1
                            last_line=l
                            if current_dataline >= blocksize:
                                break
                        if last_line=="":
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
            parser = argparse.ArgumentParser("Separate VCF into chunks of N datalines, writing comments & header to each of them.")
            parser.add_argument("file",type=validate_path,help="path to vcf")
            parser.add_argument("--block-size",type=int, help="how many data lines at most in one block")
            parser.add_argument("--prefix",type=str,help="block prefix")
            args = parser.parse_args()
            main(args.file, args.block_size, args.prefix)
        EOF
        python3  blocker.py ${batch} --block-size ${chunksize} --prefix ${prefix}
        ls -1 ${prefix}_*.gz|sort -V > chunks_written
    >>>

    runtime {
        docker: "python:3.6-slim"
        memory: "2 GB"
        cpu: 1
        disks: "local-disk 100 HDD"
        preemptible: 0
        zones: "europe-west1-b europe-west1-c europe-west1-d"
    }

    output {
        Array[File] batch_chunked = read_lines("chunks_written")
        Array[File] globbed = glob(prefix+"_*.gz")
    }
}

task merge{
    Array[File] chunks
    Int number
    String opts = if number == 0 then ""  else "--no-comments --no-header" 
    Int cpus = 10

    String docker

    command <<<
        set -euxo pipefail
        echo "Execute vcf-fusion&paste&bgzip command"
        time vcf-fusion ${opts} ${sep=" " chunks} | bgzip -@${cpus} > merged_chunk${number}.gz
    >>>

    output {
        File out = "merged_chunk" + number + ".gz"
    }

    runtime {
        docker: docker
        memory: "8 GB"
        cpu: cpus
        disks: "local-disk 100 HDD"
        preemptible: 0
        zones:"europe-west1-b europe-west1-c europe-west1-d"
    }
}

task concatenate_chunks{
    Array[File] chunks
    String filename
    String docker

    command{
        set -euxo pipefail
        #move chunks here so cromwell id etc do not mess up sorting
        mv -t ./ ${sep=" " chunks}
        
        find ./ -name "merged_chunk*.gz" | sort -V | xargs cat > ${filename} 

        tabix -p vcf ${filename}

    }

    output {
        File out = filename
        File idx = filename + ".tbi"
    }

    runtime {
        docker: docker
        memory: "3 GB"
        cpu: 1
        disks: "local-disk 600 SSD"
        preemptible: 0
        zones: "europe-west1-b europe-west1-c europe-west1-d"
    }
}