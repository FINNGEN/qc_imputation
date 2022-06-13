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