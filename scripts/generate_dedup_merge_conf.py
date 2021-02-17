#!/usr/bin/env python

import requests
import json
from typing import List
import io
import subprocess
import shlex
import argparse

cromwell_ip="0.0.0.0"

def get_qc_import_outfiles(hash):
    print(f'Asking cromwell at {cromwell_ip} for status of {hash}')
    query=(f'http://{cromwell_ip}/api/workflows/v1/{hash}/metadata'
            '?includeKey=outputs&includeKey=status&includeKey=executionStatus'
            '&includeKey=failures&includeKey=start&includeKey=end&includeKey=stdout'
            '&includeKey=inputs'
          )
    proxies = { 'http': 'socks5://localhost:5000',
                    'https': 'socks5://localhost:5000'}
    res = None

    res = requests.post(query, proxies=proxies,timeout=60)
    if not res.ok:
        res.raise_for_status()
    res_json = json.loads(res.text)
    if res_json["status"] != "Succeeded":
        raise Exception("process data only after workflows are finished. Status of " + hash + " is " +  res_json["status"])

    sample_summary_file = res_json["outputs"]["qc_imputation.joint_plots.sample_summaries"]
    per_chr_vcfs = res_json["outputs"]["qc_imputation.subset_samples.subset_vcfs"]
    return(sample_summary_file, per_chr_vcfs)

def dups_to_remove( sumstats, dupfile):
    ind_dat={}
    removals = []
    with open(dupfile, 'rt') as d:
        for l in d:
            l=l.strip("\n").split("\t")
            if (l[0].find("#")!=-1):
                ## ignore comment lines
                continue
            best_score=0
            best_id=""
            for id in l:
                if not id  in sumstats:
                    ## does not exist in data but indicate as removal just in case the sumdatfile is incomplete
                    #print("no need to remove " + id)
                    #removals.append(id)
                    continue

                score = int(sumstats[id]["n_vars_total"]) - int(sumstats[id]["n_excluded_vars"]) - int(sumstats[id]["N_MISS_in_batch_qc"])
                if score>best_score:
                    if best_id != "":
                        removals.append(best_id)
                    best_score=score
                    best_id = id
                else:
                    removals.append(id)
    return removals

def get_bucket_data_as_text(path):
    cmd=f'gsutil cat {path}'
    pr = subprocess.run(shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE,encoding="ASCII")
    if pr.returncode!=0:
        print(pr.stderr)
        raise Exception('Error occurred while accessing data from bucket: ' + path)

    return pr.stdout


def get_sumstats(sumstatpaths:List[str]):
    sumdat = {}
    for s in sumstatpaths:
        sumfile_handle = io.StringIO(get_bucket_data_as_text(s))
        header = { h:i for i,h in enumerate(sumfile_handle.readline().strip().split("\t")) }
        for l in sumfile_handle:
            l = l.strip().split("\t")
            d = { k:l[i] for k,i in header.items() }
            if d["IID"] in sumdat:
                print("Warning same ID in multiple times in summary files!" + d["IID"] +
                    "Old pipeline had all samples across all batches but modified to add only finally included samples")
            sumdat[d["IID"]] = d
    return sumdat


def run():
    parser = argparse.ArgumentParser(description="Generate conf for deduplicating and merging multiple qc-imputation runs")
    parser.add_argument('hashes', type=str, help='Comma separated list of cromwell hashes')
    parser.add_argument('dup_file',  type=str, help='File with duplicate ID samples, multiple IDs of one individual on one line separated b ytabs')
    parser.add_argument('outprefix',  type=str, help='output prefix')
    args = parser.parse_args()

    jobs = args.hashes.split(",")

    if len(jobs)<=1:
        print("No point in running single jobs depup merge")
        exit(0)

    jdats=[]
    for j in jobs:
        jdats.append(get_qc_import_outfiles(j))

    sumstats = get_sumstats([ jdat[0] for jdat in jdats])

    remove_id_list = dups_to_remove( sumstats, args.dup_file)

    per_chr_files = args.outprefix + "_all_batches_per_chr"

    remove_file = args.outprefix + "_final_sample_removals"

    with open(remove_file,"w") as rf:
        for id in remove_id_list:
            rf.write(id+"\n")

    prev=None
    ## assert that the same number of chromosomes were ran in each.
    for j in jdats:
        if prev is not None and prev != len(j[1]):
            raise Exception("somethings wrong as there are different number of output chromosomes in different qc/impute runs")

    with open(per_chr_files,"w") as batches:
        for i in range(len(jdats[0][1])):
            chr_batches = []
            for j in jdats:
                chr_batches.extend(j[1][i])
            batches.write("\t".join(chr_batches) + "\n")

    print("Ok now run wdl/merge_dedup_qc_imps.wdl with wdl/dedup.wdl dependency." +
            f"Copy these ({per_chr_files}, {remove_file}) files to Google bucket and point wd/merge_dedup_qc_imps.json to them." )

    ## zip -j wdl/merge_dedup_qc_imps_deps.zip wdl/dedup.wdl
    ## 

if __name__ == "__main__":
    run()
