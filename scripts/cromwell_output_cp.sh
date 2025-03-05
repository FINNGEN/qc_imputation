#!/bin/bash

set -euo pipefail

hash=$1
out=$(echo $2 | sed 's|/*$|/|')

for f in $(ls -1 ${hash}*)
do
    n_files=$(grep "^gs://" ${f} | wc -l)
    if [ $n_files -eq 0 ];
    then
        echo "No outputs in ${f}"
        continue
    fi
    subdir=$(echo $f | rev | cut -d "." -f 2 | rev)
    outdir=${out}${subdir}/
    echo "Transferring $n_files outputs in ${f} to ${outdir}"
    cat ${f} | gsutil -mq cp -I ${outdir}
done