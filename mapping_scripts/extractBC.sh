#!/bin/bash

if [ $# -ne 2 ]
then
    echo "Please, give: "
    echo "1) input root to fastq files"
    echo "2) path to concatenator.py"
    exit
fi

outfq=$1
protocol=celseq2
path2scripts=$2

python3 ${path2scripts}/concatenator.py --fqf ${outfq} --cbcfile ${path2scripts}/bc_celseq2.tsv --cbchd 0 --lenumi 6 --umifirst
gzip ${outfq}_cbc.fastq

