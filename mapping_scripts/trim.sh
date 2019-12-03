#!/bin/bash

if [ $# -ne 3 ]
then
  echo "Please, give (1) input fastq file; (2) path2trimgalore; (3) path2cutadapt"
  exit
fi

file2trim=$1
path2trimgalore=$2
path2cutadapt=$3

${path2trimgalore}/trim_galore --path_to_cutadapt ${path2cutadapt}/cutadapt ${file2trim}
mv ${file2trim%.fastq.gz_trimming_report.txt}_trimming_report.txt
