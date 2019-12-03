#!/bin/bash

if [ $# -ne 2 ]
then
  echo "Please, give root for input and output fastq files"
  exit
fi

fq=$1
outfq=$2
zcat ${fq}*R1* > ${outfq}_R1.fastq &
zcat ${fq}*R2* > ${outfq}_R2.fastq &
wait

gzip ${outfq}_R1.fastq &
gzip ${outfq}_R2.fastq
