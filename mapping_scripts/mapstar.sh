#!/bin/bash

if [ $# -ne 4 ]
then
    echo "Please, give:"
    echo "1) fastq file to map"
    echo "2) root for output file (no .sam or .bam extension)"
    echo "3) reference genome "
    echo "4) path to STAR"
    exit
fi

file2map=$1
outfq=$2
ref=$3
path2star=$4

${path2star}/STAR --runThreadN 12 --genomeDir $ref --readFilesIn ${file2map} --readFilesCommand zcat --outFileNamePrefix ${outfq} --outSAMtype BAM SortedByCoordinate --outSAMattributes All --outSAMstrandField intronMotif --outFilterMultimapNmax 1
rm -r ${outfq}_STARtmp
