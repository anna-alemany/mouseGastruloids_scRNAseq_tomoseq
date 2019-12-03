#!/bin/bash

if [ $# -ne 2 ]
then
    echo "Please, give following inputs:"
    echo "(1) root for R1/R1 fastq.gz files (common part of the name)"
    echo "(2) name for output files (no extension needed)"
    exit
fi

# define paths to required software. This will change for each user. 
p2s=/hpc/hub_oudenaarden/aalemany/bin/mapandgo2                 # path to the folder where mapping scripts are kept
p2trimgalore=/hpc/hub_oudenaarden/aalemany/bin/TrimGalore-0.4.3 # path to TrimGalore-0.4.3
p2cutadapt=/hpc/hub_oudenaarden/aalemany/bin                    # path to cutadapt
p2star=/hpc/hub_oudenaarden/avo/nascent/STAR-2.5.3a/bin/Linux_x86_64 # path to STAR
p2bedtools=/hpc/hub_oudenaarden/aalemany/bin/bedtools2/bin      # path to bedtools
p2samtools=/hpc/hub_oudenaarden/bdebarbanson/bin/samtools-1.3.1 # path to samtools

# if necessary, activate python3 virtual environemnt
source /hpc/hub_oudenaarden/aalemany/virtualEnvironments/venv36/bin/activate

# define user's email
email=a.alemany@hubrecht.eu

# define path to mouse reference genome and the path to bedfiles containing intron/exon locations
reference=/hpc/hub_oudenaarden/group_references/ensembl/93/mus_musculus/star_index_60
intron=/hpc/hub_oudenaarden/group_references/ensembl/93/mus_musculus/annotations_ensembl_93_mm_introns_exonsubtracted.bed
exon=/hpc/hub_oudenaarden/group_references/ensembl/93/mus_musculus/annotations_ensembl_93_mm_exons.bed

# The real script starts here.
# This will submit a sequence of scripts in an array. First, fastq lanes will be merged, then barcodes will be extracted, resulting fastq files will be trimmed and mapped, and finally reads falling into introns and exons will be identified.
# Each script will produce *.err and *.out files summarizing possible error or output messages. 
in=$1
out=$2
protocol=celseq2

# merge data
echo "${p2s}/mergeLanes.sh $in $out" | qsub -cwd -m eas -M $email -N merge-$out -e merge-${out}.err -o merge-${out}.out -l h_vmem=10G -l h_rt=15:00:00 -pe threaded 2

# extract barcodes
echo "${p2s}/extractBC.sh $out ${protocol} ${p2s}" | qsub -V -cwd -m eas -M $email -N extract-$out -e extract-${out}.err -o extract-${out}.out -l h_vmem=10G -l h_rt=15:00:00 -hold_jid merge-$out

# trim
echo "${p2s}/trim.sh ${out}_cbc.fastq.gz ${p2trimgalore} ${p2cutadapt}" | qsub -cwd -m eas -M $email -N trim-$out -e trim-${out}.err -o trim-${out}.out -l h_vmem=10G -l h_rt=15:00:00 -hold_jid extract-$out

# map with star
echo "${p2s}/mapstar.sh ${out}_cbc_trimmed.fq.gz ${out}_cbc_trimmed_star $reference ${p2star}" | qsub -cwd -m eas -M $email -N map-$out -e map-${out}.err -o map-${out}.out -l h_vmem=30G -l h_rt=15:00:00 -pe threaded 12 -hold_jid trim-$out

# create count tables from star map
echo "${p2s}/getIntronsExons.sh ${out}_cbc_trimmed_starAligned.sortedByCoord.out.bam $intron $exon ${out}_cbc_trimmed_star ${p2bedtools} ${p2samtools} ${p2s}" | qsub -V -cwd -m eas -M $email -N tab-$out -e tab-${out}.err -o tab-${out}.out -l h_vmem=20G -l h_rt=15:00:00 -pe threaded 2 -hold_jid map-$out
