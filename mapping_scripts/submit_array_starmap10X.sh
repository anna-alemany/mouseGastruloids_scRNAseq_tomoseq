#!/bin/bash

if [ $# -ne 2 ]
then
    echo "Please, give input (1) root file; (2) out"
    exit
fi

p2s=/hpc/hub_oudenaarden/aalemany/bin/mapandgo2
p2trimgalore=/hpc/hub_oudenaarden/aalemany/bin/TrimGalore-0.4.3
p2cutadapt=/hpc/hub_oudenaarden/aalemany/bin
p2star=/hpc/hub_oudenaarden/avo/nascent/STAR-2.5.3a/bin/Linux_x86_64

source /hpc/hub_oudenaarden/aalemany/virtualEnvironments/venv36/bin/activate

email=a.alemany@hubrecht.eu

in=$1
out=$2
protocol=10x
reference=mouse

# merge data
echo "${p2s}/mergeLanes.sh $in $out" | qsub -cwd -m eas -M $email -N merge-$out -e merge-${out}.err -o merge-${out}.out -l h_vmem=10G -l h_rt=15:00:00 -pe threaded 2

# extract barcodes
echo "${p2s}/fastqfile_demuxXindex_10x.py $out" | qsub -V -cwd -m eas -M $email -N extract-$out -e extract-${out}.err -o extract-${out}.out -l h_vmem=10G -l h_rt=15:00:00 -hold_jid merge-$out

# trim
echo "${p2s}/trim.sh ${out}_cbc.fastq.gz ${p2trimgalore} ${p2cutadapt}" | qsub -cwd -m eas -M $email -N trim-$out -e trim-${out}.err -o trim-${out}.out -l h_vmem=10G -l h_rt=15:00:00 -hold_jid extract-$out

# map with star
echo "${p2s}/mapstar.sh ${out}_cbc_trimmed.fq.gz ${out}_cbc_trimmed_star $reference ${p2star}" | qsub -cwd -m eas -M $email -N map-$out -e map-${out}.err -o map-${out}.out -l h_vmem=30G -l h_rt=15:00:00 -pe threaded 12 -hold_jid trim-$out

# extract cell numbers
echo "zcat ${out}_cbc_trimmed.fq.gz | awk '{if (FNR%4==1) print $0}' | awk -F ':' '{print $8}' | sort | uniq -c | sort -k1nr > ${out}_cbc_trimmed_rawCELLS.txt" | qsub -cwd -m eas -M $email -N rcell-$out -e rcell-${out}.err -o rcell-${out}.out -l h_vmem=30G -l h_rt=15:00:00 -pe threaded 12 -hold_jid map-$out
