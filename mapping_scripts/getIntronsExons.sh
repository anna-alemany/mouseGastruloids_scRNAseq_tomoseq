#!/bin/bash

if [ $# -ne 6 ]
then
    echo "Please, give 4 input files:"
    echo "1) sorted bam file (map with star)"
    echo "2) genome [mouse, human, zebrafish, spiny]"
    echo "3) output root files"
    echo "4) path to bedtools"
    echo "5) path to samtools"
    echo "6) path to countExonsIntrons.py"
    exit
fi

inbam=$1
out=$3

intron=/hpc/hub_oudenaarden/group_references/ensembl/93/mus_musculus/annotations_ensembl_93_mm_introns_exonsubtracted.bed
exon=/hpc/hub_oudenaarden/group_references/ensembl/93/mus_musculus/annotations_ensembl_93_mm_exons.bed

p2b=$4
p2s=$5

${p2b}/bamToBed -i ${inbam} -split > ${out}_bam2bed.bed
${p2b}/bedtools intersect -a ${out}_bam2bed.bed -b ${intron} -wb | awk '{if (($6==$10) && ($5==255)) print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$NF}' | uniq > ${out}_intron1.bed &
${p2b}/bedtools intersect -a ${out}_bam2bed.bed -b ${exon} -wb | awk '{if (($6==$10) && ($5==255)) print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$NF}' | uniq > ${out}_exon1.bed
wait

chroms=$(${p2s}/samtools view -H $1 | grep '^@SQ' | awk -F 'SN:|\t' '{print $3}')

if [ -f ${out}_intron.bed ]
then
    rm ${out}_intron.bed
fi
if [ -f ${out}_exon.bed ]
then
    rm ${out}_exon.bed
fi

for i in $(echo $chroms | awk '{for (i=1; i<=NF; i++) print $i}')
do
    awk -v i=$i '{if ($1==i) print $0}' ${out}_intron1.bed > ${out}_TMP_intron1_${i}.bed &
    awk -v i=$i '{if ($1==i) print $0}' ${out}_exon1.bed > ${out}_TMP_exon1_${i}.bed
    wait

    ${p2b}/bedtools intersect -a ${out}_TMP_intron1_${i}.bed -b ${out}_TMP_exon1_${i}.bed -v >> ${out}_intron.bed &
    ${p2b}/bedtools intersect -a ${out}_TMP_exon1_${i}.bed -b ${out}_TMP_intron1_${i}.bed -v >> ${out}_exon.bed
    wait

    rm  ${out}_TMP_exon1_${i}.bed ${out}_TMP_intron1_${i}.bed
done

rm ${out}_bam2bed.bed
rm ${out}_intron1.bed ${out}_exon1.bed

p2s=$6
${p2s}/countExonsIntrons.py ${out}_intron.bed ${out}_exon.bed ${out}
