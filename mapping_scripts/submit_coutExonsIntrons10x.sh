#!/bin/bash

if [ $# -ne 4 ]
then
    echo "Please, give:"
    echo "1) intron bed file"
    echo "2) exon bed file"
    echo "3) selected barcodes (cells) file"
    echo "4) root for output file"
    exit
fi

# in case you need to activate python3 virtual environment
source /hpc/hub_oudenaarden/aalemany/virtualEnvironments/venv36/bin/activate

# make sure to have the script somewhere
echo "./countExonsIntrons_10x.py $1 $2 $3 $4" | qsub -V -cwd -l h_rt=10:00:00 -l h_vmem=80G -N tab-$4 -m eas -M a.alemany@hubrecht.eu -e tab-${4}.err -o tab-${4}.out
