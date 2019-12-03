#!/usr/bin/env python3
import sys, os
import glob
import gzip
import numpy as np

try:
    fq = sys.argv[1]
except:
    sys.exit("Please, give root for R1, I1 and R2 fastq.gz files")

fqr1 = glob.glob('*'.join([fq, 'R1', 'fastq.gz']))[0]
fqr2 = glob.glob('*'.join([fq, 'R2', 'fastq.gz']))[0]

foutfile = fq + '_cbc.fastq'
fout = open(foutfile, 'w')

with gzip.open(fqr1) as fr1, gzip.open(fqr2) as fr2:
    for idx, (l1, l2) in enumerate(zip(fr1, fr2)):
        l1, l2 = [str(l.rstrip().rsplit()[0], 'utf-8') for l in [l1, l2]]
        l = np.mod(idx,4)
        if l == 0:
            n = l1
            if not l1 == l2:
                sys.exit('fastq files not synchronized @name: ' + ' '.join([l1, l2]))
        elif l == 1:
            s1, s2 = l1, l2
        elif l == 2:
            if not l1 == l2  == '+':
                sys.exit('fastq files not synchronized (+)')
        elif l == 3:
            cell = s1[:16]
            umi = s1[16:16+12]
            cellq = l1[:16]
            umiq = l1[16:16+12]
            name = ':'.join([n, cell, umi, cellq, umiq])

            fout.write('\n'.join([name,s2,'+',l2,'']))

fout.close()
os.system('gzip '+foutfile)
