#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l mem=16gb
#PBS -N FastQC
#PBS -d /beegfs/work/workspace/ws/ul_ypb85-MTM_fibro-0/fastq

module purge
echo $PWD
module load bio/fastqc/0.11.5
for i in $(ls *.fastq.gz)
do
fastqc --outdir=/beegfs/work/workspace/ws/ul_ypb85-MTM_fibro-0/fastQC/ -t 8 $i
done