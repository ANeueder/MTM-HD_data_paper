#!/bin/bash
#PBS -l nodes=1:ppn=16
#PBS -l mem=48gb
#PBS -N salmon_quant1
#PBS -d /beegfs/work/workspace/ws/ul_ypb85-MTM_RNA_first-0/mapped/batch1

module purge
module load bio/salmon/0.8.2
for fn in {387935..387994};
do
samp=`basename ${fn}`
echo "Processing sample ${samp}"
salmon quant -t /beegfs/work/ul_ypb85/GRCh38_r90_transcriptome/Homo_sapiens.GRCh38.87.transcriptome.fa -l A -a ${fn}/Aligned.toTranscriptome.out.bam -g /beegfs/work/ul_ypb85/GRCh38_r90/Homo_sapiens.GRCh38.90.gtf --seqBias --gcBias --numBootstraps 100 -o /beegfs/work/workspace/ws/ul_ypb85-MTM_RNA_first-0/quant_salmon/batch1/${samp}
done