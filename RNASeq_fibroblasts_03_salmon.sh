#!/bin/bash
#PBS -l nodes=1:ppn=16
#PBS -l mem=48gb
#PBS -N salmon
#PBS -d /beegfs/work/workspace/ws/ul_ypb85-MTM_fibro-0/mapped

module purge
for fn in 813{877..933};
do
samp=`basename ${fn}`
echo "Processing sample ${samp}"
/home/ul/ul_neurorku/ul_ypb85/salmon-0.14.1/bin/salmon quant -t /beegfs/work/ul_ypb85/GRCh38_r90_transcriptome/Homo_sapiens.GRCh38.90.transcriptome.fa -l A -a ${fn}/Aligned.toTranscriptome.out.bam -g /beegfs/work/ul_ypb85/GRCh38_r90/Homo_sapiens.GRCh38.90.gtf --seqBias --gcBias --numBootstraps 100 -o /beegfs/work/workspace/ws/ul_ypb85-MTM_fibro-0/salmon/${samp}
done