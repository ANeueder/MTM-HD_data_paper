#!/bin/bash
#PBS -l nodes=1:ppn=16
#PBS -l mem=48gb
#PBS -N Star_1
#PBS -d /beegfs/work/workspace/ws/ul_ypb85-MTM_fibro-0/fastq

module purge
module load bio/star/2.5.3a-gnu-4.9

for fn in 813{877..933};
do
mkdir -p /beegfs/work/workspace/ws/ul_ypb85-MTM_fibro-0/mapped/${fn}/
done

for fn in 813{877..895};
do
samp=`basename ${fn}`
echo "Processing sample ${samp}"
STAR --runThreadN 16 --genomeDir /beegfs/work/ul_ypb85/GRCh38_r90_star_100overhang --readFilesIn ${samp}_1.fastq.gz ${samp}_2.fastq.gz --outFileNamePrefix /beegfs/work/workspace/ws/ul_ypb85-MTM_fibro-0/mapped/${fn}/ --outSAMtype BAM Unsorted SortedByCoordinate --quantMode TranscriptomeSAM --twopassMode Basic --readFilesCommand zcat --limitBAMsortRAM 31532137230 --outReadsUnmapped Fastx --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --outFilterMultimapNmax 20 --alignIntronMin 20 --alignIntronMax 100000 --alignSJDBoverhangMin 3 --chimSegmentMin 12 --chimJunctionOverhangMin 12 --alignMatesGapMax 100000 --chimSegmentReadGapMax 3 --alignSJstitchMismatchNmax 5 -1 5 5
done

