#!/bin/bash
#PBS -l nodes=1:ppn=16
#PBS -l mem=48gb
#PBS -N Star_map_batch1
#PBS -d /beegfs/work/workspace/ws/ul_ypb85-MTM_RNA_first-0/fastq/batch1

module purge
module load bio/star/2.5.3a-gnu-4.9

for fn in {387935..387994};
do
mkdir -p /beegfs/work/workspace/ws/ul_ypb85-MTM_RNA_first-0/mapped/batch1/${fn}/
done

for fn in {387935..387994};
do
samp=`basename ${fn}`
echo "Processing sample ${samp}"
STAR --runThreadN 16 --genomeDir /beegfs/work/ul_ypb85/GRCh38_r90_star_50overhang --readFilesIn ${samp}_1.clipped.fastq.gz ${samp}_2.clipped.fastq.gz --outFileNamePrefix /beegfs/work/workspace/ws/ul_ypb85-MTM_RNA_first-0/mapped/batch1/${fn}/ --outSAMtype BAM Unsorted SortedByCoordinate --quantMode TranscriptomeSAM --twopassMode Basic --readFilesCommand zcat --limitBAMsortRAM 31532137230 --outReadsUnmapped Fastx --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --outFilterMultimapNmax 20 --alignIntronMin 20 --alignIntronMax 100000 --alignSJDBoverhangMin 3 --chimSegmentMin 12 --chimJunctionOverhangMin 12 --alignMatesGapMax 100000 --chimSegmentReadGapMax 3 --alignSJstitchMismatchNmax 5 -1 5 5
done

