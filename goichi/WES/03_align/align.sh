#!/bin/bash
#SBATCH --job-name=HL60_ctrl_L1_align
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=logs/HL60_ctrl_1_L1.%j.out
#SBATCH --error=logs/HL60_ctrl_1_L1.%j.err

set -euo pipefail

REF=~/Tools/bwa_files/GRCh38.primary_assembly.genome.fa

R1=../02_trim/trimmed/HL60_ctrl_1_CKDN260000487-1A_2335HJLT4_L1_R1.trimmed.fastq.gz
R2=../02_trim/trimmed/HL60_ctrl_1_CKDN260000487-1A_2335HJLT4_L1_R2.trimmed.fastq.gz

OUT=bams/HL60_ctrl_L1.bam

bwa mem -t 16 \
  -R "@RG\tID:HL60_ctrl_1_L1\tSM:HL60_ctrl_1\tPL:ILLUMINA" \
  $REF $R1 $R2 | \
samtools sort -@ 8 -o $OUT

samtools index $OUT
