#!/bin/bash
#SBATCH --job-name=merge_bams
#SBATCH --time=1:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --ntasks=1

# samtools merge -@ 8 \
# 	bams/HL60_ctrl.merged.bam \
# 	bams/HL60_ctrl_1_CKDN260000487-1A_2335HJLT4_L1.bam \
# 	bams/HL60_ctrl_1_CKDN260000487-1A_23KL3LLT4_L6.bam 

samtools index bams/HL60_ctrl.merged.bam