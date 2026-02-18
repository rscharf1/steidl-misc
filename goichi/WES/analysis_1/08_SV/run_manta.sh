#!/bin/bash
#SBATCH --job-name=Manta
#SBATCH --time=4:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --ntasks=1

# module load samtools 
# samtools faidx GRCh38.primary_assembly.genome.fa

configManta.py \
  --bam bams/HL60_ctrl.fixedRG.dedup.bam \
  --bam bams/PY_low_1_CKDN260000486-1A_23KL3LLT4_L4.dedup.bam \
  --bam bams/PY_high_1_CKDN260000485-1A_23KL3LLT4_L4.dedup.bam \
  --reference GRCh38.primary_assembly.genome.fa \
  --runDir manta_joint