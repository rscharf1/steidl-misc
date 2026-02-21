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
  --bam bams/HL_Ctrl_CKDN240010648-1A_227JFTLT4_L6.dedup.bam \
  --bam bams/HL_AraC_CKDN240010649-1A_227JFTLT4_L6.dedup.bam \
  --reference GRCh38.primary_assembly.genome.fa \
  --runDir manta_joint