#!/bin/bash
#SBATCH --job-name=combine
#SBATCH --time=3:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --ntasks=1

module load gatk

mkdir -p joint 

gatk CombineGVCFs \
  -R ~/Tools/bwa_files/GRCh38.primary_assembly.genome.fa \
  --variant gvcfs/HL_Ctrl_CKDN240010648-1A_227JFTLT4_L6.g.vcf.gz  \
  --variant gvcfs/HL_AraC_CKDN240010649-1A_227JFTLT4_L6.g.vcf.gz \
  -O joint/combined.g.vcf.gz