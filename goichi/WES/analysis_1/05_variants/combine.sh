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
  --variant gvcfs/HL60_ctrl.fixedRG.g.vcf.gz  \
  --variant gvcfs/PY_high_1_CKDN260000485-1A_23KL3LLT4_L4.g.vcf.gz \
  --variant gvcfs/PY_low_1_CKDN260000486-1A_23KL3LLT4_L4.g.vcf.gz \
  -O joint/combined.g.vcf.gz