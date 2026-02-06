#!/bin/bash
#SBATCH --job-name=bcf
#SBATCH --time=3:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --ntasks=1

module load bcftools

# bcftools view \
#   -f PASS \
#   -i 'QUAL>=30' \
#   ../05_variants/joint/joint.raw.vcf.gz \
#   -Oz -o joint.pass.vcf.gz

# tabix -p vcf joint.pass.vcf.gz


bcftools view \
  -i 'QUAL>=30' \
  ../05_variants/joint/joint.raw.vcf.gz \
  -Oz -o joint.qual30.vcf.gz

tabix -p vcf joint.qual30.vcf.gz