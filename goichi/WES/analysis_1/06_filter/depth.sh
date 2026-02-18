#!/bin/bash
#SBATCH --job-name=depth
#SBATCH --time=3:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --ntasks=1

module load bcftools

bcftools filter \
  -i 'MIN(FMT/DP)>=20' \
  joint.qual30.exonic.vcf.gz \
  -Oz -o joint.exonic.dp20.vcf.gz

tabix -p vcf joint.exonic.dp20.vcf.gz