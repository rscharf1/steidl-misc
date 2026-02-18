#!/bin/bash
#SBATCH --job-name=coding
#SBATCH --time=3:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --ntasks=1

module load bcftools

# awk '$3=="exon"' ~/Tools/star_files/index/gencode.v44.annotation.gtf > exons.gtf

awk '$3=="exon" {print $1"\t"$4-1"\t"$5}' \
  ~/Tools/star_files/index/gencode.v44.annotation.gtf > exons.bed

bcftools view \
  -R exons.bed \
  joint.qual30.vcf.gz \
  -Oz -o joint.qual30.exonic.vcf.gz

tabix -p vcf joint.qual30.exonic.vcf.gz