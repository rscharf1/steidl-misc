#!/bin/bash
#SBATCH --job-name=Manta
#SBATCH --time=4:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --ntasks=1

bcftools query \
  -f '%CHROM\t%POS\t%SVTYPE\t%FILTER[\t%GT]\n' \
  manta_joint/results/variants/diploidSV.vcf.gz \
  > diploidSV_table.tsv