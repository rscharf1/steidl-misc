#!/bin/bash
#SBATCH --job-name=genotype
#SBATCH --time=3:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --ntasks=1

module load gatk

gatk GenotypeGVCFs \
  -R ~/Tools/bwa_files/GRCh38.primary_assembly.genome.fa \
  -V joint/combined.g.vcf.gz \
  -O joint/joint.raw.vcf.gz