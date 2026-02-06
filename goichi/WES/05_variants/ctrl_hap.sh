#!/bin/bash
#SBATCH --job-name=ctrl_hap
#SBATCH --time=3:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=logs/haplotype_%A_%a.out
#SBATCH --error=logs/haplotype_%A_%a.err

set -euo pipefail

module load gatk

REF=~/Tools/bwa_files/GRCh38.primary_assembly.genome.fa

BAM="../04_dedup/bams/HL60_ctrl.fixedRG.dedup.bam"

OUT="gvcfs/HL60_ctrl.fixedRG.g.vcf.gz"

gatk HaplotypeCaller \
  -R ${REF} \
  -I ${BAM} \
  -O ${OUT} \
  -ERC GVCF