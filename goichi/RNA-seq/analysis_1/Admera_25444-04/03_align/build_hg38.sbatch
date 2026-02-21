#!/bin/bash
#SBATCH --job-name=bt2_build_hg38
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/bowtie2_build_hg38_%j.out

set -euo pipefail

module load bowtie2

# Move to directory with hg38.fa
cd genomes/hg38

# Sanity check
ls -lh hg38.fa

bowtie2-build \
  --threads ${SLURM_CPUS_PER_TASK} \
  hg38.fa \
  hg38
