#!/bin/bash
#SBATCH --job-name=cutadapt
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/cutadapt_%A_%a.out
#SBATCH --array=1-2

set -euo pipefail

# Get R1 for this array task
R1=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples.txt)

# Derive R2 by pattern substitution
R2=${R1/_1.fq.gz/_2.fq.gz}

# Extract sample name
SAMPLE=$(basename "$R1" _1.fq.gz)

mkdir -p trimmed logs

cutadapt \
  -a AGATCGGAAGAGC \
  -A AGATCGGAAGAGC \
  -m 30 \
  -q 20 \
  --trim-n \
  --cores ${SLURM_CPUS_PER_TASK} \
  -o trimmed/${SAMPLE}_R1.trimmed.fastq.gz \
  -p trimmed/${SAMPLE}_R2.trimmed.fastq.gz \
  "$R1" "$R2"

