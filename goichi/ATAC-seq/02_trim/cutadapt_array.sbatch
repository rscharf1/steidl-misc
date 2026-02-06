#!/bin/bash
#SBATCH --job-name=cutadapt
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/cutadapt_%A_%a.out
#SBATCH --array=0-9

set -euo pipefail

# Get R1 for this array task
# R1=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples.txt)
R1=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" samples.txt)

# Derive R2 by pattern substitution
R2=${R1/_R1_001.fastq.gz/_R2_001.fastq.gz}

# Extract sample name
SAMPLE=$(basename "$R1" _R1_001.fastq.gz)

echo "Sample: ${SAMPLE}"
echo "R1: ${R1}"
echo "R2: ${R2}"

mkdir -p trimmed logs

cutadapt \
  -a CTGTCTCTTATACACATCT \
  -A CTGTCTCTTATACACATCT \
  -m 30 \
  -q 20 \
  --trim-n \
  --cores ${SLURM_CPUS_PER_TASK} \
  -o trimmed/${SAMPLE}_R1.trimmed.fastq.gz \
  -p trimmed/${SAMPLE}_R2.trimmed.fastq.gz \
  "$R1" "$R2"

