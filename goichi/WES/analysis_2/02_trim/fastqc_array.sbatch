#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=logs/fastqc_%A_%a.out
#SBATCH --array=1-4

FASTQ=$(sed -n "${SLURM_ARRAY_TASK_ID}p" trimmed_fastq_files.txt)

mkdir -p fastqc_out
fastqc -o fastqc_out "$FASTQ"
