#!/bin/bash
#SBATCH --job-name=counts
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=logs/%A_%a.out
#SBATCH --array=1-6

set -euo pipefail

INPUT=$(sed -n "${SLURM_ARRAY_TASK_ID}p" input_bams.txt)
SAMPLE=$(basename "$INPUT" .shifted.sorted.bam)

echo $INPUT
echo $SAMPLE

mkdir -p counts

bedtools coverage \
  -sorted \
  -g hg38.genome \
  -a merged_peaks_500bp.filtered.sorted.bed \
  -b ${INPUT} \
  -counts > counts/${SAMPLE}.counts.txt