#!/bin/bash
#SBATCH --job-name=peak_call
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=16:00:00
#SBATCH --output=logs/peaks_%A_%a.out
#SBATCH --array=1-10

set -euo pipefail

# conda activate macs2

INPUT=$(sed -n "${SLURM_ARRAY_TASK_ID}p" input_bams.txt)
SAMPLE=$(basename "$INPUT" .shifted.sorted.bam)

echo $INPUT
echo $SAMPLE

macs2 callpeak \
  -t ${INPUT} \
  -f BAMPE \
  -g hs \
  --keep-dup all \
  --call-summits \
  -n ${SAMPLE} \
  --outdir peaks