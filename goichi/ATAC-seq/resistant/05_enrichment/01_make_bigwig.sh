#!/bin/bash
#SBATCH --job-name=bigwig
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --array=1-6

# conda activate deeptools

set -euo pipefail

ls bam_symlinks/*bam > input_bams.txt

mkdir -p bw

INPUT=$(sed -n "${SLURM_ARRAY_TASK_ID}p" input_bams.txt)
SAMPLE=$(basename "$INPUT" .shifted.sorted.bam)

OUT_BAM="bw/${SAMPLE}.bw"

bamCoverage -b $INPUT \
  -o ${OUT_BAM} \
  --normalizeUsing CPM \
  --binSize 10 \
  --extendReads \
  -p 8

rm input_bams.txt



# mkdir -p bw

# for bam in bam_symlinks/*bam
# do
#   base=$(basename $bam .shifted.sorted.bam)
#   bamCoverage -b $bam \
#     -o bw/${base}.bw \
#     --normalizeUsing CPM \
#     --binSize 10 \
#     --extendReads \
#     -p 8
# done