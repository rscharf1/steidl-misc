#!/bin/bash
#SBATCH --job-name=DeDup_GATK
#SBATCH --time=3:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=0-1
#SBATCH --output=logs/dedup_%A_%a.out
#SBATCH --error=logs/dedup_%A_%a.err

set -euo pipefail

module load gatk

SAMPLES=(
  "HL_Ctrl_CKDN240010648-1A_227JFTLT4_L6"
	"HL_AraC_CKDN240010649-1A_227JFTLT4_L6"
)

SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

BAM="../03_align/bams/${SAMPLE}.bam"
OUT="bams/${SAMPLE}.dedup.bam"
METRICS="metrics/${SAMPLE}.dedup.metrics.txt"

echo "BAM: ${BAM}"
echo "OUT: ${OUT}"
echo "METRICS: ${METRICS}"

mkdir -p bams metrics

gatk MarkDuplicates \
  -I ${BAM} \
  -O ${OUT} \
  -M ${METRICS} \
  --CREATE_INDEX true

# gatk MarkDuplicates \
#   -I bams/HL60_ctrl_1.merged.bam \
#   -O bams/HL60_ctrl_1.dedup.bam \
#   -M metrics/HL60_ctrl_1.dup_metrics.txt \
#   --CREATE_INDEX true

