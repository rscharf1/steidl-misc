#!/bin/bash
#SBATCH --job-name=WES_align
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=0-1
#SBATCH --output=logs/WES_align_%A_%a.out
#SBATCH --error=logs/WES_align_%A_%a.err

set -euo pipefail

REF=~/Tools/bwa_files/GRCh38.primary_assembly.genome.fa
FASTQ_DIR="$PWD/../02_trim/trimmed"
OUT="$PWD/bams"

SAMPLES=(
	"HL_Ctrl_CKDN240010648-1A_227JFTLT4_L6"
	"HL_AraC_CKDN240010649-1A_227JFTLT4_L6"
)

SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

R1="${FASTQ_DIR}/${SAMPLE}_R1.trimmed.fastq.gz"
R2="${FASTQ_DIR}/${SAMPLE}_R2.trimmed.fastq.gz"
OUT_BAM="${OUT}/${SAMPLE}.bam"

echo "Running BWA_MEM for sample: ${SAMPLE}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"

echo "R1: ${R1}"
echo "R2: ${R2}"
echo "Output BAM: ${OUT_BAM}"

module load bwa

mkdir -p bams

bwa mem -t 16 \
  -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA" \
  "$REF" "$R1" "$R2" | \
samtools sort -@ 8 -o "$OUT_BAM"

samtools index "$OUT_BAM"



