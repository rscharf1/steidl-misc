#!/bin/bash
#SBATCH --job-name=Haplotype
#SBATCH --time=3:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=0-2
#SBATCH --output=logs/haplotype_%A_%a.out
#SBATCH --error=logs/haplotype_%A_%a.err

set -euo pipefail

module load gatk

SAMPLES=(
	"HL60_ctrl.fixedRG"
	"PY_high_1_CKDN260000485-1A_23KL3LLT4_L4"
	"PY_low_1_CKDN260000486-1A_23KL3LLT4_L4"
)

SAMPLE_NAMES=(
  "HL60_ctrl_1"
  "PY_high_1"
  "PY_low_1"
)

SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

REF=~/Tools/bwa_files/GRCh38.primary_assembly.genome.fa
BAM="../04_dedup/bams/${SAMPLE}.dedup.bam"
OUT="gvcfs/${SAMPLE}.g.vcf.gz"

echo "Ref: ${REF}"
echo "BAM: ${BAM}"
echo "OUT: ${OUT}"

mkdir -p gvcfs

gatk HaplotypeCaller \
  -R ${REF} \
  -I ${BAM} \
  -O ${OUT} \
  -ERC GVCF

