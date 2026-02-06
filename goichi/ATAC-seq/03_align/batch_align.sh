#!/bin/bash
#SBATCH --job-name=bowtie_align
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=016:00:00
#SBATCH --output=logs/bowtie_%A_%a.out
#SBATCH --array=1-2

set -euo pipefail

module load bowtie2

R1=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples_rerun.txt)
R2=${R1/_R1.trimmed.fastq.gz/_R2.trimmed.fastq.gz}
SAMPLE=$(basename "$R1" _R1.trimmed.fastq.gz)

echo $R1
echo $R2
echo $SAMPLE

mkdir -p output

bowtie2 \
  -x ~/Tools/bowtie_files/GRCh38 \
  -1 ${R1} \
  -2 ${R2} \
  --very-sensitive \
  -X 2000 \
  --no-mixed \
  --no-discordant \
  -p 8 \
  2> output/${SAMPLE}.bowtie2.log \
| samtools view -bS - > output/${SAMPLE}.raw.bam

samtools sort -@ 8 -o output/${SAMPLE}.sorted.bam output/${SAMPLE}.raw.bam
samtools index output/${SAMPLE}.sorted.bam
