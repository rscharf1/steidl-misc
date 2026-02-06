#!/bin/bash
#SBATCH --job-name=ATAC_filter
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=16:00:00
#SBATCH --output=logs3/filter_%A_%a.out
#SBATCH --array=1-10

set -euo pipefail

module load samtools
module load picard

INPUT=$(sed -n "${SLURM_ARRAY_TASK_ID}p" input_bams.txt)
SAMPLE=$(basename "$INPUT" .sorted.bam)

echo $INPUT
echo $SAMPLE

# mkdir -p q30

# # # Only high confidence 
# samtools view -b -q 30 ${INPUT} > q30/${SAMPLE}.q30.bam
# samtools index q30/${SAMPLE}.q30.bam

# # # Remove mito reads 
# mkdir -p noMT

# samtools idxstats q30/${SAMPLE}.q30.bam \
#   | cut -f1 \
#   | grep -v chrM > noMT/${SAMPLE}.chroms.txt

# samtools view -b q30/${SAMPLE}.q30.bam $(cat noMT/${SAMPLE}.chroms.txt) > noMT/${SAMPLE}.noMT.bam
# samtools index noMT/${SAMPLE}.noMT.bam

# # # Remove dups
# mkdir -p noDups

# picard MarkDuplicates \
#   I=noMT/${SAMPLE}.noMT.bam \
#   O=noDups/${SAMPLE}.dedup.bam \
#   M=noDups/${SAMPLE}.dup_metrics.txt \
#   REMOVE_DUPLICATES=true

# samtools index noDups/${SAMPLE}.dedup.bam

# # Tn5 shifting
mkdir -p Tn5_shift

# alignmentSieve \
#   -b noDups/${SAMPLE}.dedup.bam \
#   --ATACshift \
#   -o Tn5_shift/${SAMPLE}.shifted.bam

# samtools view -h noDups/${SAMPLE}.dedup.bam \
# | gawk 'BEGIN{OFS="\t"}
#        /^@/ {print; next}
#        {
#          if (and($2,16)) {$4 = $4 - 5}
#          else {$4 = $4 + 4}
#          print
#        }' \
# | samtools view -b -o Tn5_shift/${SAMPLE}.shifted.bam

samtools view -h noDups/${SAMPLE}.dedup.bam \
| gawk 'BEGIN{OFS="\t"}
       /^@/ {print; next}
       {
         pos = $4
         if (and($2,16)) {
           pos = pos - 5
         } else {
           pos = pos + 4
         }

         if (pos >= 1) {
           $4 = pos
           print
         }
       }' \
| samtools view -b -o Tn5_shift/${SAMPLE}.shifted.bam

samtools sort -@ 4 \
  -o Tn5_shift/${SAMPLE}.shifted.sorted.bam \
  Tn5_shift/${SAMPLE}.shifted.bam

samtools index Tn5_shift/${SAMPLE}.shifted.sorted.bam
