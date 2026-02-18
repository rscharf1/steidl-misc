#!/bin/bash
#SBATCH --job-name=groups
#SBATCH --time=3:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --ntasks=1

module load gatk

gatk AddOrReplaceReadGroups \
  -I ../04_dedup/bams/HL60_ctrl.merged.dedup.bam \
  -O ../04_dedup/bams/HL60_ctrl.fixedRG.dedup.bam \
  -RGID HL60_ctrl_1 \
  -RGLB HL60_ctrl_1 \
  -RGPL ILLUMINA \
  -RGPU merged \
  -RGSM HL60_ctrl_1 \
  --CREATE_INDEX true