#!/bin/bash

#SBATCH --job-name=STAR
#SBATCH --time=3:00:00              
#SBATCH --mem=64G                  
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --ntasks=1

module load STAR

STAR \
  --genomeDir ~/Tools/star_files/index/GRCh38_gencode_v44 \
  --readFilesIn sub_R1.fq sub_R2.fq \
  --runThreadN 8 \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix sub_