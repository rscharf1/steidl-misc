#!/bin/bash
#SBATCH --job-name=TSS
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00

computeMatrix reference-point \
  -S bw/* \
  -R hg38_TSS_gencode_v44.bed \
  --referencePoint TSS \
  -b 1000 -a 1000 \
  --binSize 10 \
  -p 8 \
  -o TSS_matrix.gz