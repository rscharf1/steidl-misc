#!/bin/bash
#SBATCH --job-name=homer
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=homer.out

# conda activate homer

mkdir -p homer_out

# annotatePeaks.pl \
#   ../03_count_mat/merged_peaks_500bp.filtered.sorted.bed \
#   ~/Tools/homer/data/genomes/hg38 \
#   -m RUNX1.motif \
# > homer_out/peaks_RUNX1_annotated.txt

annotatePeaks.pl \
  ../03_count_mat/merged_peaks_500bp.filtered.sorted.bed \
  ~/Tools/homer/data/genomes/hg38 \
> peaks_annotated.txt