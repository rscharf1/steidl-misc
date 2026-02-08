#!/bin/bash
#SBATCH --job-name=filter
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=filter.out

bedtools intersect \
  -v \
  -a ../01_merge_peaks/merged_peaks_500bp.bed \
  -b hg38-blacklist.v2.bed > merged_peaks_500bp.filtered.bed