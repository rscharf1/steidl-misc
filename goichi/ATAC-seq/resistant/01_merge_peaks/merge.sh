#!/bin/bash
#SBATCH --job-name=merge
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=merge.out

cat peak_files/*.narrowPeak \
| sort -k1,1 -k2,2n \
| bedtools merge \
> merged_peaks.bed

awk '{
  mid = int(($2 + $3) / 2);
  start = mid - 250;
  end   = mid + 250;
  if (start < 0) start = 0;
  print $1, start, end
}' OFS="\t" merged_peaks.bed \
| sort -k1,1 -k2,2n \
| bedtools merge > merged_peaks_500bp.bed