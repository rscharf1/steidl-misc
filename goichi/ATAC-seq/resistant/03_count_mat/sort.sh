bedtools sort \
  -i ../02_blacklist/merged_peaks_500bp.filtered.bed \
  -g hg38.genome \
> merged_peaks_500bp.filtered.sorted.bed
