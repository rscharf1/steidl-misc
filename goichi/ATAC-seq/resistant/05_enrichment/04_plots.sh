#!/bin/bash
#SBATCH --job-name=plots
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00

mkdir -p plots

# plotHeatmap \
#   -m TSS_matrix.gz \
#   -out plots/TSS_heatmap1.pdf \
#   --colorMap Reds \
#   --whatToShow 'heatmap and colorbar' \
#   --refPointLabel "TSS" \
#   --regionsLabel "Genes" \
#   --samplesLabel S74 S75 S76 S77 S28 S29

# plotHeatmap \
#   -m TSS_matrix.gz \
#   -out plots/TSS_heatmap2.pdf \
#   --sortRegions descend \
#   --sortUsing mean \
#   --colorMap Reds \
#   --samplesLabel S74 S75 S76 S77 S28 S29  

# plotHeatmap \
#   -m TSS_matrix.gz \
#   -out plots/TSS_heatmap3.pdf \
#   --sortRegions descend \
#   --sortUsing mean \
#   --colorMap Reds \
#   --missingDataColor white \
#   --dpi 300 \
#   

###########

plotProfile \
  -m TSS_matrix.gz \
  -out plots/TSS_profile1.pdf \
  --perGroup \
  --refPointLabel "TSS" \
  --samplesLabel S74 S75 S76 S77 S28 S29

plotProfile \
  -m TSS_matrix.gz \
  -out plots/TSS_profile2.pdf \
  --perGroup \
  --plotHeight 5 \
  --plotWidth 6 \
  --dpi 300 \
  --samplesLabel S74 S75 S76 S77 S28 S29


  