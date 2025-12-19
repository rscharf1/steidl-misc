#!/bin/bash

#SBATCH --job-name=salmon
#SBATCH --time=6:00:00              
#SBATCH --mem=64G                  
#SBATCH --cpus-per-task=16

salmon alevin \
  --chromiumV3 \
  --index mm10_index \
  -l ISR \
  -1 ../scRNA/raw_data/H1/H1_CKDL210025541-1a-SI_TT_G6_HY3M7DSX2_S2_L002_R1_001.fastq.gz \
  -2 ../scRNA/raw_data/H1/H1_CKDL210025541-1a-SI_TT_G6_HY3M7DSX2_S2_L002_R2_001.fastq.gz \
  -p 8 \
  -o output/H1_out \
  --tgMap inputs/txp2gene_vM25_alevin.tsv \
  --whitelist inputs/H1_whitelist.txt \
  --dumpFeatures \
  --dumpBfh \
  --dumpMtx
