#!/bin/bash

#SBATCH --job-name=salmon
#SBATCH --time=6:00:00              
#SBATCH --mem=64G                  
#SBATCH --cpus-per-task=16

salmon alevin \
  --chromiumV3 \
  --index mm10_index \
  -l ISR \
  -1 ../scRNA/raw_data/H2/H2_CKDL210025542-1a-SI_TT_G7_HY3M7DSX2_S1_L002_R1_001.fastq.gz \
  -2 ../scRNA/raw_data/H2/H2_CKDL210025542-1a-SI_TT_G7_HY3M7DSX2_S1_L002_R2_001.fastq.gz \
  -p 8 \
  -o output/H2_out \
  --tgMap inputs/txp2gene_vM25_alevin.tsv \
  --whitelist inputs/H2_whitelist.txt \
  --dumpFeatures \
  --dumpBfh \
  --dumpMtx
