#!/bin/bash

#SBATCH --job-name=BWA_index
#SBATCH --time=3:00:00              
#SBATCH --mem=64G                  
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --ntasks=1

# bwa index ~/Tools/star_files/index/GRCh38.primary_assembly.genome.fa

samtools faidx ~/Tools/star_files/index/GRCh38.primary_assembly.genome.fa
