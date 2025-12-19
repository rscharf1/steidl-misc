#!/bin/bash

#SBATCH --job-name=salmon_index
#SBATCH --time=3:00:00              
#SBATCH --mem=64G                  
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --ntasks=1

salmon index --gencode -t gencode.vM25.transcripts.fa.gz -i mm10_index
