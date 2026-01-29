#!/bin/sh

#SBATCH --partition=unlimited
#SBATCH --mem-per-cpu=100G
#SBATCH --nodes=8
#SBATCH --time=07-00
#SBATCH --job-name=SRA

fasterq-dump SRR33390394
