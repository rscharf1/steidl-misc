#!/bin/bash

#SBATCH --job-name=MultiQC
#SBATCH --time=3:00:00              
#SBATCH --mem=64G                  
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --ntasks=1

multiqc .
