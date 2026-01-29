#!/bin/bash

#SBATCH --job-name=basespace
#SBATCH --time=3:00:00              
#SBATCH --mem=64G                  
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --ntasks=1

bs download project -i 486397914 -o 012726
