#!/bin/bash

#SBATCH --job-name=seqtk
#SBATCH --time=3:00:00              
#SBATCH --mem=64G                  
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --ntasks=1

module load seqtk

seqtk sample -s100 ../02_trim/trimmed/HL60_ctrl_1_CKDN260000487-1A_2335HJLT4_L1_R1.trimmed.fastq.gz 1000000 > sub_R1.fq
seqtk sample -s100 ../02_trim/trimmed/HL60_ctrl_1_CKDN260000487-1A_2335HJLT4_L1_R2.trimmed.fastq.gz 1000000 > sub_R2.fq