#!/bin/bash
#SBATCH --job-name=rclone
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00

rclone copy \
"dropbox:Goichi Sequencing Data/WES HL60/" \
fastqs \
--include "*.fq.gz" \
--progress \
--transfers 4 \
--checkers 8 \
--multi-thread-streams 4

# rclone copy \
# "dropbox:Goichi Sequencing Data/WES HL60/" \
# WES_HL60 \
# --include "*.fq.gz" \
# --dry-run