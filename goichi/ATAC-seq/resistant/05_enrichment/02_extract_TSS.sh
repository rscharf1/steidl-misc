#!/bin/bash
#SBATCH --job-name=TSS
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00

# wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz

# awk 'BEGIN{OFS="\t"} 
#      $3=="chr1" || $3=="chr2" || $3 ~ /^chr/ {
#        if ($4=="+") 
#          print $3, $5, $5+1, $13, ".", $4;
#        else 
#          print $3, $6-1, $6, $13, ".", $4;
#      }' refGene.txt > hg38_TSS.bed


awk '$3=="gene"' ~/Tools/star_files/index/gencode.v44.annotation.gtf > gencode.v44.genes.gtf

awk 'BEGIN{OFS="\t"} 
{
    split($0,a,"gene_name \"");
    split(a[2],b,"\"");
    gene=b[1];

    if ($7 == "+")
        print $1, $4-1, $4, gene, ".", $7;
    else if ($7 == "-")
        print $1, $5-1, $5, gene, ".", $7;
}' gencode.v44.genes.gtf > hg38_TSS_gencode_v44.bed