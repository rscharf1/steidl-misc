#!/bin/bash

awk -F',' '
/_Target_Genes/ {
    split($1,a,"_");
    tf=a[1];      # Extract TF name before "_Target_Genes"
    next
}
{
    for(i=1;i<=NF;i++){
        if($i!="") print tf "\t" $i
    }
}
' targets.csv > TF_gene_pairs.tsv