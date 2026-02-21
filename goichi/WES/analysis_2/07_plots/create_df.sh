#!/bin/bash

# bcftools query \
#   -f '%CHROM\t%POS\t%AD\t%DP\n' \
#   ../06_filter/joint.exonic.dp20.vcf.gz > plot1/vaf.tsv

bcftools query \
  -f '%CHROM\t%POS\t[%AD\t]\t[%DP\t]\n' \
  ../06_filter/joint.exonic.dp20.vcf.gz > vaf.tsv
