#!/usr/bin/env bash

# Define 10kb intervals over the genome
bedtools makewindows -g ../hg19.chrom.sizes -w 10000 > hg19-10kb-window.bed

zcat hg19.gc5Base.txt.gz \
| ./GCtoBedgraph.py \
| sort -k 1,1 -k2,2n \
> /tmp/hg19-GC.bedgraph #/tmp/GM12878-GC.bedgraph

# Map mean signal to 10kb windows - sum and divide by 2000 (because span = 5)
bedtools map -c 4 -o sum -F 0.5 \
    -a hg19-10kb-window.bed \
    -b /tmp/hg19-GC.bedgraph \
| awk '{if($4==".") {$4="0"}} {print $0}' \
| awk -v OFS="\t" '{print $1, $2, $3, ($4/2000)}' \
| tr [:blank:] \\t \
> GC-10kb-hg19.bedgraph
