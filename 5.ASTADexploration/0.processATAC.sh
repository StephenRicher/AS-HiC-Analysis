#!/usr/bin/env bash

# Define 10kb intervals over the genome
bedtools makewindows -g ../hg19.chrom.sizes -w 10000 > hg19-10kb-window.bed

zcat GSE47753_GM12878_ATACseq_50k_AllReps_ZINBA_pp08.bed.gz \
| sed s'/chr//' \
| awk -v OFS="\t" '{print $1, $2, $3, $8}' \
| bedtools sort \
> /tmp/GM12878-ATAC.bedgraph

# Map mean signal to 10kb windows
bedtools map -c 4 -o sum -F 0.5 \
    -a hg19-10kb-window.bed \
    -b /tmp/GM12878-ATAC.bedgraph \
| awk '{if($4==".") {$4="0"}} {print $0}' \
| tr [:blank:] \\t \
> GM12878-ATAC-10kb-hg19.bedgraph
