#!/usr/bin/env bash

# Define 10kb intervals over the genome
bedtools makewindows -g ../hg19.chrom.sizes -w 10000 > hg19-10kb-window.bed

for cell in GM12878 IMR90 H1hESC; do
    # Sum meth count within 10kb intervals
    bedtools map -c 6 -o sum -F 0.5 \
        -a hg19-10kb-window.bed \
        -b <(awk -v OFS="\t" '{$1="chr"$1; print $0}' ../ASM/"${cell}"-ASM-All-*-hg19.bed) \
    | awk '{if($4==".") {$4="0"}} {print $0}' \
    | tr [:blank:] \\t \
    > "${cell}"-ASM-10kb-hg19.bedgraph
done
