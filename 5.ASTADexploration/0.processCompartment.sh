#!/usr/bin/env bash

# Define 10kb intervals over the genome
bedtools makewindows -g ../hg19.chrom.sizes -w 10000 > hg19-10kb-window.bed

for cell in GM12878 IMR90 H1hESC; do

    cat ../../"${cell}"/fullGRCh37/dat/Cscore/chr*/10000/"${cell}"-chr*-10000-full-Cscore-reorient.bedgraph \
    | awk -v OFS="\t" '{print $1, $2, $3, ".", $4}' \
    | bedtools sort \
    > /tmp/"${cell}"-Cscore.bed

    # Map mean signal to 10kb windows
    bedtools map -c 5 -o sum -F 0.5 \
        -a hg19-10kb-window.bed \
        -b /tmp/"${cell}"-Cscore.bed \
    | awk '{if($4==".") {$4=0} {print $0}}' \
    | tr [:blank:] \\t \
    > "${cell}"-Cscore-10kb-hg19.bedgraph

done
