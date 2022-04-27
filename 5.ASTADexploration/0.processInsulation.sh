#!/usr/bin/env bash

# Define 10kb intervals over the genome
bedtools makewindows -g ../hg19.chrom.sizes -w 10000 > hg19-10kb-window.bed

for cell in GM12878 IMR90 H1hESC; do
    cat ../../"${cell}"/fullGRCh37/dat/tads/chr*/10000/"${cell}"-chr*-10000-full_tad_score.bm \
    | grep -v '#' \
    | awk -v OFS="\t" '{score=($4+$5+$6+$7)/4; print $1, $2, $3, score}' \
    | bedtools sort \
    > "${cell}"-Insulation.bedgraph

    # Map mean signal to 10kb windows
    bedtools map -c 4 -o sum -F 0.5 \
        -a hg19-10kb-window.bed \
        -b "${cell}"-Insulation.bedgraph \
    | awk '{if($4==".") {$4=0} {print $0}}' \
    | tr [:blank:] \\t \
    > "${cell}"-Insulation-10kb-hg19.bedgraph
done
