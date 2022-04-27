#!/usr/bin/env bash

# Process WGBS sequencing data to mean signal per 10kb in hg19.
# Convert bigwig to bedgraph.
# Map sum signal to 10kb. Sum is appropriate because window size is constant.
# Remap to GRCh37. - ideally we would remap first but leads to memory error.
# Then take mean signal across hg19 windows (to fix any problems during remapping)

# Define 10kb intervals over the genome
bedtools makewindows -g ../hg38.chrom.sizes -w 10000 > hg38-10kb-window.bed

# Map mean signal to 10kb windows
bedtools map -c 4 -o sum -F 0.5 \
    -a hg38-10kb-window.bed \
    -b <(bigWigToBedGraph GSM2308632_ENCFF796NFQ_signal_GRCh38.bigWig /dev/stdout | sed 's/chr//') \
| awk '{if($4==".") {$4=0} {print $0}}' \
> GM12878-WGBS-10kb-GRCh38.bedgraph

# Remap to hg19 coordinates
CrossMap.py bed ../hg38ToHg19.over.chain.gz \
    GM12878-WGBS-10kb-GRCh38.bedgraph \
    GM12878-WGBS-10kb-hg19-remap.bedgraph

# Sort
bedtools sort -i GM12878-WGBS-10kb-hg19-remap.bedgraph \
> GM12878-WGBS-10kb-hg19-remap.sorted.bedgraph

# Repeat window mapping to new build
bedtools makewindows -g ../hg19.chrom.sizes -w 10000 > hg19-10kb-window.bed

bedtools map -c 4 -o mean -F 0.5 \
    -a hg19-10kb-window.bed \
    -b GM12878-WGBS-10kb-hg19-remap.sorted.bedgraph \
| awk '{if($4==".") {$4="0"}} {print $0}' \
| tr [:blank:] \\t \
> GM12878-WGBS-10kb-hg19.bedgraph
