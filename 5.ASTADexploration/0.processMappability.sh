#!/usr/bin/env bash

# Define 10kb intervals over the genome
bedtools makewindows -g ../hg19.chrom.sizes -w 10000 > hg19-10kb-window.bed

mappability=../../annotation/wgEncodeCrgMapabilityAlign100mer.bigWig

# Map mean signal to 10kb windows
bedtools map -c 4 -o sum -F 0.5 \
    -a hg19-10kb-window.bed \
    -b <(bigWigToBedGraph "${mappability}" /dev/stdout | sed 's/chr//') \
| awk '{if($4==".") {$4=0} {print $0}}' \
| tr [:blank:] \\t \
> Mappability-10kb-hg19.bedgraph
