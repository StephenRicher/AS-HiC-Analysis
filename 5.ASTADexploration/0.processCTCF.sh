#!/usr/bin/env bash

# Define 10kb intervals over the genome
bedtools makewindows -g ../hg19.chrom.sizes -w 10000 > hg19-10kb-window.bed

for cell in GM12878 IMR90 H1hESC; do
    if [ "${cell}" = "GM12878" ]; then
        CTCF=../../annotation/GSM749704_hg19_wgEncodeUwTfbsGm12878CtcfStdRawRep1-GM12878.bigWig
    elif [ "${cell}" = "IMR90" ]; then
        CTCF=../../annotation/GSM935404_hg19_wgEncodeSydhTfbsImr90CtcfbIggrabSig-IMR90.bigWig
    else
        CTCF=../../annotation/wgEncodeBroadHistoneH1hescCtcfStdSig-H1hESC.bigWig
    fi

    # Map mean signal to 10kb windows
    bedtools map -c 4 -o sum -F 0.5 \
        -a hg19-10kb-window.bed \
        -b <(bigWigToBedGraph "${CTCF}" /dev/stdout) \
    | awk '{if($4==".") {$4=0} {print $0}}' \
    | tr [:blank:] \\t \
    > "${cell}"-CTCF-10kb-hg19.bedgraph
done
