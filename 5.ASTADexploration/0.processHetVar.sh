#!/usr/bin/env bash

# Define 10kb intervals over the genome
bedtools makewindows -g ../hg19.chrom.sizes -w 10000 > hg19-10kb-window.bed

for cell in GM12878 IMR90 H1hESC; do
    for var in snps indels; do
        if [ "${var}" = "snps" ]; then
            name=Het_SNPs
        else
            name=Het_Indels
        fi

        bcftools view -g het -v "${var}" ../0.processVCFs/"${cell}"-all.vcf.gz \
        | awk -v OFS="\t" '{print $1, $2, $2+1, 1}' \
        | bedtools sort \
        > /tmp/"${cell}"-hetOnly.bed

        # Sum SNP counts within 10kb intervals
        bedtools map -c 4 -o sum -F 0.5 \
            -a hg19-10kb-window.bed \
            -b /tmp/"${cell}"-hetOnly.bed \
            | awk '{if($4==".") {$4="0"}} {print $0}' \
            | tr [:blank:] \\t \
            > "${cell}"-"${name}"-10kb-hg19.bedgraph
    done
done
