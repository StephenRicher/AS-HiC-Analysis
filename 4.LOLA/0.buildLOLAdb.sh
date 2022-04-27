#!/usr/bin/env bash

mkdir -p data/LOLAregionDB/hg19/{GM12878,IMR90,H1hESC}/all/regions

for cell in GM12878 IMR90 H1hESC; do
    path=data/LOLAregionDB/hg19/"${cell}"/all
    index="${path}"/index.txt
    echo 'filename,cellType,description,dataSource' > "${index}"

    # Chrom HMM
    for chromHMM in ../2.chromHMM/"${cell}"-chromHMM-*.bed; do
        cp "${chromHMM}" "${path}"/regions/"${chromHMM##*/}"
    done

    echo "${cell}"'-chromHMM-ActiveTSS.bed,'"${cell}"',Active TSS,NA' >> "${index}"
    echo "${cell}"'-chromHMM-FlankingActiveTSS.bed,'"${cell}"',Flanking Active TSS,NA' >> "${index}"
    echo "${cell}"'-chromHMM-TranscriptionAt5_3.bed,'"${cell}"',Transcription At 5_3,NA' >> "${index}"
    echo "${cell}"'-chromHMM-StrongTranscription.bed,'"${cell}"',Strong Transcription,NA' >> "${index}"
    echo "${cell}"'-chromHMM-WeakTranscription.bed,'"${cell}"',Weak Transcription,NA' >> "${index}"
    echo "${cell}"'-chromHMM-GenicEnhancer.bed,'"${cell}"',Genic Enhancer,NA' >> "${index}"
    echo "${cell}"'-chromHMM-Enhancers.bed,'"${cell}"',Enhancers,NA' >> "${index}"
    echo "${cell}"'-chromHMM-ZNFgenesRepeats.bed,'"${cell}"',ZNF / Genes Repeats,NA' >> "${index}"
    echo "${cell}"'-chromHMM-Heterochromatin.bed,'"${cell}"',Heterochromatin,NA' >> "${index}"
    echo "${cell}"'-chromHMM-BivalentPoisedTSS.bed,'"${cell}"',Bivalent Poised TSS,NA' >> "${index}"
    echo "${cell}"'-chromHMM-FlankingBivalentTSS.bed,'"${cell}"',Flanking Bivalent TSS,NA' >> "${index}"
    echo "${cell}"'-chromHMM-BivalentEnhancer.bed,'"${cell}"',Bivalent Enhancer,NA' >> "${index}"
    echo "${cell}"'-chromHMM-RepressedPolycomb.bed,'"${cell}"',Repressed Polycomb,NA' >> "${index}"
    echo "${cell}"'-chromHMM-WeakRepressedPolycomb.bed,'"${cell}"',Weak Repressed Polycomb,NA' >> "${index}"
    echo "${cell}"'-chromHMM-Quiescent.bed,'"${cell}"',Quiescent,NA' >> "${index}"

    # A compartment
    cat  ../../"${cell}"/fullGRCh37/dat/Cscore/*/10000/"${cell}"-chr*-10000-full-Cscore-reorient.bedgraph \
    | awk '$4>0' | bedtools sort \
    > "${path}"/regions/"${cell}"-Cscore-A.bedgraph
    echo "${cell}"'-Cscore-A.bedgraph,'"${cell}"',A compartment,NA'>> "${index}"

    # B compartment
    cat  ../../"${cell}"/fullGRCh37/dat/Cscore/*/10000/"${cell}"-chr*-10000-full-Cscore-reorient.bedgraph \
    | awk '$4<0' | bedtools sort \
    > "${path}"/regions/"${cell}"-Cscore-B.bedgraph
    echo "${cell}"'-Cscore-B.bedgraph,'"${cell}"',B compartment,NA' \
    >> "${index}"

    # lncRNA
    grep 'gene_type=lncRNA' ../../annotation/gencode.v38lift37-genes.gff3 \
    | awk -v OFS="\t" '{print "chr"$1, $4, $5}' \
    > "${path}"/regions/gencode.v38lift37.lncRNA.bed
    echo 'gencode.v38lift37.lncRNA.bed,'"${cell}"',lncRNA,NA' >> "${index}"

    # Protein Coding
    grep 'gene_type=protein_coding' ../../annotation/gencode.v38lift37-genes.gff3 \
    | awk -v OFS="\t" '{print "chr"$1, $4, $5}' \
    > "${path}"/regions/gencode.v38lift37.proteinCoding.bed
    echo 'gencode.v38lift37.proteinCoding.bed,'"${cell}"',Protein coding,NA' \
    >> "${index}"

    #for t in snps indels; do
    #    bcftools view -v "${t}" -g het ../0.processVCFs/"${cell}"-all.vcf.gz \
    #    | grep -v '#' \
    #    | awk -v OFS="\t" '{print $1, $2, $2+1}' \
    #    | bedtools sort > "${path}"/regions/"${cell}"-Het-"${t}".bed
    #    echo "${cell}"-Het-"${t}".bed,"${cell}",Het "${t}",NA >> "${index}"
    #done

    tail -n +2 ../../"${cell}"/"${cell}"-CNV-UCSC.bed \
    | awk -v OFS="\t" '$4=="Loss"' | cut -f 1-3 \
    > "${path}"/regions/"${cell}"-Loss.bed
    echo "${cell}"-Loss.bed,"${cell}",Loss,NA >> "${index}"

    tail -n +2 ../../"${cell}"/"${cell}"-CNV-UCSC.bed \
    | awk -v OFS="\t" '$4=="Deletion"' | cut -f 1-3 \
    > "${path}"/regions/"${cell}"-Deletion.bed
    echo "${cell}"-Deletion.bed,"${cell}",Deletion,NA >> "${index}"

    tail -n +2 ../../"${cell}"/"${cell}"-CNV-UCSC.bed \
    | awk -v OFS="\t" '$4=="Gain" || $4=="Amplification"' | cut -f 1-3 \
    > "${path}"/regions/"${cell}"-Gain.bed
    echo "${cell}"-Gain.bed,"${cell}",Gain,NA >> "${index}"

    tail -n +2 ../../"${cell}"/"${cell}"-CNV-UCSC.bed \
    | awk -v OFS="\t" '$4=="Normal"' | cut -f 1-3  \
    > "${path}"/regions/"${cell}"-Normal.bed
    echo "${cell}"-Normal.bed,"${cell}",Normal,NA >> "${index}"

    cp ../../../HiC-Data/annotation/GRCh37-CpG.bed "${path}"/regions/"${cell}"-CpG.bed
    echo "${cell}"-CpG.bed,"${cell}",CpG,NA >> "${index}"

done
