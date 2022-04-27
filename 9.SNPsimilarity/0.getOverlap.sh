declare -A pairs=( [GM12878]=H1hESC [H1hESC]=IMR90 [IMR90]=GM12878 )

for cell1 in "${!pairs[@]}"; do
    cell2=${pairs["${cell1}"]}
    echo "${cell1}" "${cell2}"

    # Get all regions that are BOTH ASTADs
    bedtools intersect \
        -a <(awk '$4=="ASTAD"' ../0.getDiffDomain/"${cell1}"-TADs-autosome.bed | bedtools sort) \
        -b <(awk '$4=="ASTAD"' ../0.getDiffDomain/"${cell2}"-TADs-autosome.bed | bedtools sort) \
    | bedtools sort \
    | bedtools merge \
    > "${cell1}"-"${cell2}"-ASTAD-overlap.bed

    # Get complement of those regions
    bedtools complement \
        -i "${cell1}"-"${cell2}"-ASTAD-overlap.bed \
        -g hg19.chrom.sizes \
    > /tmp/not-ASTAD-overlap.bed

    # Find regions that are both TADs (these may intesect the first list)
    bedtools intersect \
        -a <(bedtools sort -i ../0.getDiffDomain/"${cell1}"-TADs-autosome.bed) \
        -b <(bedtools sort -i ../0.getDiffDomain/"${cell2}"-TADs-autosome.bed) \
    | bedtools sort \
    | bedtools merge \
    > /tmp/"${cell1}"-"${cell2}"-TAD-overlap.bed

    # Filter to get only unique TAD intervals that don't overlap ASTADs
    bedtools intersect \
        -a /tmp/not-ASTAD-overlap.bed \
        -b /tmp/"${cell1}"-"${cell2}"-TAD-overlap.bed \
    | bedtools sort \
    | bedtools merge \
    > "${cell1}"-"${cell2}"-nonASTAD-overlap.bed

    # Perform isec with SNPs
    bcftools isec \
        -p /tmp/ASTADisec -c none --regions-file "${cell1}"-"${cell2}"-ASTAD-overlap.bed \
        ../0.processVCFs/"${cell1}"-hetOnly.vcf.gz ../0.processVCFs/"${cell2}"-hetOnly.vcf.gz

    bcftools isec \
        -p /tmp/nonASTADisec -c none --regions-file "${cell1}"-"${cell2}"-nonASTAD-overlap.bed \
        ../0.processVCFs/"${cell1}"-hetOnly.vcf.gz ../0.processVCFs/"${cell2}"-hetOnly.vcf.gz

    # Get total unique SNP positions with the ASTAD region
    cat \
        <(bcftools view -H --regions-file "${cell1}"-"${cell2}"-ASTAD-overlap.bed ../0.processVCFs/"${cell1}"-hetOnly.vcf.gz) \
        <(bcftools view -H --regions-file "${cell1}"-"${cell2}"-ASTAD-overlap.bed ../0.processVCFs/"${cell2}"-hetOnly.vcf.gz) \
    | cut -f 1-2 \
    | sort \
    | uniq \
    > /tmp/ASTAD-totalSNP.txt

    # Get total unique SNP positions with the ASTAD region
    cat \
        <(bcftools view -H --regions-file "${cell1}"-"${cell2}"-nonASTAD-overlap.bed ../0.processVCFs/"${cell1}"-hetOnly.vcf.gz) \
        <(bcftools view -H --regions-file "${cell1}"-"${cell2}"-nonASTAD-overlap.bed ../0.processVCFs/"${cell2}"-hetOnly.vcf.gz) \
    | cut -f 1-2 \
    | sort \
    | uniq \
    > /tmp/nonASTAD-totalSNP.txt

    nonASTADtotal=$(wc -l /tmp/nonASTAD-totalSNP.txt | awk '{print $1}')
    nonASTADidentical=$(bcftools view /tmp/nonASTADisec/0002.vcf | wc -l)
    nonASTADunique=$((nonASTADtotal-nonASTADidentical))
    ASTADtotal=$(wc -l /tmp/ASTAD-totalSNP.txt | awk '{print $1}')
    ASTADidentical=$(bcftools view /tmp/ASTADisec/0002.vcf | wc -l)
    ASTADunique=$((ASTADtotal-ASTADidentical))

    # Write cross-tabulation
    echo -e "\t"non-ASTAD"\t"ASTAD > "${cell1}"-"${cell2}"-SNP-overlapSummary.txt
    echo -e Identical"\t""${nonASTADidentical}""\t""${ASTADidentical}" >> "${cell1}"-"${cell2}"-SNP-overlapSummary.txt
    echo -e Unique"\t""${nonASTADunique}""\t""${ASTADunique}" >> "${cell1}"-"${cell2}"-SNP-overlapSummary.txt

done
