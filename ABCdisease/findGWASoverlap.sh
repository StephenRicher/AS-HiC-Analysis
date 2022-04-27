
for cell in GM12878 IMR90 H1hESC; do
    mkdir -p "${cell}"

    # Only consider Het Variants overlapping TADs
    bedtools intersect -wa -u \
        -a <(bcftools view -g het ../0.processVCFs/"${cell}"-hetOnly.vcf.gz) \
        -b <(bedtools merge -i ../0.getDiffDomain/"${cell}"-TADs-autosome.bed) \
    | awk -v OFS="\t" '{print $1, $2, $2+1}' \
    > /tmp/"${cell}"-Het-overlapTAD.bed

    # Only consider GWAS variants Heterozygous in cell line AND overlapping TADs
    bedtools intersect -wa -u \
        -a <(zcat ../annotation/gwas_catalog_v1.0-associations_e104_r2021-11-22.hg19.bed.gz) \
        -b /tmp/"${cell}"-Het-overlapTAD.bed \
    > /tmp/"${cell}"-HetGWAS.bed

    if [ "${cell}" = H1hESC ]; then
        name=H1-hESC
    else
        name="${cell}"
    fi

    # Filter ABC by cell type and those overlapping TADs
    zcat AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz \
    | grep "${name}" \
    | bedtools sort \
    | bedtools intersect -wa -u -a - -b <(bedtools merge -i ../0.getDiffDomain/"${cell}"-TADs-autosome.bed) \
    > /tmp/"${cell}"-ABC-predictions.bed

    # Save all valid ABC enhancers overlapping at ASTAD
    bedtools intersect -wa -u \
        -a /tmp/"${cell}"-ABC-predictions.bed \
        -b ../0.getDiffDomain/"${cell}"-ASTADs-autosome.bed \
    > "${cell}"/"${cell}"-ABC_in_ASTAD.bed

    # Save all other ABC enhancers NOT in ASTAD
    comm -32 \
        <(sort /tmp/"${cell}"-ABC-predictions.bed) \
        <(sort "${cell}"/"${cell}"-ABC_in_ASTAD.bed) \
    > "${cell}"/"${cell}"-ABC_Not_in_ASTAD.bed

    # Filter ABC enhancers to include only those with Het GWAS variant
    bedtools intersect -wa -u \
        -a /tmp/"${cell}"-ABC-predictions.bed \
        -b /tmp/"${cell}"-HetGWAS.bed \
    > /tmp/"${cell}"-ABC_with_GWAS.bed

    # Get Het GWAS variants overlapping an ABC enhancer
    bedtools intersect -wa -u \
    	-a /tmp/"${cell}"-HetGWAS.bed \
    	-b  /tmp/"${cell}"-ABC_with_GWAS.bed \
    > /tmp/"${cell}"-GWAS_with_ABC.bed

    # Get all Het variants overlapping ABC enhancer
    bedtools intersect -wa -u \
    	-a /tmp/"${cell}"-Het-overlapTAD.bed \
    	-b /tmp/"${cell}"-ABC-predictions.bed \
    > /tmp/"${cell}"-Het-overlapABC.bed

    # Get all other Het variant that DO NOT overlap an ABC GWAS
    comm -32 \
    	<(sort /tmp/"${cell}"-Het-overlapTAD.bed) \
    	<(sort /tmp/"${cell}"-Het-overlapABC.bed) \
    > /tmp/"${cell}"-Het-NotoverlapABC.bed

    # Get all Het variants overlapping ABC GWAS in ASTAD
    bedtools intersect -wa -u \
    	-a /tmp/"${cell}"-Het-overlapABC.bed \
    	-b ../0.getDiffDomain/"${cell}"-ASTADs-autosome.bed \
    > "${cell}"/"${cell}"-Het-overlapABC-ASTAD.bed

    # Get all Het variants overllaping ABC NOT in ASTAD
    comm -32 \
    	<(sort /tmp/"${cell}"-Het-overlapABC.bed) \
    	<(sort "${cell}"/"${cell}"-Het-overlapABC-ASTAD.bed) \
    > "${cell}"/"${cell}"-Het-overlapABC-nonASTAD.bed

    # Get all Het variants NOT overlapping ABC GWAS in ASTAD
    bedtools intersect -wa -u \
    	-a /tmp/"${cell}"-Het-NotoverlapABC.bed \
    	-b ../0.getDiffDomain/"${cell}"-ASTADs-autosome.bed \
    > "${cell}"/"${cell}"-Het-NotoverlapABC-ASTAD.bed

    # Get all Het variants overllaping ABC NOT in ASTAD
    comm -32 \
    	<(sort /tmp/"${cell}"-Het-NotoverlapABC.bed) \
    	<(sort "${cell}"/"${cell}"-Het-NotoverlapABC-ASTAD.bed) \
    > "${cell}"/"${cell}"-Het-NotoverlapABC-nonASTAD.bed

    for group in ASTAD nonASTAD; do
        # Retrieve the GWAS IDs for the Het SNPs overlapping ABC in ASTADs / non-ASTADs
        bedtools intersect -wa -u \
            -a /tmp/"${cell}"-GWAS_with_ABC.bed \
            -b "${cell}"/"${cell}"-Het-overlapABC-"${group}".bed \
        > /tmp/"${cell}"-GWASids-"${group}".bed

        # Final list of all ABC enhancers that contain a Het GWAS variant in an ASTAD
        bedtools intersect -wa -wb \
            -a /tmp/"${cell}"-ABC_with_GWAS.bed \
            -b /tmp/"${cell}"-GWASids-"${group}".bed \
        | cut -f 1-24,28 \
        | awk -v OFS="\t" -v group="${group}" '{print $0, group}'
    done > "${cell}"/"${cell}"-ABC-GWAS-allStats.bed

    cat \
    	<(awk '{print "ASTAD-ABC"}' "${cell}"/"${cell}"-Het-overlapABC-ASTAD.bed) \
    	<(awk '{print "ASTAD-notABC"}' "${cell}"/"${cell}"-Het-NotoverlapABC-ASTAD.bed) \
    	<(awk '{print "notASTAD-notABC"}' "${cell}"/"${cell}"-Het-NotoverlapABC-nonASTAD.bed) \
    	<(awk '{print "notASTAD-ABC"}' "${cell}"/"${cell}"-Het-overlapABC-nonASTAD.bed) \
    | uniq -c \
    | awk -v cell="${cell}" '{print cell, $1, $2}'
done > GWAS-ABC-overlaps.txt
