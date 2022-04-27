# Merge adjacent ASTADs and find intervals with both Imprinted and ASEG

for cell in GM12878 IMR90 H1hESC; do
    # Remove old copies of files
    rm "${cell}"-ASEG-inImprintedASTAD-GRCh37.bed

    # Merge ASTADs to single contiguous ASTAD domains
    bedtools merge \
        -i ../0.processDiffTAD/results/"${cell}"/alleleGRCh37/"${cell}"_a1-vs-"${cell}"_a2-all-20000-SNPsplit_diff_tad-ASTAD.bed \
        > /tmp/"${cell}"-ASTADs-merged.bed

    # Get ASTAD domains overlapping imprinted
    bedtools intersect -wa -u \
      -a /tmp/"${cell}"-ASTADs-merged.bed \
      -b Imprinted-GeneImprint-GRCh37.bed \
      > /tmp/"${cell}"-ASTADs-imprinted.bed


    # Get ASTAD domains overlapping mono-allelically expressed
    bedtools intersect -wa -u \
    -a /tmp/"${cell}"-ASTADs-merged.bed \
    -b "${cell}"-ASEG-*-GRCh37.bed  \
    > /tmp/"${cell}"-ASTADs-ASEG.bed

    # Get ASTAD domains similiar to Imprinted and Mono-allelic
    bedtools intersect -r -f 1 \
        -a /tmp/"${cell}"-ASTADs-ASEG.bed \
        -b /tmp/"${cell}"-ASTADs-imprinted.bed \
    > "${cell}"-ASTADdomains-Imprinted_and_ASEG.bed

    # Find all ASEG genes overlapping the relevant ASTADs
    bedtools intersect -wa \
        -a "${cell}"-ASEG-*-GRCh37.bed \
        -b "${cell}"-ASTADdomains-Imprinted_and_ASEG.bed \
    > /tmp/"${cell}"-ASEG-inImprintedASTAD-GRCh37.bed

    # Exclude Imprinted genes
    bedtools intersect -v \
        -a /tmp/"${cell}"-ASEG-inImprintedASTAD-GRCh37.bed \
        -b Imprinted-GeneImprint-GRCh37.bed \
    > "${cell}"-ASEG-inImprintedASTAD-GRCh37.bed

done
