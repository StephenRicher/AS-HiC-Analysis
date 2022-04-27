
echo -e cell"\t"tad"\t"variant"\t"genotype"\t"ASTAD"\t"total \
> variant-enrichment.txt
for cell in GM12878 IMR90 H1hESC; do
    for tad in TADs ASTADs; do
        bedtools merge -i ../0.getDiffDomain/"${cell}"-"${tad}"-autosome.bed \
        | sed s'/chr//' \
        > /tmp/"${cell}"-"${tad}"-merged.bed
    done
    for var in snps indels; do
        for g in het hom; do
            bcftools view -g "${g}" -v "${var}" \
                ../0.processVCFs/"${cell}"-all.vcf.gz \
            | sed s'/chr//' \
            > /tmp/"${cell}"-"${g}"-"${v}".bed

            # Get total number of variant overlapping tads
            bedtools intersect -wa -u \
                -a /tmp/"${cell}"-"${g}"-"${v}".bed \
                -b /tmp/"${cell}"-TADs-merged.bed \
            > /tmp/"${cell}"-"${g}"-"${v}"-all.bed
            total=$(wc -l /tmp/"${cell}"-"${g}"-"${v}"-all.bed | awk '{print $1}')

            # Count number of variants that overlap one of these domains
            bedtools intersect -wa -u \
                -a /tmp/"${cell}"-"${g}"-"${v}".bed \
                -b /tmp/"${cell}"-ASTADs-merged.bed \
            > /tmp/"${cell}"-"${g}"-"${v}"-ASTADs-overlap.bed
            overlap=$(wc -l /tmp/"${cell}"-"${g}"-"${v}"-ASTADs-overlap.bed | awk '{print $1}')

            echo -e "${cell}""\t""${tad}""\t""${var}""\t""${g}""\t""${overlap}""\t""${total}" \
            >> variant-enrichment.txt
        done
    done
    rm /tmp/"${cell}"*.bed
done
