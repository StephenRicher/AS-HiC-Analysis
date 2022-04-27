cell1=GM12878
cell2=H1hESC
cell3=IMR90

diffTAD1=../0.getDiffDomain/"${cell1}"-TADs-autosome.bed
diffTAD2=../0.getDiffDomain/"${cell2}"-TADs-autosome.bed
diffTAD3=../0.getDiffDomain/"${cell3}"-TADs-autosome.bed

# Get all regions that are BOTH ASTADs
bedtools multiinter -i \
    <(awk '$4 == "ASTAD"' "${diffTAD1}") \
    <(awk '$4 == "ASTAD"' "${diffTAD2}") \
    <(awk '$4 == "ASTAD"' "${diffTAD3}") \
| awk '$4==3' | awk -v OFS="\t" '{print $1, $2, $3}' \
| bedtools sort \
| bedtools merge \
> 3way-ASTAD-overlap.bed

# Get complement of those regions
bedtools complement \
    -i 3way-ASTAD-overlap.bed \
    -g ../hg19.chrom.sizes \
> /tmp/not-3way-ASTAD-overlap.bed

# Find regions that are both non-ASTADs (these may intesect the first list)
bedtools multiinter -i \
    "${diffTAD1}" \
    "${diffTAD2}" \
    "${diffTAD3}" \
| awk '$4==3' | awk -v OFS="\t" '{print $1, $2, $3}' \
| bedtools sort \
| bedtools merge \
> /tmp/3way-TAD-overlap.bed

# Filter to get only unique TAD intervals that don't overlap ASTADs
bedtools intersect \
    -a /tmp/not-3way-ASTAD-overlap.bed \
    -b /tmp/3way-TAD-overlap.bed \
| bedtools sort \
| bedtools merge \
> 3way-nonASTAD-overlap.bed

# Perform isec with SNPs
bcftools isec -n=3 \
    -p /tmp/ASTADisec -c none --regions-file 3way-ASTAD-overlap.bed \
    ../0.processVCFs/"${cell1}"-hetOnly.vcf.gz ../0.processVCFs/"${cell2}"-hetOnly.vcf.gz \
    ../0.processVCFs/"${cell3}"-hetOnly.vcf.gz

# Get total unique SNP positions within the ASTAD overlap region
cat \
    <(bcftools view -H --regions-file 3way-ASTAD-overlap.bed ../0.processVCFs/"${cell1}"-hetOnly.vcf.gz) \
    <(bcftools view -H --regions-file 3way-ASTAD-overlap.bed ../0.processVCFs/"${cell2}"-hetOnly.vcf.gz) \
    <(bcftools view -H --regions-file 3way-ASTAD-overlap.bed ../0.processVCFs/"${cell3}"-hetOnly.vcf.gz) \
| cut -f 1-2 \
| sort \
| uniq \
> /tmp/ASTAD-totalSNP.txt

ASTADtotal=$(wc -l /tmp/ASTAD-totalSNP.txt | awk '{print $1}')
ASTADidentical=$(wc -l /tmp/ASTADisec/sites.txt | awk '{print $1}')
ASTADunique=$((ASTADtotal-ASTADidentical))

# Find 3way intervals that DO NOT overlap any identical SNP
bedtools intersect -v \
    -a 3way-ASTAD-overlap.bed \
    -b <(awk -v OFS="\t" '{print $1, $2, $2+1}' /tmp/ASTADisec/sites.txt) \
> 3way-ASTAD-noIdentical.bed

# Repeat for non-ASTAD regions
bcftools isec -n=3 \
    -p /tmp/nonASTADisec -c none --regions-file 3way-nonASTAD-overlap.bed \
    ../0.processVCFs/"${cell1}"-hetOnly.vcf.gz ../0.processVCFs/"${cell2}"-hetOnly.vcf.gz \
    ../0.processVCFs/"${cell3}"-hetOnly.vcf.gz

# Get total unique SNP positions with the region
cat \
    <(bcftools view -H --regions-file 3way-nonASTAD-overlap.bed ../0.processVCFs/"${cell1}"-hetOnly.vcf.gz) \
    <(bcftools view -H --regions-file 3way-nonASTAD-overlap.bed ../0.processVCFs/"${cell2}"-hetOnly.vcf.gz) \
    <(bcftools view -H --regions-file 3way-nonASTAD-overlap.bed ../0.processVCFs/"${cell3}"-hetOnly.vcf.gz) \
| cut -f 1-2 \
| sort \
| uniq \
> /tmp/nonASTAD-totalSNP.txt

nonASTADtotal=$(wc -l /tmp/nonASTAD-totalSNP.txt | awk '{print $1}')
nonASTADidentical=$(wc -l /tmp/nonASTADisec/sites.txt | awk '{print $1}')
nonASTADunique=$((nonASTADtotal-nonASTADidentical))

# Write cross-tabulation
echo -e "\t"non-ASTAD"\t"ASTAD > 3way-SNP-overlapSummary.txt
echo -e Identical"\t""${nonASTADidentical}""\t""${ASTADidentical}" >> 3way-SNP-overlapSummary.txt
echo -e Unique"\t""${nonASTADunique}""\t""${ASTADunique}" >> 3way-SNP-overlapSummary.txt
