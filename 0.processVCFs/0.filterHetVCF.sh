# Move copies of VCF to here and standardise naming syntax
cp ../../IMR90/genome/IMR90-all.filt.vcf.gz IMR90-all.vcf.gz
cp ../../H1hESC/genome/H1hESC-all.filt.vcf.gz H1hESC-all.vcf.gz
bgzip -c ../../GM12878/genome/NA12878-hg19-chrPrefix.vcf > GM12878-all.vcf.gz

# Extract heterozygous variants only
for cell in GM12878 IMR90 H1hESC; do
    tabix "${cell}"-all.vcf.gz
    bcftools view -g het "${cell}"-all.vcf.gz \
    | bgzip > "${cell}"-hetOnly.vcf.gz
    tabix "${cell}"-hetOnly.vcf.gz
done

# Write UCSC compatible format for visualisation
for cell in GM12878 IMR90 H1hESC; do
    echo track name="Variants (${cell}) description="Green: SNP, Purple: Indel" db=hg19 visibility=1 itemRgb="On"" \
    > /tmp/"${cell}"-variants-UCSC.bed
    bcftools view -H -v indels -g het "${cell}"-all.vcf.gz \
    | awk -v OFS="\t" '{print $1, $2, $2+1, "Indel-Het", "0", ".",  $2, $2+1, "190,174,212"}' \
    >> /tmp/"${cell}"-variants-UCSC.bed
    bcftools view -H -v snps -g het "${cell}"-all.vcf.gz \
    | awk -v OFS="\t" '{print $1, $2, $2+1, "SNP-Het", "0", ".",  $2, $2+1, "127,201,127"}' \
    >> /tmp/"${cell}"-variants-UCSC.bed
    bcftools view -H -v indels -g hom "${cell}"-all.vcf.gz \
    | awk -v OFS="\t" '{print $1, $2, $2+1, "Indel-Hom", "0", ".",  $2, $2+1, "253,192,134"}' \
    >> /tmp/"${cell}"-variants-UCSC.bed
    bcftools view -H -v snps -g hom "${cell}"-all.vcf.gz \
    | awk -v OFS="\t" '{print $1, $2, $2+1, "SNP-Hom", "0", ".",  $2, $2+1, "255,255,153"}' \
    >> /tmp/"${cell}"-variants-UCSC.bed
    bedtools sort -i /tmp/"${cell}"-variants-UCSC.bed \
    > "${cell}"-variants-UCSC.bed
done

for cell in GM12878 IMR90 H1hESC; do
    for mode in phased all; do
        if [ "${mode}" = phased ]; then
            f="--phased"
        else
            f=""
        fi
        bcftools view -H -m 2 -M 2 -v snps -g het $f "${cell}"-all.vcf.gz \
        | awk '{print $1}' | uniq -c \
        | awk -v OFS="\t" -v cell="${cell}" -v mode="${mode}" '{print cell, mode, $2, $1}'
    done
done > allPhasingData.tsv

# Generate genome wide haplotype sequence
for cell in GM12878 IMR90 H1hESC; do
    for a in 1 2; do
        zcat ../../annotation/GATKresourceBundle/ucsc.hg19.fasta.gz \
        | bcftools consensus --haplotype "${a}"pIu "${cell}"-all.vcf.gz \
        | bgzip > "${cell}"-a"${a}".fasta.gz
        samtools faidx "${cell}"-a"${a}".fasta.gz
    done
done

# Intersect heterozygous variants to find which sites
# are similar across all cell lines
bcftools isec -p intersection/ --nfiles +1 --collapse none \
    GM12878-hetOnly.vcf.gz H1hESC-hetOnly.vcf.gz IMR90-hetOnly.vcf.gz

# Reformat intersection output
awk -v OFS="\t" \
    'BEGIN {print "#chrom", "pos", "REF", "ALT", "GM12878", "H1hESC", "IMR90"} \
    {print $1, $2, $3, $4, substr($5,1,1), substr($5,2,1), substr($5,3,1)}' \
sites.txt > sharedVCF.tsv

# Intersect with GWAS variants
bedtools intersect -wa \
    -a ../annotation/gwas_catalog_v1.0-associations_e104_r2021-11-22.hg19.bed.gz \
    -b <(awk -v OFS="\t" '{print $1, $2, $2+1}' sharedVCF.tsv) \
> gwasOverlapVar.bed
