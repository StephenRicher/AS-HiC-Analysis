for cell in GM12878 IMR90 H1hESC; do
    python makeCNV-UCSC.py ../../"${cell}"/"${cell}"-CNV.txt \
    > "${cell}"-CNV-UCSC.bed
    grep -v Normal "${cell}"-CNV-UCSC.bed | bedtools sort | bedtools merge > "${cell}"-CNV-nonNormal.bed
done
