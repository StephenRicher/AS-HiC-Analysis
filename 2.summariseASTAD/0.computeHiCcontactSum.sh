> allHiC-contactSum.txt
hicInfo=/media/stephen/Elements/envs/b24ae601f4d79c98e2559781f0fef128/bin/hicInfo
for chrom in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22; do
    echo "${chrom}"
    "${hicInfo}" --matrices ../../GM12878/fullGRCh37/dat/matrix/chr"${chrom}"/base/raw/GM12878-chr"${chrom}".10000-full.h5 >> allHiC-contactSum.txt
    "${hicInfo}" --matrices ../../IMR90/fullGRCh37/dat/matrix/chr"${chrom}"/base/raw/IMR90-chr"${chrom}".10000-full.h5 >> allHiC-contactSum.txt
    "${hicInfo}" --matrices ../../H1hESC/fullGRCh37/dat/matrix/chr"${chrom}"/base/raw/H1hESC-chr"${chrom}".10000-full.h5 >> allHiC-contactSum.txt
    "${hicInfo}" --matrices ../../GM12878/alleleGRCh37/dat/matrix/chr"${chrom}"/base/raw/GM12878_a[12]-chr"${chrom}".10000-SNPsplit.h5 >> allHiC-contactSum.txt
    "${hicInfo}" --matrices ../../IMR90/alleleGRCh37/dat/matrix/chr"${chrom}"/base/raw/IMR90_a[12]-chr"${chrom}".10000-SNPsplit.h5 >> allHiC-contactSum.txt
    "${hicInfo}" --matrices ../../H1hESC/alleleGRCh37/dat/matrix/chr"${chrom}"/base/raw/H1hESC_a[12]-chr"${chrom}".10000-SNPsplit.h5 >> allHiC-contactSum.txt
done
