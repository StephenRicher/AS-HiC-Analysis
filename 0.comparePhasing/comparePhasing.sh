c=11

estimate=../../GM12878/testChr"${c}"/alleleEstimate/snpsplit/GM12878-snpsplit.txt
truth=../../GM12878/testChr"${c}"/alleleTrue/snpsplit/GM12878-snpsplit.txt

if [ "${c}" = 22 ]; then
    paste <(cut -f 1-4 $estimate) <(rev $estimate | awk '{print $1}') > /tmp/estimate.txt
    estimate=/tmp/estimate.txt
fi

comm -12 <(sort "${truth}") <(sort "${estimate}") > phasing-Common-chr"${c}".txt

comm -23 <(sort "${truth}") <(sort "${estimate}") > phasing-TruthSetOnly-chr"${c}".txt

comm -13 <(sort "${truth}") <(sort "${estimate}") > phasing-EstimatedSetOnly-chr"${c}".txt

comm -12 \
    <(cut -f 2-3 "${truth}" | sort) \
    <(cut -f 2-3 "${estimate}" | sort) \
> phasing-MatchedPositions-chr"${c}".txt

comm -12 \
    <(cut -f 2-3 phasing-TruthSetOnly-chr"${c}".txt | sort) \
    <(cut -f 2-3 phasing-EstimatedSetOnly-chr"${c}".txt | sort) \
> phasing-MismatchedPositions-chr"${c}".txt
