# Split chromatin state marks by in sorted BED files
# Source: https://egg2.wustl.edu/roadmap/web_portal/imputed.html#chr_imp
# Core 15-state model (5 marks, 127 epigenomes)

for ID in E003 E017 E116; do
    if [ "$ID" == "E003" ]; then
        cell=H1hESC
    elif [ "$ID" == "E017" ]; then
        cell=IMR90
    else
        cell=GM12878
    fi
    zcat "${ID}"_15_coreMarks_mnemonics.bed.gz \
      | awk '$4=="1_TssA"' \
      | bedtools sort > "${cell}"-chromHMM-ActiveTSS.bed
    zcat "${ID}"_15_coreMarks_mnemonics.bed.gz \
      | awk '$4=="2_TssAFlnk"' \
      | bedtools sort > "${cell}"-chromHMM-FlankingActiveTSS.bed
    zcat "${ID}"_15_coreMarks_mnemonics.bed.gz \
      | awk '$4=="3_TxFlnk"' \
      | bedtools sort > "${cell}"-chromHMM-TranscriptionAt5_3.bed
    zcat "${ID}"_15_coreMarks_mnemonics.bed.gz \
      | awk '$4=="4_Tx"' \
      | bedtools sort > "${cell}"-chromHMM-StrongTranscription.bed
    zcat "${ID}"_15_coreMarks_mnemonics.bed.gz \
      | awk '$4=="5_TxWk"' \
      | bedtools sort > "${cell}"-chromHMM-WeakTranscription.bed
    zcat "${ID}"_15_coreMarks_mnemonics.bed.gz \
      | awk '$4=="6_EnhG"' \
      | bedtools sort > "${cell}"-chromHMM-GenicEnhancer.bed
    zcat "${ID}"_15_coreMarks_mnemonics.bed.gz \
      | awk '$4=="7_Enh"' \
      | bedtools sort > "${cell}"-chromHMM-Enhancers.bed
    zcat "${ID}"_15_coreMarks_mnemonics.bed.gz \
      | awk '$4=="8_ZNF/Rpts"' \
      | bedtools sort > "${cell}"-chromHMM-ZNFgenesRepeats.bed
    zcat "${ID}"_15_coreMarks_mnemonics.bed.gz \
      | awk '$4=="9_Het"' \
      | bedtools sort > "${cell}"-chromHMM-Heterochromatin.bed
    zcat "${ID}"_15_coreMarks_mnemonics.bed.gz \
      | awk '$4=="10_TssBiv"' \
      | bedtools sort > "${cell}"-chromHMM-BivalentPoisedTSS.bed
    zcat "${ID}"_15_coreMarks_mnemonics.bed.gz \
      | awk '$4=="11_BivFlnk"' \
      | bedtools sort > "${cell}"-chromHMM-FlankingBivalentTSS.bed
    zcat "${ID}"_15_coreMarks_mnemonics.bed.gz \
      | awk '$4=="12_EnhBiv"' \
      | bedtools sort > "${cell}"-chromHMM-BivalentEnhancer.bed
    zcat "${ID}"_15_coreMarks_mnemonics.bed.gz \
      | awk '$4=="13_ReprPC"' \
      | bedtools sort > "${cell}"-chromHMM-RepressedPolycomb.bed
    zcat "${ID}"_15_coreMarks_mnemonics.bed.gz \
      | awk '$4=="14_ReprPCWk"' \
      | bedtools sort > "${cell}"-chromHMM-WeakRepressedPolycomb.bed
    zcat "${ID}"_15_coreMarks_mnemonics.bed.gz \
      | awk '$4=="15_Quies"' \
      | bedtools sort > "${cell}"-chromHMM-Quiescent.bed
done
