awk '$1=="chr12"' sharedVCF-withGWAS.tsv | awk '$2>=52430000' | awk '$2<=53490000' > KRT1-conservedTAD.tsv
awk '$1=="chr11"' sharedVCF-withGWAS.tsv | awk '$2>=1500000' | awk '$2<=2300000' > H19-conservedTAD.tsv
