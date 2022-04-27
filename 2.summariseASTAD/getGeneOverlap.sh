bedtools intersect -wa \
  -a results/GM12878/alleleGRCh37/GM12878_a1-vs-GM12878_a2-all-20000-overlapASTAD.gff3 \
  -b results/IMR90/alleleGRCh37/IMR90_a1-vs-IMR90_a2-all-20000-overlapASTAD.gff3 \
| bedtools sort \
| uniq \
> gencode.v38lift37-genes-ASTADoverlap_IMR90-GM12878.gff3
