name=ASTADsummary.svg
svg_stack.py --direction=h \
    <(python ../modifySVGLabel.py a ASTAD-Frequency-By-Chrom.svg) \
    <(python ../modifySVGLabel.py b PhasedSNP_Density.svg) \
> /tmp/"${name}"-ab.svg

svg_stack.py --direction=h \
    <(python ../modifySVGLabel.py c FullContactDensity.svg) \
    <(python ../modifySVGLabel.py d AlleleContactDensity.svg) \
> /tmp/"${name}"-cd.svg

svg_stack.py --direction=v \
    /tmp/"${name}"-ab.svg /tmp/"${name}"-cd.svg \
> "${name}"
