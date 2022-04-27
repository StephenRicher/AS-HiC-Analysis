name=vennOverlap.svg
svg_stack.py --direction=h \
    <(python ../modifySVGLabel.py a TAD-overlap.svg) \
    <(python ../modifySVGLabel.py b TAD-overlap-percent.svg) \
> /tmp/"${name}"-ab.svg

svg_stack.py --direction=h \
    <(python ../modifySVGLabel.py c ASTAD-overlap.svg) \
    <(python ../modifySVGLabel.py d ASTAD-overlap-percent.svg) \
> /tmp/"${name}"-cd.svg

svg_stack.py --direction=v \
    /tmp/"${name}"-ab.svg /tmp/"${name}"-cd.svg \
> "${name}"


name=obsExpTADOverlap.svg
svg_stack.py --direction=h \
    <(python ../modifySVGLabel.py a GM_H1-TADOverlap.svg) \
    <(python ../modifySVGLabel.py b GM_IM-TADOverlap.svg) \
> /tmp/"${name}"-ab.svg

svg_stack.py --direction=h \
    <(python ../modifySVGLabel.py c IM_H1-TADOverlap.svg) \
    <(python ../modifySVGLabel.py d 3way-TADOverlap.svg) \
> /tmp/"${name}"-cd.svg

svg_stack.py --direction=v \
    /tmp/"${name}"-ab.svg /tmp/"${name}"-cd.svg \
> "${name}"

name=obsExpBaseOverlap.svg
svg_stack.py --direction=h \
    <(python ../modifySVGLabel.py a GM_H1-BaseOverlap.svg) \
    <(python ../modifySVGLabel.py b GM_IM-BaseOverlap.svg) \
> /tmp/"${name}"-ab.svg

svg_stack.py --direction=h \
    <(python ../modifySVGLabel.py c IM_H1-BaseOverlap.svg) \
    <(python ../modifySVGLabel.py d 3way-BaseOverlap.svg) \
> /tmp/"${name}"-cd.svg

svg_stack.py --direction=v \
    /tmp/"${name}"-ab.svg /tmp/"${name}"-cd.svg \
> "${name}"
