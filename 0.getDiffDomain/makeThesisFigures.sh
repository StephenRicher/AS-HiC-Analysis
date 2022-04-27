name=subtractionExample.svg
svg_stack.py --direction=h Subtraction-a.svg Subtraction-b.svg \
> /tmp/"${name}"-ab.svg

svg_stack.py --direction=h Subtraction-c.svg Subtraction-d.svg \
> /tmp/"${name}"-cd.svg

svg_stack.py --direction=v \
    /tmp/"${name}"-ab.svg /tmp/"${name}"-cd.svg \
> "${name}"

rsvg-convert -f pdf -o "${name/.svg/.pdf}" "${name}"

name=ExampleASTAD.svg
svg_stack.py --direction=h ExampleASTAD-a.svg ExampleASTAD-b.svg \
> /tmp/"${name}"-ab.svg

svg_stack.py --direction=h ExampleASTAD-c.svg ExampleASTAD-d.svg \
> /tmp/"${name}"-cd.svg

svg_stack.py --direction=v \
    /tmp/"${name}"-ab.svg /tmp/"${name}"-cd.svg \
> "${name}"

rsvg-convert -f pdf -o "${name/.svg/.pdf}" "${name}"

