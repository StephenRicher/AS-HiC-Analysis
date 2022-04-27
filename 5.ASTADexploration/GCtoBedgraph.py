#!/usr/bin/env python3

import sys

validChroms = [str(x) for x in range(1, 23)]
span = 5
for line in sys.stdin:
    line = line.strip()
    if not line:
        continue
    elif line.startswith('variableStep'):
        chrom = line.split()[1].split('=')[1].strip('chr')
        span = int(line.split()[2].split('=')[1])
    elif chrom in validChroms:
        start, gc = line.split()
        print(chrom, start, int(start) + span, gc, sep='\t')
