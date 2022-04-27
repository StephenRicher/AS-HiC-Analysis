#!/usr/bin/env python3

""" Convert get_CNV output to BED format with colour for UCSC """

import os
import sys
import pandas as pd

def setColour(x):
    if x == 2:
        return '255,0,0'
    elif x == 1:
        return '255,152,152'
    elif x == 0:
        return '255,255,255'
    elif x == -1:
        return '152,152,255'
    elif x == -2:
        return '0,0,255'
    else:
        # Unknown - should not happen
        return '0,0,0'

def setName(x):
    if x == 2:
        return 'Amplification'
    elif x == 1:
        return 'Gain'
    elif x == 0:
        return 'Normal'
    elif x == -1:
        return 'Loss'
    elif x == -2:
        return 'Deletion'
    else:
        # Unknown - should not happen
        return 'Unknown'

file = sys.argv[1]
cell = os.path.basename(file).split('-')[0]
names = {'chrom': str, 'start': int, 'end': int, 'CNV': int}
cnv = pd.read_csv(file, usecols=[1,2,3,4], skiprows=1, names=names.keys(), dtype=names, sep='\t')
cnv['start'] = cnv['start'] - 1 # Convert to 0-based
cnv['name'] = cnv['CNV'].apply(setName)
cnv['blank'] = '.'
cnv['score'] = 0
cnv['colour'] = cnv['CNV'].apply(setColour)
names = ['chrom', 'start', 'end', 'name', 'score', 'blank', 'start', 'end', 'colour']

print(f'track name="CNV ({cell})" description="Dark Blue: deletion, Light blue: Loss, White: Normal, Light Red: Gain, Dark Red: Amplification" db=hg19 visibility=1 itemRgb="On"')
print('browser position chr11')
cnv[names].to_csv(sys.stdout, index=False, header=False, sep='\t')
