import pandas as pd
import seaborn as sns
from pybedtools import BedTool


def readDiffTAD(
        cell, diffOnly=False, X=True, normalOnly=False,
        excludeLowMap=False, excludeHighSignal=False, pyBed=False):
    assert cell in ['GM12878', 'IMR90', 'H1hESC', 'IMR-90', 'H1-hESC']
    if cell == 'IMR-90':
        cell = 'IMR90'
    elif cell == 'H1-hESC':
        cell = 'H1hESC'
    path = f'../0.getDiffDomain/{cell}-allTAD.pkl'
    names = ([
        'chrom', 'start', 'end', 'type', 'z', 'CNV',
        'High_Signal_Region', 'Low_Mappability'])
    df = pd.read_pickle(path)[names]
    if not X:
        df = df.loc[df['chrom'] != 'chrX']
    # Only retain differential TADs
    if diffOnly:
        df = df.loc[df['type'] == 'ASTAD']
    if normalOnly:
        df = df.loc[df['CNV'] == 'Normal']
    if excludeLowMap:
        df = df.loc[df['Low_Mappability'] == 0]
    if excludeHighSignal:
        df = df.loc[df['High_Signal_Region'] == 0]
    if pyBed:
        df = BedTool.from_dataframe(df)
    return df


def formatP(p):
    """ Return formatted p for title """
    if p > 0.999:
        pformat = '> 0.999'
    elif p < 0.001:
        pformat = '< .001'
    else:
        pformat = '= ' + f'{p:.3f}'[1:]
    if p < 0.001:
        return pformat + ' ***'
    elif p < 0.01:
        return pformat + ' **'
    elif p < 0.05:
        return pformat + ' *'
    else:
        return pformat



def formatCell(cell):
    """ Return formatted cell names """
    names = {'H1hESC': 'H1-hESC', 'IMR90': 'IMR-90', 'GM12878': 'GM12878'}
    if cell in names:
        return names[cell]
    else:
        return cell


def defaultPlotting(size=9, width=180, ratio=0.5):
    mm = 1 / 25.4 # mm in an inch
    colour = '#444444'
    sns.set(rc={
        'font.size': size, 'axes.titlesize': size, 'axes.labelsize': size,
        'xtick.labelsize': size, 'ytick.labelsize': size,
        'legend.fontsize': size, 'legend.title_fontsize': size,
        'font.family': 'sans-serif', 'lines.linewidth': 1.5,
        'axes.labelcolor': colour, 'xtick.color': colour,
        'ytick.color': colour,
        'figure.figsize': (width * mm, width * ratio * mm),
        'axes.spines.top': False, 'axes.spines.right': False,
    }, style='ticks')


def sampleByGroup(df, group, samplesPerGroup, replace=False):
    """ Perform random sampling with different numbers
        per group. i.e. maintain frequency by chromosome. """
    grouped = df.groupby(group)
    return (grouped.apply(lambda x: x.sample(samplesPerGroup[x.name], replace=replace)))


def readBlacklist():
    dtypes = {'chr': str, 'start': int, 'end': int, 'blacklist': str}
    return pd.read_csv(
        '../annotation/hg19-blacklist.bed',
        dtype=dtypes, names=dtypes.keys(), sep='\t')


def processBlacklist(domains):
    domains = domains.copy()
    domains['ID'] = range(0, len(domains))
    mergeBy = ['chrom', 'start', 'end', 'ID']
    domains_bed = BedTool.from_dataframe(domains[mergeBy])
    blacklist = readBlacklist()
    allOverlap = []
    noOverlap = []
    for reason, df in blacklist.groupby('blacklist'):
        names = ([
            'chrom', 'start', 'end', 'ID', 'chrom2',
            'start2', 'end2', reason, 'overlap'])
        bed = BedTool.from_dataframe(df)
        overlap = domains_bed.intersect(bed, wo=True).to_dataframe(names=names)
        if overlap.empty:
            noOverlap.append(reason)
            continue
        # Sum overlap per TAD
        overlap = (
            overlap.groupby(mergeBy)['overlap']
            .sum().reset_index().rename({'overlap': reason}, axis=1))
        overlap[reason] = overlap[reason] / (overlap['end'] - overlap['start'])
        allOverlap.append(overlap)

    if len(allOverlap) == 2:
        allOverlap = pd.merge(
            allOverlap[0], allOverlap[1],
            left_on=mergeBy, right_on=mergeBy, how='outer').fillna(0)
    elif len(allOverlap) == 1:
        allOverlap = allOverlap[0]
        allOverlap[noOverlap[0]] = 0
    else:
        domains[noOverlap[0]] = 0
        domains[noOverlap[1]] = 0
        return domains

    domains = (pd.merge(
        domains, allOverlap, left_on=mergeBy,
        right_on=mergeBy, how='left').drop('ID', axis=1))
    domains['High_Signal_Region'] = domains['High_Signal_Region'].fillna(0)
    domains['Low_Mappability'] = domains['Low_Mappability'].fillna(0)
    return domains
