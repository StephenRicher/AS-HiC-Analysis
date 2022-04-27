#!/usr/bin/env Rscript

library('LOLA')
library('data.table')
library('GenomicRanges')

# Set path to path of script (in RStudio)
args = commandArgs(trailingOnly=TRUE)
setwd(args)

for (cell in c('IMR90', 'GM12878', 'H1hESC')) {
  # Load region database
  regionDB = loadRegionDB(paste0('data/LOLAregionDB/hg19/', cell, '/'))


  # Load Universe and disjoin to remove overlaps
  activeDHS = readBed(paste0('../0.getDiffDomain/', cell, '-TADs-autosome.bed'))
  activeDHS = disjoin(activeDHS)

  # Load diff domains and redefine for more appropriate enrichment
  regionSet = readBed(paste0('../0.getDiffDomain/', cell, '-ASTADs-autosome.bed'))
  userSets = GRangesList(regionSet)
  userSetsRedefined = redefineUserSets(userSets, activeDHS)

  # Ensure universe is appropriate
  checkUniverseAppropriateness(userSetsRedefined, activeDHS)

  # Run analysis
  res = runLOLA(userSetsRedefined, activeDHS, regionDB, cores=1, direction='enrichment')

  # Recompute p-value from -log10(p)
  res[, 'p':=10^(-pValueLog)]

  # Write to csv
  fwrite(res, file=paste0(cell, '-fullVCF-LOLA-results.csv'))
}
