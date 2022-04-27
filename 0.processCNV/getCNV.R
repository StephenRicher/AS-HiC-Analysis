#!/usr/bin/env Rscript

library(HiCcompare)

args = commandArgs(trailingOnly=TRUE)
bamPath = args[1]
outCNV = args[2]

cnv <- get_CNV(path2bam = bamPath, out.file = outCNV,
               bin.size = 20, genome = 'hg19', CNV.level = 2)
