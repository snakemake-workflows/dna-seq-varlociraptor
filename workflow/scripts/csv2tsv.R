#!/usr/bin/env Rscript

library(readr)

args = commandArgs(trailingOnly=TRUE)
infile = args[[1]]
outfile = args[[2]]

d = read_csv(infile)
d

write_tsv(d, outfile)

