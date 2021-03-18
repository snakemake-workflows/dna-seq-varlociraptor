#!/usr/bin/env Rscript

library(readr)
library(tidyverse)


args = commandArgs(trailingOnly=TRUE)
##args = c("~/vol/huge/exome_pools/dna-seq-varlociraptor/results/candidates/hg19/M55774.hg19_candidates.tsv",
##         "~/vol/huge/exome_pools/dna-seq-varlociraptor/results/candidates/hg19/M55774.hg19_candidates.core.tsv")
infile = args[[1]]
outfile = args[[2]]

d = read_tsv(infile)
head(d)
dim(d)

d2 = d %>%
  select(Chr, Start, End, Ref, Alt)

write_tsv(d2, outfile)

