#!/usr/bin/env Rscript

library(readr)
library(tidyverse)
library(ColorPalette)
library(grDevices)


args = commandArgs(trailingOnly=TRUE)
##args = c("~/vol/huge/exome_pools/download/M84902.pF3.pM3.exome_pools.tsv",
##         "~/vol/huge/exome_pools/download/M84902.pF3.pM3.exome_pools.bed")
##args = c("~/vol/pico/exome_pools_testing/download/M55774.pF2.pM2.exome_pools.tsv",
##         "~/vol/pico/exome_pools_testing/download/M55774.pF2.pM2.exome_pools.bed")
infile = args[[1]]
outfile = args[[2]]

d = read_tsv(infile)
head(d)
dim(d)


prob2rgb <- function(probability_list, hue=150, boost=1.5) {
  probability_list = sapply(probability_list, function(x) min(x * boost, 1))
  sapply(probability_list, function(x) gsub(" ", "", toString(col2rgb(hsv2rgb(hue,x,x))) )  )
}


df = d %>%
  mutate(Start = Start -1) %>%
  ##mutate(Name = sprintf("%s:%s-%s_%s_%s", Chr, Start, End, Ref, Alt)) %>%
  ##mutate(Name = sprintf("FATHERS:N_COV=%s,N_VAR=%s,N_VAF=%.3f,%s>%s",FATHERS_N_COV, FATHERS_N_VAR, FATHERS_VAF,Ref, Alt)) %>%
  mutate(Name = sprintf("FATHERS_ONLY__%s>%s__prob=%.3f",Ref, Alt, PROB_FATHERS_ONLY)) %>%
  mutate(Score = sprintf("%.3f", PROB_FATHERS_ONLY)) %>%
  mutate(Strand = ".") %>%
  mutate(ThickStart = Start) %>%
  mutate(ThickEnd = End) %>%
  mutate(ItemRgb = prob2rgb(PROB_FATHERS_ONLY, hue=200)) %>%
  select(Chr, Start, End, Name, Score, Score, Strand, ThickStart, ThickEnd, ItemRgb)
df

dm = d %>%
  mutate(Start = Start -1) %>%
  ##mutate(Name = sprintf("%s:%s-%s_%s_%s", Chr, Start, End, Ref, Alt)) %>%
  ##mutate(Name = sprintf("MOTHERS:N_COV=%s,N_VAR=%s,N_VAF=%.3f,%s>%s",MOTHERS_N_COV, MOTHERS_N_VAR, MOTHERS_VAF,Ref, Alt)) %>%
  mutate(Name = sprintf("MOTHERS_ONLY__%s>%s__prob=%.3f",Ref, Alt, PROB_MOTHERS_ONLY)) %>%
  mutate(Score = sprintf("%.3f", PROB_MOTHERS_ONLY)) %>%
  mutate(Strand = ".") %>%
  mutate(ThickStart = Start) %>%
  mutate(ThickEnd = End) %>%
  mutate(ItemRgb = prob2rgb(PROB_MOTHERS_ONLY, hue=300)) %>%
  select(Chr, Start, End, Name, Score, Score, Strand, ThickStart, ThickEnd, ItemRgb)
dm

dfm = d %>%
  mutate(Start = Start -1) %>%
  ##mutate(Name = sprintf("%s:%s-%s_%s_%s", Chr, Start, End, Ref, Alt)) %>%
  mutate(Name = sprintf("FATHERS_AND_MOTHERS__%s>%s__prob=%.3f",Ref, Alt, PROB_FATHERS_AND_MOTHERS)) %>%
  mutate(Score = sprintf("%.3f", PROB_FATHERS_AND_MOTHERS)) %>%
  mutate(Strand = ".") %>%
  mutate(ThickStart = Start) %>%
  mutate(ThickEnd = End) %>%
  mutate(ItemRgb = prob2rgb(PROB_FATHERS_AND_MOTHERS, hue=80)) %>%
  select(Chr, Start, End, Name, Score, Score, Strand, ThickStart, ThickEnd, ItemRgb)
dfm


write_tsv(rbind(dfm, dm, df), outfile, col_names=F)
