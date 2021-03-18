#!/usr/bin/env Rscript

library(readr)
library(tidyverse)
library(ColorPalette)
library(grDevices)


args = commandArgs(trailingOnly=TRUE)
##args = c("~/vol/huge/exome_pools/download/M84902.pF3.pM3.exome_pools.tsv",
##         "~/vol/huge/exome_pools/download/M84902.pF3.pM3.exome_pools.bed")
infile = args[[1]]
outfile = args[[2]]

d = read_tsv(infile)
head(d)
dim(d)


prob2rgb <- function(hue, probability_list) {
  sapply(probability_list, function(x) gsub(" ", "", toString(col2rgb(hsv2rgb(hue,x,x))) )  )
}
##prob2rgb(150, 0.6)

df = d %>%
  mutate(Start = Start -1) %>%
  ##mutate(Name = sprintf("%s:%s-%s_%s_%s", Chr, Start, End, Ref, Alt)) %>%
  ##mutate(Name = sprintf("FATHERS:N_COV=%s,N_VAR=%s,N_VAF=%.3f,%s>%s",FATHERS_N_COV, FATHERS_N_VAR, FATHERS_VAF,Ref, Alt)) %>%
  mutate(Name = sprintf("FATHERS_ONLY__%s>%s__prob=%.3f",Ref, Alt, PROB_FATHERS_ONLY)) %>%
  mutate(Score = sprintf("%.3f", PROB_FATHERS_ONLY)) %>%
  mutate(Strand = ".") %>%
  mutate(ThickStart = Start) %>%
  mutate(ThickEnd = End) %>%
  mutate(ItemRgb = prob2rgb(200, PROB_FATHERS_ONLY)) %>%
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
  mutate(ItemRgb = prob2rgb(300, PROB_MOTHERS_ONLY)) %>%
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
  mutate(ItemRgb = prob2rgb(80, PROB_FATHERS_AND_MOTHERS)) %>%
  select(Chr, Start, End, Name, Score, Score, Strand, ThickStart, ThickEnd, ItemRgb)
dfm


write_tsv(rbind(dfm, dm, df), outfile, col_names=F)
