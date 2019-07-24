#!/usr/bin/env Rscript 
# Parse Omar/Danish's sift_mat_all.Rdata into a clean table

library(tidyverse)
library(siftr)

load('data/mutfunc/human/conservation/sift_mat_all.Rdata')

tbl <- sapply(sift_mat_all, siftr::filterPredictions, simplify = FALSE, score_thresh=Inf, ic_thresh=3.25, residue_thresh=2) %>%
  bind_rows(.id = 'acc') %>%
  as_tibble()

write_tsv(tbl, 'data/mutfunc/human/conservation/sift_all.tab')
