#!/usr/bin/env Rscript
# Run siftr predictor on input alignments. Expects alignments to be named as uniprotID.alignedfasta

library(siftr)
library(tidyverse)

alignment_files <- commandArgs(trailingOnly = TRUE)

score_mats <- sapply(alignment_files, siftr::predictFromAlignmentFreqNormed)
names(score_mats) <- str_split(basename(alignment_files), '\\.', simplify = TRUE)[,1]

filtered_table <- sapply(score_mats, siftr::filterPredictions, score_thresh=Inf, simplify = FALSE) %>%
  bind_rows(.id = 'uniprot_id')

writeLines(format_tsv(filtered_table))
