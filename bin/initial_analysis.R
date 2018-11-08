#!/usr/bin/env Rscript 
# Preliminary analysis of deep mutagenesis data

library(tidyverse)
library(magrittr)

# Not implemented yet as no data need
#source('bin/analysis_functions.R')
#source('bin/load_processed_data.R')

# Currently process_data simply imports and minimally processes test studies
source('bin/process_data.R')

## Import temp sift/foldx values for initial examinations
hsp90_sift <- read_tsv('data/yeast_hsp90_sift.tsv') %>%
  mutate(mut_id = paste(acc, pos, alt, sep = '_'))
hsp90_foldx <- read_tsv('data/yeast_hsp90_foldx.tsv') %>%
  mutate(mut_id = paste(uniprot_id, uniprot_pos, aa_mt, sep = '_'))

ubi_sift <- read_tsv('data/yeast_ubi_sift.tsv') %>%
  mutate(mut_id = paste(acc, pos, alt, sep = '_'))
ubi_foldx <- read_tsv('data/yeast_ubi_foldx.tsv') %>%
  mutate(mut_id = paste(uniprot_id, uniprot_pos, aa_mt, sep = '_'))

## Add Sift/FoldX Scores to Data
deep_mut_data$hietpas_2011_hsp90 %<>% left_join(., select(hsp90_sift, ref, score, median_ic, n_aa, n_seq, mut_id), by = 'mut_id') %>%
  left_join(., select(hsp90_foldx, aa_wt, ddG, ddG_sd, mut_id), by='mut_id')

deep_mut_data$roscoe_2013_ubi %<>% left_join(., select(ubi_sift, ref, score, median_ic, n_aa, n_seq, mut_id), by = 'mut_id') %>%
  left_join(., select(ubi_foldx, aa_wt, ddG, ddG_sd, mut_id), by='mut_id')

## Sift/FoldX Plots
p_hietpas_2011_sift <- ggplot(deep_mut_data$hietpas_2011_hsp90, aes(x=score, y=selection_coefficient)) +
  geom_point() +
  geom_density2d() + 
  xlab('SIFT Score') + 
  ylab('Selection Coef')

p_hietpas_2011_foldx <- ggplot(deep_mut_data$hietpas_2011_hsp90, aes(x=ddG, y=selection_coefficient)) +
  geom_point() +
  geom_density2d() + 
  xlab('FoldX ddG') + 
  ylab('Selection Coef')

p_roscoe_2013_sift <- ggplot(deep_mut_data$roscoe_2013_ubi, aes(x=score, y=selection_num)) +
  geom_point() +
  geom_density2d() + 
  xlab('SIFT Score') + 
  ylab('Selection Coef')

p_roscoe_2013_foldx <- ggplot(deep_mut_data$roscoe_2013_ubi, aes(x=ddG, y=selection_num)) +
  geom_point() +
  geom_density2d() + 
  xlab('FoldX ddG') + 
  ylab('Selection Coef')
