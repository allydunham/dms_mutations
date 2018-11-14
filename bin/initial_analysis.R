#!/usr/bin/env Rscript 
# Preliminary analysis of deep mutagenesis data

library(tidyverse)
library(magrittr)
library(ggpubr)

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

pab1_sift <- read_tsv('data/yeast_pab1_sift.tsv') %>%
  mutate(mut_id = paste(acc, pos, alt, sep = '_'))
pab1_foldx <- read_tsv('data/yeast_pab1_foldx.tsv') %>%
  mutate(mut_id = paste(uniprot_id, uniprot_pos, aa_mt, sep = '_'))

braf_sift <- read_tsv('data/human_braf_sift.tsv') %>%
  mutate(mut_id = paste(acc, pos, alt, sep = '_'))
braf_foldx <- read_tsv('data/human_braf_foldx.tsv') %>%
  mutate(mut_id = paste(uniprot_id, uniprot_pos, aa_mt, sep = '_'))


## Add Sift/FoldX Scores to Data
hietpas_2011 <- left_join(deep_mut_data$hietpas_2011_hsp90, select(hsp90_sift, ref, score, median_ic, n_aa, n_seq, mut_id), by = 'mut_id') %>%
  left_join(., select(hsp90_foldx, aa_wt, ddG, ddG_sd, mut_id), by='mut_id')

roscoe_2013 <- left_join(deep_mut_data$roscoe_2013_ubi, select(ubi_sift, ref, score, median_ic, n_aa, n_seq, mut_id), by = 'mut_id') %>%
  left_join(., select(ubi_foldx, aa_wt, ddG, ddG_sd, mut_id), by='mut_id')

jiang_2013 <- left_join(deep_mut_data$jiang_2013_hsp90, select(hsp90_sift, ref, score, median_ic, n_aa, n_seq, mut_id), by = 'mut_id') %>%
  left_join(., select(hsp90_foldx, aa_wt, ddG, ddG_sd, mut_id), by='mut_id')

melamed_2013 <- left_join(deep_mut_data$melamed_2013_pab1, select(pab1_sift, ref, score, median_ic, n_aa, n_seq, mut_id), by = 'mut_id') %>%
  left_join(., select(pab1_foldx, aa_wt, ddG, ddG_sd, mut_id), by='mut_id')

wagenaar_2014 <- left_join(deep_mut_data$wagenaar_2014_braf, select(braf_sift, ref, score, median_ic, n_aa, n_seq, mut_id), by = 'mut_id') %>%
  left_join(., select(braf_foldx, aa_wt, ddG, ddG_sd, mut_id), by='mut_id')

## Sift/FoldX Plots
plot_metrics <- function(tbl, colname, ylabel=NULL){
  if (!all(is.na(tbl$score))){
    p_sift <- ggplot(tbl, aes_string(x='score', y=colname)) +
      geom_point() +
      geom_density2d() + 
      xlab('SIFT Score') + 
      ylab(ifelse(is.null(ylabel), colname, ylabel))
  } else {
    p_sift <- NA
  }
  
  if (!all(is.na(tbl$ddG))){
    p_foldx <- ggplot(tbl, aes_string(x='ddG', y=colname)) +
      geom_point() +
      geom_density2d() + 
      xlab('ddG') + 
      ylab('')
  } else {
    p_foldx <- NA
  }
  
  plots <- list(p_sift, p_foldx)
  plots <- plots[!is.na(plots)]
  if (length(plots) > 0){
    return(ggarrange(plotlist = plots))
  } else {
    warning('Neither SIF scores nor FoldX ddG scores avilable')
    return(NULL)
  }
}

p_hietpas_2011 <- plot_metrics(hietpas_2011, 'selection_coefficient')
p_roscoe_2013 <- plot_metrics(roscoe_2013, 'selection_num')
p_jiang_2013 <- plot_metrics(jiang_2013, 'average_num')
p_melamed_2013 <- plot_metrics(melamed_2013, 'enrichment_ratio')
p_wagenaar_2014 <- plot_metrics(wagenaar_2014, 'median_enrichment')
