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
import_sift <- function(path){
  return(
    read_tsv(path) %>% mutate(mut_id = paste(acc, pos, alt, sep = '_'))
  )
}

import_foldx <- function(path){
  return(
    read_tsv(path) %>% mutate(mut_id = paste(uniprot_id, uniprot_pos, aa_mt, sep = '_'))
  )
}

sift <- sapply(c('data/yeast_hsp90_sift.tsv',
                 'data/yeast_ubi_sift.tsv',
                 'data/yeast_pab1_sift.tsv',
                 'data/human_braf_sift.tsv'),
               import_sift, simplify = FALSE) %>%
  set_names(gsub('(data/|\\_sift\\.tsv)', '', names(.)))

foldx <- sapply(c('data/yeast_hsp90_foldx.tsv',
                  'data/yeast_ubi_foldx.tsv',
                  'data/yeast_pab1_foldx.tsv',
                  'data/human_braf_foldx.tsv'),
                import_foldx, simplify = FALSE) %>%
  set_names(gsub('(data/|\\_foldx\\.tsv)', '', names(.)))


## Add Sift/FoldX Scores to Data
join_metrics <- function(tbl, gene){
  return(
    left_join(deep_mut_data[[tbl]], select(sift[[gene]], ref, score, median_ic, n_aa, n_seq, mut_id), by = 'mut_id') %>%
      left_join(., select(foldx[[gene]], aa_wt, ddG, ddG_sd, mut_id), by='mut_id')
  )
}

studies <- mapply(join_metrics,
                  c('hietpas_2011_hsp90', 'roscoe_2013_ubi', 'jiang_2013_hsp90', 'melamed_2013_pab1', 'wagenaar_2014_braf'),
                  c('yeast_hsp90', 'yeast_ubi', 'yeast_hsp90', 'yeast_pab1', 'human_braf'),
                  SIMPLIFY = FALSE)

## Sift/FoldX Plots
plot_metrics <- function(tbl, colname, ylabel=NULL, title=''){
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
    return(ggarrange(plotlist = plots) %>% annotate_figure(., top=text_grob(title)))
  } else {
    warning('Neither SIF scores nor FoldX ddG scores avilable')
    return(NULL)
  }
}

p_metrics <- mapply(function(x, y){plot_metrics(studies[[x]], y, title=x)},
                    names(studies),
                    c('selection_coefficient', 'selection_num', 'average_num', 'enrichment_ratio', 'median_enrichment'),
                    SIMPLIFY = FALSE)

for (n in names(p_metrics)){
  ggsave(paste0('figures/initial_analysis/', n, '_sift_foldx.pdf'), p_metrics[[n]], width = 7, height = 5)
}
