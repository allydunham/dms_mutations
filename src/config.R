#!/usr/bin/env Rscript 
# Load libraries and configurations for all mutations analyses, to keep consistency within project

# Load all libraries (keeps namespace consistent)
library(Biostrings)
library(tidyverse)
library(rlang)
library(magrittr)
library(broom)
library(readxl)
library(ggpubr)
library(GGally)
library(scales)
library(pdist)
library(caret)
library(Rtsne)

# Source custom functions
source('src/util/misc_util.R')
source('src/util/autoplot.R')
source('src/data_processing/deep_mut_data.R')

#### Project config ####
# Constants
MUT_SCORE_NAME <- 'Standardised Mutagenesis Score'
RAW_MUT_SCORE_NAME <- 'Raw Mutagenesis Score'

MUT_CATEGORIES <- list(deleterious='deleterious', neutral='neutral', beneficial='beneficial')

# Categories of secondary structure
SS_REDUCED_HASH <- c(C='None', S='Turn', H='Helix', T='Turn', E='Strand', G='Helix', B='Strand', I='Helix')

# Categories of amino acid
# From Sigma-Aldrich website (https://www.sigmaaldrich.com/life-science/metabolomics/learning-center/amino-acid-reference-chart.html)
AA_REDUCED_CLASSES <- list(Aliphatic=c('A', 'I', 'L', 'M', 'V'),
                           Aromatic=c('F', 'W', 'Y'),
                           Polar=c('N', 'C', 'Q', 'S', 'T'),
                           Basic=c('R', 'H', 'K'),
                           Acidic=c('D', 'E'),
                           Glycine='G',
                           Proline='P')
AA_REDUCED_HASH <- structure(rep(names(AA_REDUCED_CLASSES), times=sapply(AA_REDUCED_CLASSES, length)),
                             names=unlist(AA_REDUCED_CLASSES))

## AA colours, based roughly on groups
AA_COLOURS <- c(A='red', I='salmon', L='firebrick', M='orange', V='tomato',
             F='gold', W='yellow3', Y='khaki3',
             N='cadetblue1', C='cornflowerblue', Q='cyan', S='blue', T='darkslateblue',
             R='green', H='green4', K='seagreen1',
             D='purple', E='pink',
             G='antiquewhite2', P='black', X='grey')


# Manual deleterious thresholds per dms study
MANUAL_THRESHOLDS <- c(araya_2012_yap1=-2, ashenberg_2016_nucleoprotein=-1, brenan_2016_erk2=-0.5, findlay_2018_brca1=-1,
                       firnberg_2014_tem1=-1.5, giacomelli_2018_tp53.p53_null_etoposide=-1, giacomelli_2018_tp53.p53_null_nutlin3=-0.75,
                       giacomelli_2018_tp53.p53_wt_nutlin3=-1, hietpas_2011_hsp90=-0.5, jiang_2013_hsp90=-0.5, kitzman_2015_gal4=-2.5,
                       lee_2018_ha=-3.75, melamed_2013_pab1=-1.25, melnikov_2014_aph3_ii.ami=-1, melnikov_2014_aph3_ii.g418=-0.75,
                       melnikov_2014_aph3_ii.kan=-2, melnikov_2014_aph3_ii.neo=-2, melnikov_2014_aph3_ii.paro=-2,
                       melnikov_2014_aph3_ii.ribo=-1.25, mishra_2016_hsp90=-0.5, olson_2014_proteing=-5, roscoe_2013_ubi4=-0.25,
                       roscoe_2014_ubi4.excess_e1=-1, roscoe_2014_ubi4.limiting_e1=-1.5, starita_2013_ube4b=-2.5,
                       starita_2015_brca1.bard1_ring_binding=-0.25, starita_2015_brca1.e3_activity=-1, wagenaar_2014_braf=0,
                       weile_2017_sumo1=-3, weile_2017_ube2i=-2)

# ggplot theme (currently unset)
theme_set(theme_gray()) 

# Function to format study names
format_study <- function(x, max_width=60){
  spl <- str_split(x, '[_.]')[[1]]
  
  # Year
  year <- spl[2]
  
  # Authour
  if (all(spl[1:2] == c('roscoe', '2014'))){
    authour <- 'Roscoe & Bolon'
  } else {
    authour <- str_c(str_to_title(spl[1]), ' et al.')
  }
  
  # Gene
  gene <- spl[3]
  if (gene %in% c('nucleoprotein', 'proteing')){
    gene <- c(nucleoprotein='NP', proteing='ProteinG')[gene]
  } else if (gene == 'aph3'){
    gene <- 'APH3-II'
    spl <- spl[c(1, 2, 3, 5:length(spl))] # remove erroneously split ii
  } else {
    gene <- str_to_upper(gene)
  }
  
  # Condition
  condition <- NA
  if (!is.na(spl[4])){
    condition <- str_c(spl[4:length(spl)], collapse = ' ')
  }
  
  result <- str_c(authour, ' ', year, ' (', gene, ifelse(is.na(condition), '', str_c(', ', condition)), ')') %>%
    str_wrap(width = max_width)
  
  return(result)
}

pretty_p_values <- function(p, breaks = c(0.001, 0.01, 0.05), markdown_exp=FALSE, prefix_p=FALSE){
  breaks <- sort(breaks, decreasing = TRUE)
  break_str <- as.character(breaks)
  gt <- '>'
  lt <- '<'
  
  if (markdown_exp){
    break_str <- str_replace(break_str, "^1?e(-?[0-9]*\\.?[0-9]*)", "10<sup>\\1</sup>")
    gt <- '&gt;'
    lt <- '&lt;'
  }
  
  p_out <- rep(str_c(gt, ' ', break_str[1]), length(p))
  for (i in 1:length(breaks)){
    p_out[p < breaks[i]] <- str_c(lt, ' ', break_str[i])
  }
  break_levels <- c(str_c(gt, ' ', break_str[1]), unlist(map(break_str, ~str_c(lt, ' ', .))))
  
  if (prefix_p){
    p_out <- str_c('p ', p_out)
    break_levels <- str_c('p ', break_levels)
  }
  
  return(factor(p_out, levels = break_levels))
}

blank_plot <- function(text = ''){
  ggplot(tibble(x=c(0, 1)), aes(x=x, y=x)) +
    geom_blank() +
    annotate(geom = 'text', x = 0.5, y = 0.5, label = text) +
    theme(panel.grid.major.y = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank())
}