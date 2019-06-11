#!/usr/bin/env Rscript 
# Script to load libraries and configurations for all mutations analyses, to keep consistency

# Load all libraries (keeps namespace consistent)
library(Biostrings)
library(tidyverse)
library(magrittr)
library(broom)
library(readxl)
library(ggpubr)
library(pdist)

# Source custom functions
source('bin/r_modules/utility_functions.R')
source('bin/r_modules/autoplot.R')
source('bin/r_modules/dm_functions.R')
source('bin/r_modules/variant_loading_functions.R')
source('bin/r_modules/variant_analysis_functions.R')
source('bin/r_modules/prediction_analysis_functions.R')

#### Project config ####
# Constants
MUT_SCORE_NAME <- 'Standardised Mutagenesis Score'
RAW_MUT_SCORE_NAME <- 'Raw Mutagenesis Score'

MUT_CATEGORIES <- list(deleterious='deleterious', neutral='neutral', beneficial='beneficial')

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
