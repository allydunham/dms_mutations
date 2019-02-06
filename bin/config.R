#!/usr/bin/env Rscript 
# Script to load libraries and configurations for all mutations analyses, to keep consistency

# Load all libraries (keeps namespace consistent)
library(Biostrings)
library(tidyverse)
library(magrittr)
library(broom)
library(readxl)
library(ggpubr)

# Source custom functions
source('bin/autoplot.R')
source('bin/dm_functions.R')
source('bin/variant_functions.R')

### Project configs ###
# Constants
MUT_SCORE_NAME <- 'Standardised Mutagenesis Score'
RAW_MUT_SCORE_NAME <- 'Raw Mutagenesis Score'

MUT_CATEGORIES <- list(deleterious='deleterious', neutral='neutral', beneficial='beneficial')

# ggplot theme (currently unset)
theme_set(theme_gray()) 
