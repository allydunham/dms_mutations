#!/usr/bin/env Rscript 
# Script to load and process deep mutagenesis study data for grouped analysis
library(tidyverse)
library(readxl)

deep_mut_data <- list()

#### Hietpas 2011 Hsp90 ####
deep_mut_data$hietpas_2011_hsp90 <- read_csv('data/raw/processed/hietpas_2011_pdz_ligands_fitness.csv') %>%
  mutate(species = 'saccharomyces_cerevisiae',
         gene = 'hsp90')
