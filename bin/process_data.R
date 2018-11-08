#!/usr/bin/env Rscript 
# Script to load and process deep mutagenesis study data for grouped analysis
library(tidyverse)
library(readxl)

source('bin/functions.R')

deep_mut_data <- list()

#### Hietpas 2011 Hsp90 ####
deep_mut_data$hietpas_2011_hsp90 <- read_csv('data/raw/processed/hietpas_2011_pdz_ligands_fitness.csv') %>%
  mutate(species = 'saccharomyces_cerevisiae',
         name = 'hsp90',
         gene = 'hsc82',
         uniprot_acc = 'P02829',
         mut_id = gen_mut_id(uniprot_acc, position, aa))

#### Araya 2012 hYAP65 ####
deep_mut_data$araya_2012_hYAP65 <- read_tsv('data/raw/processed/araya_2012_hYAP65_ww.tsv', na = 'na')

#### Roscoe 2013 Ubiquitin ####
deep_mut_data$roscoe_2013_ubi <- read_xlsx('data/raw/processed/roscoe_2013_ubi_fitness.xlsx', skip = 4) %>%
  rename(position = Position,
         alt = `Amino Acid`,
         selection_chr = Apparent,
         sd_chr = `Quantified Synonyms`) %>%
  mutate(selection_num = as.numeric(selection_chr),
         sd_num = as.numeric(sd_chr),
         name = 'ubiquitin',
         gene = 'ubc',
         uniprot_acc = 'P0CH08',
         species = 'saccaromyces_cerevisiae',
         mut_id = gen_mut_id(uniprot_acc, position, alt))
