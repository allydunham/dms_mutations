#!/usr/bin/env Rscript 
# Script to load processed deep mutagenesis data for analysis

source('src/config.R')
source('src/data_processing/variant_loading.R')

#### Load Data ####
root <- 'data/standardised'

deep_datasets <- str_c(root, '/', grep('*.dm', dir(root, recursive = TRUE), value = TRUE))

per_codon_datasets <- c('hietpas_2011_hsp90', 'weile_2017_sumo1', 'weile_2017_ube2i', 'findlay_2018_brca1', 'firnberg_2014_tem1')

per_codon <- sapply(deep_datasets, function(x){str_split(x, '/', simplify = TRUE)[1, 3] %in% per_codon_datasets}) %>% unname()

deep_variant_data <- mapply(import_dm_predictions_dataset, dm_file=deep_datasets, per_codon=per_codon)
names(deep_variant_data) <- sapply(deep_variant_data, function(x){make_dm_dataset_name(x$dm)})

# Save generated dataset
write_rds(deep_variant_data, 'data/rdata/processed_variant_data.RDS')

