#!/usr/bin/env Rscript 
# Script to transform Hietpas 2011 HSP90 per codon scoring into per AA scores by averaging over all codons

source('bin/config.R')

# DM File (future runs can use this and not have to correct)
dm <- read_deep_mut('data/standardised/hietpas_2011_hsp90/per_codon_scores/P02829_HSP90.dm')
dm$variant_data <- dm$variant_data %>%
  group_by(variants) %>%
  summarise(score = mean(score),
            raw_score = mean(raw_score)) %>%
  mutate(position = as.integer(str_sub(variants, start = 4, end = -2))) %>%
  arrange(position) %>%
  select(-position)

write_deep_mut(dm, 'data/standardised/hietpas_2011_hsp90/P02829_HSP90.dm')

# FoldX (only bother with averages and individual list files)
fx_head <- read_lines('data/standardised/hietpas_2011_hsp90/per_codon_scores/2CG9/Average_2CG9_Repair.fxout', n_max = 8)

for (pdb in str_split(dm$pdb_id, ':')){
  pdb_id <- pdb[1]
  pdb_offset <- ifelse(is.na(pdb[3]), 0, as.integer(pdb[3]))
  
  muts <- read_lines(str_c('data/standardised/hietpas_2011_hsp90/per_codon_scores/', pdb_id, '/individual_list_', pdb_id, '.txt')) %>%
    str_replace(., ';', '') %>%
    lapply(., format_pdb_variants, pdb_offset=pdb_offset) %>%
    unlist()
  
  fx <- read_tsv(str_c('data/standardised/hietpas_2011_hsp90/per_codon_scores/', pdb_id, '/Average_', pdb_id, '_Repair.fxout'), skip = 8,
                 col_types = cols(.default=col_double(), Pdb=col_character())) %>%
    mutate(Pdb=muts) %>%
    group_by(Pdb) %>%
    summarise_all(.funs = mean)
  
  # Write muts list
  write(fx$Pdb, str_c('data/standardised/hietpas_2011_hsp90/', pdb_id, '/individual_list_', pdb_id, '.txt'))
  
  # Write averages table
  fx <- mutate(fx, Pdb=str_c(pdb_id, '_Repair_', 1:dim(fx)[1]))
  fx_av_file <- str_c('data/standardised/hietpas_2011_hsp90/', pdb_id, '/Average_', pdb_id, '_Repair.fxout')
  write_lines(fx_head, fx_av_file)
  write_tsv(fx, fx_av_file, append = TRUE, col_names = TRUE)
}

# EVCouplings
ev <- read_csv('data/standardised/hietpas_2011_hsp90/per_codon_scores/ev/mutate/ev_dataset_predicted.csv',
               col_types = cols(mutant = col_character(), .default = col_double())) %>%
  group_by(mutant) %>%
  summarise_all(.funs = funs(mean)) %>%
  mutate(position = as.integer(str_sub(mutant, start = 2, end = -2))) %>%
  arrange(position) %>%
  select(-position)
write_csv(ev, 'data/standardised/hietpas_2011_hsp90/ev/mutate/ev_dataset_predicted.csv')
