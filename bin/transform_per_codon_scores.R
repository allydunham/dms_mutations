#!/usr/bin/env Rscript 
# Script to transform Hietpas 2011 HSP90 per codon scoring into per AA scores by averaging over all codons

source('bin/config.R')

# FoldX Header Lines
FX_HEAD <- read_lines('data/standardised/hietpas_2011_hsp90/per_codon_scores/2CG9/Average_2CG9_Repair.fxout', n_max = 8)

# Average duplicate variants and write dm File
adjust_dm_file <- function(dm, path){
  dm$variant_data <- dm$variant_data %>%
  group_by(variants) %>%
  summarise(score = mean(score),
            raw_score = mean(raw_score)) %>%
  mutate(position = as.integer(str_sub(variants, start = 4, end = -2))) %>%
  arrange(position) %>%
  select(-position)

  write_deep_mut(dm, path)
}

# Average mutant runs and adjust individual list
# pdb is a vec of c(pdb_id, chain, (offset), (subset))
adjust_foldx <- function(pdb, root){
  pdb_id <- pdb[1]
  pdb_chain <- pdb[2]
  pdb_offset <- ifelse(is.na(pdb[3]), 0, as.integer(pdb[3]))
  
  muts <- read_lines(str_c(root, '/per_codon_scores/', pdb_id, '/individual_list_', pdb_id, '.txt'))
  
  fx <- read_tsv(str_c(root, '/per_codon_scores/', pdb_id, '/Average_', pdb_id, '_Repair.fxout'), skip = 8,
                 col_types = cols(.default=col_double(), Pdb=col_character())) %>%
    mutate(Pdb=muts) %>%
    group_by(Pdb) %>%
    summarise_all(.funs = mean) %>%
    mutate(position = as.integer(str_sub(Pdb, start = 3, end = -3))) %>%
    arrange(position) %>%
    select(-position)
  
  # Write muts list
  write_lines(fx$Pdb, str_c(root, '/', pdb_id, '/individual_list_', pdb_id, '.txt'))
  
  # Write averages table
  fx <- mutate(fx, Pdb=str_c(pdb_id, '_Repair_', 1:dim(fx)[1]))
  fx_av_file <- str_c(root, '/', pdb_id, '/Average_', pdb_id, '_Repair.fxout')
  write_lines(FX_HEAD, fx_av_file)
  write_tsv(fx, fx_av_file, append = TRUE, col_names = TRUE)
}

# Adjust EVCouplings output file
adjust_ev <- function(inpath, outpath){
  ev <- read_csv(inpath,
                 col_types = cols(mutant = col_character(), .default = col_double())) %>%
    group_by(mutant) %>%
    summarise_all(.funs = funs(mean)) %>%
    mutate(position = as.integer(str_sub(mutant, start = 2, end = -2))) %>%
    arrange(position) %>%
    select(-position)
  write_csv(ev, outpath)
}

#### Hietpas 2011 hsp90 ####
# DM File (future runs can use this and not have to correct)
hietpas_dm <- read_deep_mut('data/standardised/hietpas_2011_hsp90/per_codon_scores/P02829_HSP90.dm')
adjust_dm_file(hietpas_dm, 'data/standardised/hietpas_2011_hsp90/P02829_HSP90.dm')

# FoldX (only bother with averages and individual list files)

for (pdb in str_split(hietpas_dm$pdb_id, ':')){
  adjust_foldx(pdb, 'data/standardised/hietpas_2011_hsp90')
}

# EVCouplings
adjust_ev('data/standardised/hietpas_2011_hsp90/per_codon_scores/ev/mutate/ev_dataset_predicted.csv',
          'data/standardised/hietpas_2011_hsp90/ev/mutate/ev_dataset_predicted.csv')


#### Weile 2017 SUMO1 ####
# DM File (future runs can use this and not have to correct)
weile_sumo1_dm <- read_deep_mut('data/standardised/weile_2017_sumo1/per_codon_scores/P63165_SUMO1.dm')
adjust_dm_file(weile_sumo1_dm, 'data/standardised/weile_2017_sumo1/P63165_SUMO1.dm')

# FoldX
for (pdb in str_split(weile_sumo1_dm$pdb_id, ':')){
  adjust_foldx(pdb, 'data/standardised/weile_2017_sumo1')
}

# EVCouplings
adjust_ev('data/standardised/weile_2017_sumo1/per_codon_scores/ev/mutate/ev_dataset_predicted.csv',
          'data/standardised/weile_2017_sumo1/ev/mutate/ev_dataset_predicted.csv')

#### Weile 2017 Ube2i ####
weile_ube2i_dm <- read_deep_mut('data/standardised/weile_2017_ube2i/per_codon_scores/P63279_UBE2I.dm')
adjust_dm_file(weile_ube2i_dm, 'data/standardised/weile_2017_ube2i/P63279_UBE2I.dm')

# FoldX
for (pdb in str_split(weile_ube2i_dm$pdb_id, ':')){
  adjust_foldx(pdb, 'data/standardised/weile_2017_ube2i')
}

# EVCouplings
adjust_ev('data/standardised/weile_2017_ube2i/per_codon_scores/ev/mutate/ev_dataset_predicted.csv',
          'data/standardised/weile_2017_ube2i/ev/mutate/ev_dataset_predicted.csv')
