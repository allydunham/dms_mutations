#!/usr/bin/env Rscript 
# Script to load and process deep mutagenesis study data for grouped analysis
# Fields each dataset needs:
## accession, species, common name, gene name, domain name, mut_id, position(s), variant(s)

library(Biostrings)
library(tidyverse)
library(magrittr)
library(readxl)

source('bin/functions.R')

deep_mut_data <- list()
formatted_deep_data <- list()
raw_seqs <- list()

#### Import Seqs ####
fasta_files <- dir('meta/seq', full.names = TRUE)
raw_seqs <- sapply(fasta_files, readAAStringSet)
names(raw_seqs) <- sapply(str_split(names(raw_seqs),'[\\/\\.]'), function(x){x[3]})

#### Hietpas 2011 Hsp90 ####
deep_mut_data$hietpas_2011_hsp90 <- read_csv('data/raw/processed/hietpas_2011_pdz_ligands_fitness.csv') %>%
  rename(alt_aa = aa) %>%
  mutate(species = 'saccharomyces_cerevisiae',
         name = 'hsp90',
         gene = 'hsp82',
         uniprot_acc = 'P02829',
         domain = NA,
         mut_id = gen_mut_id(uniprot_acc, NA, alt_aa, position))

df <- deep_mut_data$hietpas_2011_hsp90 %>%
  rename(score = selection_coefficient) %>%
  mutate(raw_score = score) %>%
  select(position, alt_aa, score, raw_score, codon) %>%
  mutate(ref_aa = as.vector(raw_seqs$s_cerevisiae_hsp82[[1]])[position],
         variants = str_c('p.', ref_aa, position, alt_aa)) %>%
  select(variants, score, raw_score, alt_codon = codon)

formatted_deep_data$hietpas_2011_hsp90 <- DeepMut(df, gene_name = 'Hsp90', species = 'Saccharomyces cerevisiae', domain = '-',
                                                  uniprot_id='P02829', authour='Hietpas et al.', year=2011,
                                                  transform = 'None', ref_seq = as.character(raw_seqs$s_cerevisiae_hsp82[[1]]),
                                                  misc=list(ensembl_gene_id='YPL240C', doi='https://doi.org/10.1073/pnas.1016024108',
                                                            url='http://www.pnas.org/content/108/19/7896', pubmed_id='21464309',
                                                            title='Experimental illumination of a fitness landscape', alt_name = 'Hsp82'))

#### Araya 2012 hYAP65 ####
deep_mut_data$araya_2012_hYAP65 <- read_tsv('data/raw/processed/araya_2012_hYAP65_ww.tsv', na = 'na',
                                            col_types = cols(positions=col_character())) %>%
  mutate(species = 'homo_sapiens',
         name = 'hYAP65',
         gene = 'yap65',
         domain = 'ww',
         uniprot_acc = 'P46937')

#### Starita 2013 Ube4b ####
deep_mut_data$starita_2013_ube4b <- read_xlsx('data/raw/processed/starita_2013_ube4b_ubox.xlsx', na = c('NA', ''))

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
         mut_id = gen_mut_id(uniprot_acc, NA, alt, position))

#### Jiang 2013 hsp90 ####
deep_mut_data$jiang_2013_hsp90 <- read_xlsx('data/raw/processed/jiang_2013_hsp90.xlsx', skip = 2) %>%
  select(-X__1) %>%
  rename_all(tolower) %>%
  rename(alt_aa = `amino acid`,
         sd = `standard deviation`) %>%
  mutate(uniprot_acc = 'P02829',
         mut_id = gen_mut_id(uniprot_acc, NA, alt_aa, position),
         average_num = as.numeric(average))

#### Forsyth 2013 igg ####
deep_mut_data$forsyth_2013_igg <- read_xlsx('data/raw/processed/forsyth_2013_igg_cdr.xlsx', na = 'NA') %>%
  rename(ref_codon = `WT codon`) %>%
  gather(key = 'alt_codon', value = 'enrichment', -position, -ref_codon, -distance) %>%
  separate(position, c('chain', 'position'), remove = TRUE) %>%
  mutate(ref_aa = str_sub(position, 1, 1),
         position = str_sub(position, 2))

#### Melamed 2013 pab1 ####
deep_mut_data$melamed_2013_pab1 <- read_xlsx('data/raw/processed/melamed_2013_pab1_rrm_enrichment_ratios.xlsx') %>%
  rename(ref_aa = WT_aa) %>%
  gather(key = 'alt_aa', value = 'enrichment_ratio', -position, -ref_aa) %>%
  mutate(species = 'saccharomyces_cerevisiae',
         uniprot_acc = 'P04147',
         mut_id = gen_mut_id(uniprot_acc, NA, alt_aa, position))

#### Wagenaar 2014 braf ####
## Only includes position/aa combos deemed significant
deep_mut_data$wagenaar_2014_braf <- read_xls('data/raw/processed/wagenaar_2014_braf.xls', skip = 3) %>%
  rename(position = Position,
         alt_aa = acid,
         median_enrichment = Median,
         rep1_codon1 = `Replicate 1`,
         rep1_codon2 = X__1,
         rep1_codon3 = X__2,
         rep1_codon4 = X__3,
         rep1_codon5 = X__4,
         rep1_codon6 = X__5,
         rep2_codon1 = `Replicate 2`,
         rep2_codon2 = X__6,
         rep2_codon3 = X__7,
         rep2_codon4 = X__8,
         rep2_codon5 = X__9,
         rep2_codon6 = X__10,
         ic50_vs_brafV600E = BRAFV600E,
         individually_tested = `mutant?`,
         possible_by_single_sub = `substitution?`) %>%
  filter(!is.na(rep1_codon1) & !rep1_codon1 == 'Replicate 1') %>%
  mutate_at(vars(-alt_aa, -individually_tested, -possible_by_single_sub, -ic50_vs_brafV600E), as.numeric) %>%
  mutate(species = 'homo_sapiens',
         uniprot_acc = 'P15056',
         mut_id = gen_mut_id(uniprot_acc, NA, alt_aa, position))

#### Firnberg 2014 tem1 ####
deep_mut_data$firnberg_2014_tem1 <- read_xlsx('data/raw/processed/firnberg_2014_tem1.xlsx') %>%
  rename(position = `Ambler Position`,
         ref_codon = `WT codon`,
         alt_codon = `Mutant codon`,
         ref_aa = `WT AA`,
         alt_aa = `Mutant AA`,
         base_changes = `Base Changes`,
         seq_counts_0.25 = `Sequencing Counts`,
         seq_counts_0.5 = X__1,
         seq_counts_1 = X__2,
         seq_counts_2 = X__3,
         seq_counts_4 = X__4,
         seq_counts_8 = X__5,
         seq_counts_16 = X__6,
         seq_counts_32 = X__7,
         seq_counts_64 = X__8,
         seq_counts_128 = X__9,
         seq_counts_256 = X__10,
         seq_counts_512 = X__11,
         seq_counts_1024 = X__12,
         total_seq_count = `Total Counts`,
         fitness = Fitness,
         fitness_err = `Estimated error in fitness`)
#### Roscoe 2014 Ubiquitin (E1 reactivity) ####
deep_mut_data$roscoe_2014_ubi_limiting_e1 <- read_xlsx('data/raw/processed/roscoe_2014_ubi_limiting_E1_reactivity.xlsx', skip = 3) %>%
  rename(position = Position,
         alt_aa = `Amino Acid`,
         log2_e1_react_vs_display = `log2 (E1react/display)`,
         rel_e1_reactivity = `Relative E1-reactivity (avg WT=1, avg STOP=0)`,
         sd_in_symonoymous_codons = `Standard deviation among synonymous codons`,
         notes = Notes)

deep_mut_data$roscoe_2014_ubi_excess_e1 <- read_xlsx('data/raw/processed/roscoe_2014_ubi_excess_E1_reactivity.xlsx', skip = 2) %>%
  rename(position = Position,
         alt_aa = `Amino Acid`,
         log2_e1_react_vs_display = `log2 (E1react/display)`,
         rel_e1_reactivity = `Relative Reactivity with Excess E1 (avg WT=1; avg STOP=0)`,
         sd_in_symonoymous_codons = `Standard deviation among synonymous codons`)

#### Melnikov 2014 APH(3')II ####
## Large folder of different conditions, needs more processing
melnikov_count_files <- dir('data/raw/processed/melnikov_2014_counts/') %>%
  grep('\\.aacounts\\.txt', ., value = TRUE)

# Function to read aa count tables from melnikov et al. 2014
read_melnikov_table <- function(fi){
  tbl <- read_tsv(paste0('data/raw/processed/melnikov_2014_counts/', fi), skip = 1, col_names = FALSE) %>%
    t() %>%
    set_rownames(NULL) %>%
    as_tibble() %>%
    set_colnames(.[1,]) %>%
    filter(!Position == 'Position') %>%
    rename(position = Position,
           ref_aa = `Wild-type`) %>%
    mutate_at(vars(-ref_aa), as.numeric)
  return(tbl)
}

melnikov_counts <- sapply(melnikov_count_files, read_melnikov_table, simplify = FALSE) %>%
  set_names(gsub('(KKA2\\_|\\.aacounts\\.txt)', '', names(.)))

melnikov_bkg_counts <- melnikov_counts[c('Bkg1', 'Bkg2')]

melnikov_counts <- melnikov_counts[which(!names(melnikov_counts) %in% c('Bkg1', 'Bkg2'))]

melnikov_e_scores <- mapply(melnikov_fitness, melnikov_counts, names(melnikov_counts),
                            MoreArgs = list(bkg=melnikov_bkg_counts), SIMPLIFY = FALSE) %>%
  bind_rows(.id = 'experiment') %>%
  separate(experiment, c('round', 'drug', 'library'), sep='_')

# Currently poorly normalised against wt performance so no comparable
deep_mut_data$melnikov_2014_aph3ii <- melnikov_e_scores

#### Findlay 2014 BRCA1 ####
deep_mut_data$findlay_2014_brca1_exon18 <- read_xlsx('data/raw/processed/findlay_2014_brca1_exon18_counts.xlsx', skip = 3) %>%
  rename(exon_position = `Exon Position`,
         alt_aa = Variant,
         average_effect_size = `Average effect size (both reps of L and R)`,
         mutPredSplice_score = `MutPredSplice score`,
         mutPredSplice_output = `MutPredSplice output`,
         mut_type = `Mutation Type`,
         lib_R1_rep1_effect_size = `Library R1 Replicate 1 effect size`,
         lib_R1_rep2_effect_size = `Library R1 Replicate 2 effect size`,
         lib_R_rep1_effect_size = `Library R Replicate 1 effect size`,
         lib_R_rep2_effect_size = `Library R Replicate 2 effect size`,
         lib_L_rep1_effect_size = `Library L Replicate 1 effect size`,
         lib_L_rep2_effect_size = `Library L Replicate 2 effect size`)

deep_mut_data$findlay_2014_brca1_hexamers <- read_xlsx('data/raw/processed/findlay_2014_brca1_hexamer_counts.xlsx', skip=3) %>%
  rename(hexamer = Hexamer,
         percent_in_input = `% in input library`,
         ESRseq_score = `Ke et al. 2011 ESRseq score`,
         nonsense = `nonsense hexamer?`,
         wt_hamming_dist = `hamming distance to WT`,
         log2_enrichment_score = `Log2 enrichment score`)

#### Findlay 2014 DBR1 ####
deep_mut_data$findlay_2014_dbr1 <- read_xlsx('data/raw/processed/findlay_2014_dbr1_exon2_counts.xlsx', skip = 3) %>%
  rename(seq = Sequence,
         log2_enrichment_score_day11_rep1 = `Day 11 log2 enrichment score (replicate 1)`,
         log2_enrichment_score_day11_rep2 = `Day 11 log2 enrichment score (replicate 2)`,
         position = `Affected codon`,
         ref_aa = `Reference AA`,
         alt_aa = `Substituted AA`,
         mut_type = `Variant Class`)

#### Olson 2014 Protein G ####
single_muts <- read_xlsx('data/raw/processed/olson_2014_protein_g_counts.xlsx', range = cell_limits(ul = c(3, 14), lr = c(NA, 18))) %>%
  rename(ref_aa1 = `WT amino acid`,
         pos1 = `Position`,
         alt_aa1 = `Mutation`,
         input_count = `Input Count`,
         selection_count = `Selection Count`)

double_muts <- read_xlsx('data/raw/processed/olson_2014_protein_g_counts.xlsx', range = cell_limits(ul = c(3, 2), lr = c(NA, 11))) %>%
  rename(ref_aa1 = `Mut1 WT amino acid`,
         pos1 = `Mut1 Position`,
         alt_aa1 = `Mut1 Mutation`,
         ref_aa2 = `Mut2 WT amino acid`,
         pos2 = `Mut2 Position`,
         alt_aa2 = `Mut2 Mutation`,
         input_count = `Input Count`,
         selection_count = `Selection Count`,
         fitness1 = `Mut1 Fitness`,
         fitness2 = `Mut2 Fitness`)

wt <- read_xlsx('data/raw/processed/olson_2014_protein_g_counts.xlsx', range = "U3:V4") %>%
  rename(input_count = `Input Count`,
         selection_count = `Selection Count`)

# Simplistic combining of data to get everying in, needs more complex processing
deep_mut_data$olson_2014_protein_g <- bind_rows(wt, single_muts, double_muts)

#### Starita 2015 Brca1 ####
deep_mut_data$starita_2015_brac1 <- read_xls('data/raw/processed/starita_2015_brca1_ring.xls', na = 'NA') %>%
  rename_all(tolower)

#### Kitzman 2015 Gal4 ####
kitzman_2015_path <- 'data/raw/processed/kitzman_2015_gal4_enrichment.xlsx'
read_kitzman_sheet <- function(sheet){
  tbl <- read_xlsx(kitzman_2015_path, skip = 1, na = 'ND', sheet = sheet) %>%
    rename(position = `Residue #`) %>%
    mutate(ref_aa = apply(., 1, function(x, nam){nam[x == 'wt' & !is.na(x)]}, nam = names(.)),
           label = sheet) %>%
    gather(key = 'alt_aa', value = 'log2_enrichment', -position, -ref_aa, -label) %>%
    mutate(log2_enrichment = if_else(log2_enrichment == 'wt', '0', log2_enrichment)) %>% # set wt to 0 log2 enrichment ratio
    mutate(log2_enrichment = as.numeric(log2_enrichment))
  return(tbl)
}

deep_mut_data$kitzman_2015_gal4 <- lapply(excel_sheets(kitzman_2015_path),
                                          read_kitzman_sheet) %>%
  bind_rows(.) %>%
  spread(key = 'label', value = 'log2_enrichment')

#### Mishra 2016 Hsp90 ####
## Each sheet also contains meta info that might be useful later
mishra_2016_path <- 'data/raw/processed/mishra_2016_hsp90_enrichment.xlsx'
read_mishra_sheet <- function(sheet){
  tbl <- read_xlsx(mishra_2016_path, sheet = sheet, col_names = FALSE)
  
  # Check sheet type
  if (tbl[1,1] == 'Stop counts'){
    ## Process sheets with a single replicate
    nom <- tbl[7,] %>% unlist(., use.names = FALSE)
    tbl <- tbl[8:length(tbl),] %>%
      set_names(nom) %>%
      rename_at(vars(-position, -aa), funs(paste0('rep1_', .))) %>%
      mutate_at(vars(-aa), as.numeric) %>%
      rename(alt_aa = aa) %>%
      mutate(avg_norm_ratiochange = rep1_norm_ratiochange)
    
  } else {
    ## Process sheets with replicates
    # Get first row of sub-tables
    top_row <- which(tbl$X__1 == 'position') + 1
    
    # Get bottom row of sub-tables
    bot_row <- sapply(top_row, find_next_empty_row, tbl=tbl) - 1
    
    # Extract sub-table names
    rep_nom <- tbl[top_row[1] - 1,] %>% unlist(., use.names = FALSE)
    ave_nom <- tbl[top_row[length(top_row)] - 1,] %>% unlist(., use.names = FALSE)
    ave_nom <- ave_nom[!is.na(ave_nom)]
    
    # Extract Subtables and add names
    rep1 <- tbl[top_row[1]:bot_row[1],] %>% 
      set_names(rep_nom) %>%
      rename_at(vars(-position, -aa), funs(paste0('rep1_', .)))
    
    rep2 <- tbl[top_row[2]:bot_row[2],] %>%
      set_names(rep_nom) %>%
      rename_at(vars(-position, -aa), funs(paste0('rep2_', .)))
    
    ave <- tbl[top_row[3]:bot_row[3],] %>%
      select_if(colSums(!is.na(.)) > 0) %>%
      set_names(ave_nom) %>%
      select(-s1, -s2) %>%
      rename(aa = `amino acid`)
    
    tbl <- full_join(rep1, rep2, by=c('position', 'aa')) %>%
      full_join(., ave, by=c('position', 'aa')) %>%
      mutate_at(vars(-aa), as.numeric) %>%
      rename(alt_aa = aa,
             avg_norm_ratiochange = avg)
  }
  return(tbl)
}

deep_mut_data$mishra_2016_hsp90 <- map(excel_sheets(mishra_2016_path), read_mishra_sheet) %>%
  bind_rows()

#### Sarkisyan 2016 GFP ####
# A file with just AAs is also available, but includes less info
deep_mut_data$sarkisyan_2016_gfp <- read_tsv('data/raw/processed/sarkisyan_2016_gfp_nucleotides.tsv')

#### Brenan 2016 Erk2 ####
deep_mut_data$brenan_2016_erk2 <- read_xlsx('data/raw/processed/brenan_2016_erk2.xlsx', sheet = 'Supplemental_Table_1') %>%
  rename_all(funs(gsub(' ', '_', tolower(.))))

#### Ashenberg 2016 Flu Nucleoprotein ####
deep_mut_data$ashenberg_2016_np <- read_csv('data/raw/processed/ashenberg_2017_flu_np.csv')

#### Weile 2017 ube2l ####
deep_mut_data$weile_2017_ube2l <- read_csv('data/raw/processed/weile_2017_ube2l_score_comp.csv', na = c('NA','','None'))

#### Weile 2017 sumo1 ####
deep_mut_data$weile_2017_sumo1 <- read_csv('data/raw/processed/weile_2017_sumo1_score_comp.csv', na = c('NA','','None'))

#### Weile 2017 tpk1 ####
deep_mut_data$weile_2017_tpk1 <- read_csv('data/raw/processed/weile_2017_tpk1_score_comp.csv', na = c('NA','','None'))

#### Weile 2017 calm ####
deep_mut_data$weile_2017_calm <- read_csv('data/raw/processed/weile_2017_calm_score_comp.csv', na = c('NA','','None'))

#### Findlay 2018 Brca1 ####
deep_mut_data$findlay_2018_brca1 <- read_xlsx('data/raw/processed/findlay_2018_brca1_ring_brct.xlsx', skip = 2, na = 'NA') %>%
  rename_all(funs(gsub('[\\/ \\(\\)]+', '_', .))) %>%
  rename(ref_nuc = reference,
         alt_nuc = alt,
         ref_aa = aa_ref,
         alt_aa = aa_alt)

#### Lee 2018 Flu haemagglutinin ####
# Xlsx also has sheets with preferences for each repeat, but importing just the average result for now
# Reports per site AA preference, which is fairly dissimilar to other metrics
deep_mut_data$lee_2018_flu_ha <- read_xlsx('data/raw/processed/lee_2018_influenza_ha.xlsx', sheet = 'avg_prefs') %>%
  rename(position = site) %>%
  gather(key = 'alt_aa', value = 'aa_pref', -position, -entropy, -neffective)

#### Giacomelli 2018 tp53 ####
deep_mut_data$giacomelli_2018_tp53 <- read_xlsx('data/raw/processed/giacomelli_2018_tp53.xlsx', skip=1) %>%
  rename_all(funs(tolower(gsub('[\\/ \\(\\)\\-]+', '_', gsub('\\*', '_2', .))))) %>%
  rename(ref_aa = aa_wt,
         alt_aa = aa_variant)

#### Save processed output ####
saveRDS(deep_mut_data, 'data/raw_deep_mut_data.RDS')
saveRDS(formatted_deep_data, 'data/formatted_deep_mut_data.RDS')

for (i in names(formatted_deep_data)){
  write_deep_mut(formatted_deep_data[[i]], str_c('data/standardised/', i, '.dm'))
}

dataset_size <- sapply(deep_mut_data, function(x){dim(x)[1]}) %>% unlist()
