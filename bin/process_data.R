#!/usr/bin/env Rscript 
# Script to load and process deep mutagenesis study data for grouped analysis
# Fields each dataset needs:
## accession, species, common name, gene name, mut_id, position(s), variant(s)

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
         mut_id = gen_mut_id(uniprot_acc, position, alt))

#### Jiang 2013 hsp90 ####
deep_mut_data$jiang_2013_hsp90 <- read_xlsx('data/raw/processed/jiang_2013_hsp90.xlsx', skip = 2) %>%
  select(-X__1)

#### Forsyth 2013 igg ####
deep_mut_data$forsyth_2013_igg <- read_tsv('data/raw/processed/forsyth_2013_igg_cdr.tsv') %>%
  rename(ref_codon = `WT codon`) %>%
  gather(key = 'alt_codon', value = enrichment, -position, -ref_codon, -distance)
#### Melamed 2013 pab1 ####
deep_mut_data$melamed_2013_pab1 <- read_xlsx('data/raw/processed/melamed_2013_pab1_rrm_enrichment_ratios.xlsx') %>%
  rename(ref_aa = WT_aa) %>%
  gather(key = 'alt_aa', value = 'enrichment_ratio', -position, -ref_aa)

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
         possible_by_single_sub = `substitution?`)

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
deep_mut_data$melnikov_2014_aph3ii <- NA

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

#### Starita 2015 Brca1 ####

#### Kitzman 2015 Gal4 ####

#### Mishra 2016 Hsp90 ####

#### Sarkisyan 2016 GFP ####

#### Brenan 2016 Erk2 ####

#### Ashenberg 2016 Flu Nucleoprotein ####

#### Weile 2017 ube2l ####

#### Weile 2017 sumo1 ####

#### Weile 2017 tpk1 ####

#### Weile 2017 calm ####

#### Findlay 2018 Brca1 ####

#### Lee 2018 Flu haemagglutinin ####

#### Giacomelli 2018 tp53 ####

#### Save processed output ####
saveRDS(deep_mut_data, 'data/deep_mut_data.RDS')

