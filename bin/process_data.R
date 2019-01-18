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

species_options <- list(cerevisiae='Saccharomyces cerevisiae', sapiens='Homo sapiens', musculus='Mus musculus', coli='Escherichia coli',
                        strep='Streptococcus', flu='Influenza')

#### Import Seqs ####
fasta_files <- dir('meta/seq', full.names = TRUE)
raw_seqs <- sapply(fasta_files, readAAStringSet) %>% sapply(., function(x){as.character(as.vector(x[[1]]))})
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
  mutate(ref_aa = raw_seqs$s_cerevisiae_hsp82[position],
         variants = str_c('p.', ref_aa, position, alt_aa)) %>%
  select(variants, score, raw_score, alt_codon = codon)

formatted_deep_data$hietpas_2011_hsp90 <- DeepMut(df, gene_name = 'HSP90', species = species_options$cerevisiae,
                                                  uniprot_id='P02829', authour='Hietpas et al.', pdb_id = c('2CG9:A', '2CGE:A'),
                                                  transform = 'None', aa_seq = str_c(raw_seqs$s_cerevisiae_hsp82, collapse = ''),
                                                  ensembl_gene_id='YPL240C', doi='10.1073/pnas.1016024108', year=2011,
                                                  url='http://www.pnas.org/content/108/19/7896', pmid='21464309',
                                                  title='Experimental illumination of a fitness landscape', alt_name = 'HSP82')

#### Araya 2012 hYAP65 ####
deep_mut_data$araya_2012_hYAP65 <- read_tsv('data/raw/processed/araya_2012_hYAP65_ww.tsv', na = 'na',
                                            col_types = cols(positions=col_character())) %>%
  mutate(species = 'homo_sapiens',
         name = 'hYAP65',
         gene = 'yap65',
         domain = 'ww',
         uniprot_acc = 'P46937')

# extract ww ref seq
df <- deep_mut_data$araya_2012_hYAP65 %>%
  select(slope, fitness, rsquared, mutations, positions, amino.acids) %>%
  mutate(mutations = str_split(mutations, ','),
         alt_aas = str_split(amino.acids, ','),
         # position from Araya et al. appears to be offset by 9 compared to their seq, so add 160 not 170 to reach domain start
         positions = sapply(str_split(positions, ','), function(x){as.numeric(x) + 160}),
         ref_aas = sapply(positions, function(x){raw_seqs$h_sapiens_yap1[x]}),
         variants = mapply(function(pos, ref, alt){str_c('p.', ref, pos, alt, collapse = ',')},
                           positions, ref_aas, alt_aas),
         score = log2(fitness)) %>%
  rename(raw_score = fitness) %>%
  select(variants, score, raw_score, slope, rsquared)

formatted_deep_data$araya_2012_hYAP65 <- DeepMut(variant_data = df, gene_name = 'YAP1', domain = 'WW', species = species_options$sapiens,
                                                 aa_seq = str_c(raw_seqs$h_sapiens_yap1, collapse = ''), transform = 'log2', uniprot_id = 'P46937',
                                                 authour = 'Araya at al.', year = 2012, pdb_id = c('4REX:A', '1JMQ:A:160'),
                                                 alt_name = 'hYAP65', title='A fundamental protein property, thermodynamic stability, revealed solely from large-scale measurements of protein function',
                                                 url='http://www.pnas.org/content/109/42/16858', doi='10.1073/pnas.1209751109',
                                                 pmid='23035249', note='Authours WW domain differs slightly from uniprot, 170-203 vs 171-204')

#### Starita 2013 Ube4b ####
deep_mut_data$starita_2013_ube4b <- read_xlsx('data/raw/processed/starita_2013_ube4b_ubox.xlsx', na = c('NA', ''))

df <- deep_mut_data$starita_2013_ube4b %>%
  rename(score = nscor_log2_ratio) %>%
  mutate(raw_score = score) %>%
  separate(seqID, into = c('position', 'alt_aa'), sep='-') %>%
  mutate(position = str_split(position, ','),
         position = sapply(position, function(x){as.numeric(x) + 1072})) %>% # tested region starts at +1072 according to Starita (slightly before uniprot UBOX) This does lead to ref seq aligning
  filter(sapply(position, function(x){!any(is.na(x))})) %>%
  mutate(alt_aa = str_split(alt_aa, ','),
         ref_aa = sapply(position, function(x){raw_seqs$m_musculus_ube4b[x]}),  
         variants = mapply(function(pos, ref, alt){str_c('p.', ref, pos, alt, collapse = ',')}, position, ref_aa, alt_aa)) %>%
  select(variants, score, raw_score, log2_ratio) %>%
  rename(score_without_nscor_adjust = log2_ratio)

formatted_deep_data$starita_2013_ube4b <- DeepMut(variant_data = df, gene_name = 'UBE4B', domain = 'UBOX', species = species_options$musculus,
                                                  aa_seq = str_c(raw_seqs$m_musculus_ube4b, collapse = ''), transform = 'None', uniprot_id = 'Q9ES00',
                                                  authour = 'Starita et al.', year = 2013,
                                                  title='Activity-enhancing mutations in an E3 ubiquitin ligaseidentified by high-throughput mutagenesis',
                                                  url='http://www.pnas.org/content/110/14/E1263.long', doi='10.1073/pnas.1303309110', pmid='23509263')

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

df <- deep_mut_data$roscoe_2013_ubi %>%
  mutate(ref = raw_seqs$s_cerevisiae_ubi1[position],
         variants = str_c('p.', ref, position, alt),
         raw_score = selection_num) %>%
  rename(score = selection_num,
         score_chr = selection_chr) %>%
  select(variants, score, raw_score, sd_num, score_chr, sd_chr)

formatted_deep_data$roscoe_2013_ubi <- DeepMut(variant_data = df, gene_name = 'UBI4', species = species_options$cerevisiae, uniprot_id = 'P0CG63',
                                               aa_seq = str_c(raw_seqs$s_cerevisiae_ubi4, collapse = ''), authour = 'Roscoe et al.', year = 2013,
                                               alt_name = 'Ubiquitin', pmid = '23376099', pdb_id=c('3CMM:B', '3OLM:D'),
                                               title = 'Analyses of the effects of all ubiquitin point mutants on yeast growth rate',
                                               url='https://www.sciencedirect.com/science/article/pii/S0022283613000636',
                                               doi='10.1016/j.jmb.2013.01.032')

#### Jiang 2013 hsp90 ####
deep_mut_data$jiang_2013_hsp90 <- read_xlsx('data/raw/processed/jiang_2013_hsp90.xlsx', skip = 2) %>%
  select(-X__1) %>%
  rename_all(tolower) %>%
  rename(alt_aa = `amino acid`,
         sd = `standard deviation`) %>%
  mutate(uniprot_acc = 'P02829',
         mut_id = gen_mut_id(uniprot_acc, NA, alt_aa, position),
         average_num = as.numeric(average))

df <- deep_mut_data$jiang_2013_hsp90 %>%
  mutate(ref = raw_seqs$s_cerevisiae_hsp82[position],
         variants = str_c('p.',ref,position,alt_aa),
         score = log2(average_num)) %>%
  rename(raw_score = average_num) %>%
  select(variants, score, raw_score, sd, gpd, tef, tefdter, cyc, adh, cycdter, adhdter)

formatted_deep_data$jiang_2013_hsp90 <- DeepMut(variant_data = df, gene_name = 'HSP90', domain = 'Putative substrate binding loop',
                                                species = species_options$cerevisiae, transform = 'log2', uniprot_id = 'P02829',
                                                aa_seq = str_c(raw_seqs$s_cerevisiae_hsp82, collapse = ''), authour = 'Jiang et al.', year = 2013,
                                                alt_name='HSP82', url='https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1003600',
                                                doi='10.1371/journal.pgen.1003600', pmid='23825969', pdb_id = c('2CG9:A', '2CGE:A'),
                                                title='Latent Effects of Hsp90 Mutants Revealed at Reduced Expression Levels')

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

df <- deep_mut_data$melamed_2013_pab1 %>%
  mutate(variants = str_c('p.', ref_aa, position, alt_aa),
         score = enrichment_ratio) %>%
  rename(raw_score = enrichment_ratio) %>%
  select(variants, score, raw_score)

formatted_deep_data$melamed_2013_pab1 <- DeepMut(variant_data = df, gene_name = 'PAB1', domain = 'RRM', species = species_options$cerevisiae,
                                                 aa_seq = str_c(raw_seqs$s_cerevisiae_pab1, collapse = ''), uniprot_id = 'P04147',
                                                 authour = 'Melamed et al.', year = 2013, url='https://rnajournal.cshlp.org/content/19/11/1537',
                                                 doi='10.1261/rna.040709.113', pmid='24064791', 
                                                 title='Deep mutational scanning of an RRM domain of the Saccharomyces cerevisiae poly(A)-binding protein')

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

df <- deep_mut_data$wagenaar_2014_braf %>%
  mutate(ref = raw_seqs$h_sapiens_braf_v600e[position],
         variants = str_c('p.', ref, position, alt_aa),
         score = log2(median_enrichment),
         raw_score = median_enrichment) %>%
  select(variants, score, raw_score, ic50_vs_brafV600E, individually_tested, possible_by_single_sub)

formatted_deep_data$wagenaar_2014_braf <- DeepMut(variant_data = df, gene_name = 'BRAF', species = species_options$sapiens,
                                                  aa_seq = str_c(raw_seqs$h_sapiens_braf_v600e, collapse = ''), transform = 'log2', uniprot_id = 'P15056',
                                                  authour = 'Wagenaar et al.', year = 2014, pdb_id=c('1UWH:A', '3C4C:A'),
                                                  notes='Only retained variants they deemed significantly different from wt',
                                                  title='Resistance to vemurafenib resulting from a novel mutation in the BRAFV600E kinase domain',
                                                  pmid='24112705', doi='10.1111/pcmr.12171', url='https://onlinelibrary.wiley.com/doi/full/10.1111/pcmr.12171')

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

df <- deep_mut_data$firnberg_2014_tem1 %>%
  filter(!ref_aa == '*') %>%
  mutate(position = rep(1:length(raw_seqs$e_coli_tem1), each=64), # Given pos are weirdly numbered, but match WT seq exactly with all codons (64) at each pos
         variants = str_c('p.', ref_aa, position, alt_aa),
         score=log2(fitness),
         genomic_variants = str_c('c.', (position-3)*3 + 1, '_', (position-2)*3, 'del', ref_codon,  'ins', alt_codon)) %>%
  select(variants, score, raw_score=fitness, raw_score_err=fitness_err, genomic_variants, base_changes, total_seq_count) %>%
  filter(!is.na(variants))

formatted_deep_data$firnberg_2014_tem1 <- DeepMut(variant_data = df, gene_name = 'TEM1', species = species_options$coli,
                                                  aa_seq = str_c(raw_seqs$e_coli_tem1, collapse = ''), transform = 'log2', uniprot_id = 'Q6SJ61',
                                                  authour = 'Firnberg et al.', year = 2014,
                                                  title='A Comprehensive, High-Resolution Map of a Gene’s Fitness Landscape',
                                                  doi='10.1093/molbev/msu081', pmid='24567513',
                                                  url='https://academic.oup.com/mbe/article/31/6/1581/2925654')

#### Roscoe 2014 Ubiquitin (E1 reactivity) ####
deep_mut_data$roscoe_2014_ubi_limiting_e1 <- read_xlsx('data/raw/processed/roscoe_2014_ubi_limiting_E1_reactivity.xlsx', skip = 3) %>%
  rename(position = Position,
         alt_aa = `Amino Acid`,
         log2_e1_react_vs_display = `log2 (E1react/display)`,
         rel_e1_reactivity = `Relative E1-reactivity (avg WT=1, avg STOP=0)`,
         sd_in_symonoymous_codons = `Standard deviation among synonymous codons`,
         notes = Notes) %>%
  mutate_at(.vars = vars(log2_e1_react_vs_display, rel_e1_reactivity, sd_in_symonoymous_codons), as.numeric)

deep_mut_data$roscoe_2014_ubi_excess_e1 <- read_xlsx('data/raw/processed/roscoe_2014_ubi_excess_E1_reactivity.xlsx', skip = 2) %>%
  rename(position = Position,
         alt_aa = `Amino Acid`,
         log2_e1_react_vs_display = `log2 (E1react/display)`,
         rel_e1_reactivity = `Relative Reactivity with Excess E1 (avg WT=1; avg STOP=0)`,
         sd_in_symonoymous_codons = `Standard deviation among synonymous codons`) %>%
  mutate_at(.vars = vars(log2_e1_react_vs_display, rel_e1_reactivity, sd_in_symonoymous_codons), as.numeric)

df1 <- deep_mut_data$roscoe_2014_ubi_limiting_e1 %>%
  mutate(ref_aa = raw_seqs$s_cerevisiae_ubi1[position],
         variants = str_c('p.', ref_aa, position, alt_aa),
         raw_score = log2_e1_react_vs_display) %>%
  select(variants, score = log2_e1_react_vs_display, raw_score, rel_e1_reactivity, sd_in_symonoymous_codons, notes)

df2 <- deep_mut_data$roscoe_2014_ubi_excess_e1 %>%
  mutate(ref_aa = raw_seqs$s_cerevisiae_ubi1[position],
         variants = str_c('p.', ref_aa, position, alt_aa),
         raw_score = log2_e1_react_vs_display) %>%
  select(variants, score = log2_e1_react_vs_display, raw_score, rel_e1_reactivity, sd_in_symonoymous_codons)


formatted_deep_data$roscoe_2014_ubi <- DeepMutSet(list(limiting_e1=DeepMut(variant_data = df1, gene_name = 'UBI4', species = species_options$cerevisiae,
                                                                           aa_seq = str_c(raw_seqs$s_cerevisiae_ubi4, collapse = ''), transform = 'None',
                                                                           uniprot_id = 'P0CG63', authour = 'Roscoe and Bolon', year = 2014,
                                                                           pdb_id=c('3CMM:B', '3OLM:D'),
                                                                           alt_name='Ubiquitin', pmid='24862281', doi='10.1016/j.jmb.2014.05.019',
                                                                           url='https://www.sciencedirect.com/science/article/pii/S0022283614002587?via%3Dihub',
                                                                           title='Systematic Exploration of Ubiquitin Sequence, E1 Activation Efficiency, and Experimental Fitness in Yeast'),
                                                       excess_e1=DeepMut(variant_data = df2, gene_name = 'UBI4', species = species_options$cerevisiae,
                                                                         aa_seq = str_c(raw_seqs$s_cerevisiae_ubi4, collapse = ''), transform = 'None',
                                                                         uniprot_id = 'P0CG63', authour = 'Roscoe and Bolon', year = 2014,
                                                                         pdb_id=c('3CMM:B', '3OLM:D'),
                                                                         alt_name='Ubiquitin', pmid='24862281', doi='10.1016/j.jmb.2014.05.019',
                                                                         url='https://www.sciencedirect.com/science/article/pii/S0022283614002587?via%3Dihub',
                                                                         title='Systematic Exploration of Ubiquitin Sequence, E1 Activation Efficiency, and Experimental Fitness in Yeast')))

#### Melnikov 2014 APH(3')II ####
## Large folder of different conditions, needs more processing
melnikov_count_files <- dir('data/raw/processed/melnikov_2014_counts/') %>%
  grep('\\.aacounts\\.txt', ., value = TRUE)

melnikov_count_files <- melnikov_count_files[!grepl('(S[12]\\_Ami|S3\\_Kan)', melnikov_count_files)]

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

  
# Currently poorly normalised against wt performance so not necessarily comparable
deep_mut_data$melnikov_2014_aph3ii <- melnikov_e_scores

# TODO - deal with the different library/round/drug data in a more satisfactory manner
# Round and library contain the same information (plus round notes which needed a re-test at different MIC)
#  -> discard round since we hae already selected the right MIC experiments
df <- deep_mut_data$melnikov_2014_aph3ii %>%
  mutate(variants = str_c('p.', ref_aa, position, alt_aa),
         score=log2(e_score)) %>%
  select(variants, score, raw_score='e_score', drug, library) %>%
  select(-raw_score) %>%
  spread(key = 'library', value = 'score') %>%
  mutate(diff = abs(L1 - L2))

# Broad correlation between the two libraries, discard outliers and take mean
#p_melnikov_lib_test <- ggplot(df, aes(x=L1, y=L2, col=drug)) + geom_point() 
#ggsave('figures/initial_analysis/melnikov_library_correlation.pdf', p_melnikov_lib_test, width = 8, height = 6)

max_diff <- sd(df$diff, na.rm = TRUE) * 3

df %<>% filter(diff < max_diff) %>%
  mutate(score = (L1 + L2)/2) %>%
  drop_na(score) %>%
  select(variants, drug, score) %>% 
  mutate(rel_conc = 1/as.numeric(str_sub(drug, -1)),
         drug = str_sub(drug, 1, -3))

## Currently take the unweighted mean across drug concs, which might not be the most scientific method
# although will give a measure that relates severity of variant
df_mean <- group_by(df, variants, drug) %>%
  summarise(mscore = mean(score)) %>%
  full_join(., df, by = c("variants", "drug")) %>%
  spread(key = rel_conc, value = score, sep = '_') %>%
  mutate(raw_score = rel_conc_1) %>%
  rename(score=mscore) %>%
  ungroup()

lst <- sapply(unique(df_mean$drug),
              function(x){DeepMut(variant_data = filter(df_mean, drug == x) %>% select(variants, score, raw_score, rel_conc_0.125, rel_conc_0.25, rel_conc_0.5, rel_conc_1),
                                  gene_name = "APH3-II", alt_name="APH(3')II", species = species_options$coli,
                                  aa_seq = str_c(raw_seqs$e_coli_aph_3prime_II, collapse = ''),
                                  transform = 'mean of different drug concs (raw is at 1:1 MIC)',
                                  uniprot_id = 'Q58HT3', authour = 'Melnikov et al.', year = 2014,
                                  title='Comprehensive mutational scanning of a kinase in vivo reveals substrate-dependent fitness landscapes',
                                  doi='10.1093/nar/gku511', pmid='24914046',
                                  url='https://academic.oup.com/nar/article/42/14/e112/1266940',
                                  notes='gene_name is transformed to be computer readable, main name in alt_name field. Concs are given relative to esstimated drug minimum inhibitory concentration (MIC)',
                                  drug=x)},
              simplify = FALSE)

formatted_deep_data$melnikov_2014_aph3ii <- DeepMutSet(lst)

#### Findlay 2014 BRCA1 ####
deep_mut_data$findlay_2014_brca1_exon18 <- read_xlsx('data/raw/processed/findlay_2014_brca1_exon18_counts.xlsx', skip = 3) %>%
  rename(exon_position = `Exon Position`,
         alt = Variant,
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

s_dna <- deep_mut_data$findlay_2014_brca1_exon18 %>%
  group_by(exon_position) %>%
  summarise(ref = Biostrings::DNA_BASES[!Biostrings::DNA_BASES %in% alt]) %>%
  pull(ref)

# Found Exon18 location using EBI Water tool, based on exon18 being in +1 frame and aligning translated sequence below
# s_prot <- translate(DNAString(str_c('G', str_c(s_dna, collapse = ''), 'GG')))
# s_prot_leader <- raw_seqs$h_sapiens_brca1[[1]][1:1691] # aa seq before 
# s_prot_follower <- raw_seqs$h_sapiens_brca1[[1]][1719:1863]

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

E_wt <- wt$selection_count/wt$input_count

# Simplistic combining of data to get everying in, needs more complex processing
deep_mut_data$olson_2014_protein_g <- bind_rows(single_muts, double_muts) %>%
  mutate(er = selection_count/input_count,
         f = er/E_wt,
         score = log2(f))

df <- deep_mut_data$olson_2014_protein_g %>% # Position is offset 1 from uniprot ref seq
  mutate(variants = ifelse(is.na(pos2), str_c('p.', ref_aa1, pos1-1, alt_aa1),str_c('p.', ref_aa1, pos1-1, alt_aa1,',p.', ref_aa2, pos2-1, alt_aa2))) %>%
  select(variants, score, raw_score=er)
  
formatted_deep_data$olson_2014_protein_g <- DeepMut(variant_data = df, gene_name = 'ProteinG', species = species_options$strep,
                                                    aa_seq = str_c(raw_seqs$strep_protein_g, collapse = ''), transform = 'log2',
                                                    uniprot_id = 'P19909', authour = 'Olson et al.', year = 2014,
                                                    notes='Uniprot ID is for the Protein G precurssor which has a slightly different sequence',
                                                    doi='10.1016/j.cub.2014.09.072', url='https://www.sciencedirect.com/science/article/pii/S0960982214012688',
                                                    pmid='25455030', title='A Comprehensive Biophysical Description of Pairwise Epistasis throughout an Entire Protein Domain')

#### Starita 2015 Brca1 ####
deep_mut_data$starita_2015_brca1 <- read_xls('data/raw/processed/starita_2015_brca1_ring.xls', na = 'NA') %>%
  rename_all(tolower) %>%
  mutate(ref = raw_seqs$h_sapiens_brca1[pos], # Ref seq given by study has a mysterious, undocumented R at pos 175 where normal refs have K, using K here since the change is not explained in the paper and appears erroneous
         variants = str_c('p.',ref,pos,mut))

# TODO Better ways of making scores initially non negative
df_y2h <- deep_mut_data$starita_2015_brca1 %>%
  filter(!variant_id == 'NA-NA') %>%
  mutate(score = log2((y2h_score + 1)/2)) %>%
  select(variants, score, raw_score = y2h_score)

df_e3 <- deep_mut_data$starita_2015_brca1 %>%
  filter(!variant_id == 'NA-NA') %>%
  mutate(tmp = if_else(e3_score > 0, e3_score, min(abs(e3_score), na.rm = TRUE)),
         score = log2(tmp)) %>%
  select(variants, score, raw_score = e3_score)

formatted_deep_data$starita_2015_brca1 <- DeepMutSet(list(
  bard1_ring_binding = DeepMut(df_y2h, gene_name = 'BRCA1', domain = 'RING', species = species_options$sapiens,
                               aa_seq = str_c(raw_seqs$h_sapiens_brca1, collapse = ''), transform = 'log2',
                               uniprot_id = 'P38398', authour = 'Starita et al.', year = 2015,
                               notes='Y2H experiment measured binding to BARD1 RING domain',
                               title='Massively Parallel Functional Analysis of BRCA1 RING Domain Variants',
                               url='http://www.genetics.org/content/200/2/413',
                               doi='10.1534/genetics.115.175802', pmid='25823446'),
  e3_activity = DeepMut(variant_data = df_e3, gene_name = 'BRCA1', domain = 'RING', species = species_options$sapiens,
                        aa_seq = str_c(raw_seqs$h_sapiens_brca1, collapse = ''), transform = 'log2(set < 0 to min > 0)',
                        uniprot_id = 'P38398', authour = 'Starita et al.', year = 2015,
                        notes='Phage display assay measures E3 ubiquitin ligase activity',
                        title='Massively Parallel Functional Analysis of BRCA1 RING Domain Variants',
                        url='http://www.genetics.org/content/200/2/413',
                        doi='10.1534/genetics.115.175802', pmid='25823446')))

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
  spread(key = 'label', value = 'log2_enrichment') %>%
  mutate(score = (SEL_A_40h + SEL_B_40h + SEL_C_40h)/3,
         raw_score=score,
         variants = str_c('p.', ref_aa, position, alt_aa)) %>%
  filter(!alt_aa == 'delInFrame')

formatted_deep_data$kitzman_2015_gal4 <- DeepMut(variant_data = select(deep_mut_data$kitzman_2015_gal4, variants, score, raw_score,
                                                                       NONSEL_24h, SEL_A_24h, SEL_A_40h, SEL_B_40h, SEL_C_40h, SEL_C_64h),
                                                 transform = 'Log2, Average of A,B,C at 40h', gene_name = 'GAL4', species = species_options$cerevisiae,
                                                 authour = 'Kitzman et al.', year = 2015, aa_seq = str_c(raw_seqs$s_cerevisiae_gal4, collapse = ''),
                                                 uniprot_id = 'P04386', pdb_id = c('1D66:A'),
                                                 title='Massively parallel single amino-acid mutagenesis', pmid='25559584',
                                                 doi='10.1038/nmeth.3223', url='https://www.nature.com/articles/nmeth.3223')

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

df <- deep_mut_data$mishra_2016_hsp90 %>%
  mutate(ref_aa = raw_seqs$s_cerevisiae_hsp82[position],
         variants = str_c('p.', ref_aa, position, alt_aa),
         score = avg_norm_ratiochange,
         raw_score = score)

formatted_deep_data$mishra_2016_hsp90 <- DeepMut(variant_data = df, gene_name = 'HSP90', species = species_options$cerevisiae,
                                                  authour = 'Mishra et al.', year = 2016, aa_seq = str_c(raw_seqs$s_cerevisiae_hsp82, collapse = ''),
                                                 transform = 'None', uniprot_id = 'P02829', pdb_id = c('2CG9:A', '2CGE:A'),
                                                 alt_name='HSP82', title='Systematic Mutant Analyses Elucidate General and Client-Specific Aspects of Hsp90 Function',
                                                 doi='10.1016/j.celrep.2016.03.046', pmid='27068472',
                                                 url='https://www.sciencedirect.com/science/article/pii/S2211124716303175')

#### Sarkisyan 2016 GFP ####
# A file with just AAs is also available, but includes less info
tmp <- read_tsv('data/raw/processed/sarkisyan_2016_gfp_AAs.tsv')
deep_mut_data$sarkisyan_2016_gfp <- read_tsv('data/raw/processed/sarkisyan_2016_gfp_nucleotides.tsv')

#### Brenan 2016 Erk2 ####
deep_mut_data$brenan_2016_erk2 <- read_xlsx('data/raw/processed/brenan_2016_erk2.xlsx', sheet = 'Supplemental_Table_1') %>%
  rename_all(funs(gsub(' ', '_', tolower(.))))

df <- deep_mut_data$brenan_2016_erk2 %>%
  mutate(variants = str_c('p.', erk2_mutant),
         score = `lfc_(etp_vs._dox)`) %>%
  select(variants, score, raw_score=`lfc_(etp_vs._dox)`, nuc_acid_changes, lfc_etp_vs_sch=`lfc_(etp_vs._sch)`, lfc_etp_vs_vrt=`lfc_(etp_vs._vrt)`,
         dox_rank, sch_rank, vrt_rank, vrt_specific_allele, sch_specific_allele)

formatted_deep_data$brenan_2016_erk2 <- DeepMut(variant_data = df, gene_name = 'ERK2', species = species_options$sapiens, transform = 'None',
                                                authour = 'Brenan et al. 2016', year = 2016, aa_seq = str_c(raw_seqs$h_sapiens_mapk1, collapse = ''),
                                                uniprot_id = 'P28482', pdb_id = c('1PME:A', '1TVO:A'),
                                                alt_name='MAPK1', doi='10.1016/j.celrep.2016.09.061', pmid='27760319',
                                                url='https://www.sciencedirect.com/science/article/pii/S2211124716313171',
                                                notes='Also includes scores in two other drug conditions. Used A375 cells with BRAFV600E.',
                                                title='Phenotypic Characterization of a Comprehensive Set of MAPK1/ERK2 Missense Mutants')

#### Ashenberg 2016 Flu Nucleoprotein ####
deep_mut_data$ashenberg_2016_np <- read_csv('data/raw/processed/ashenberg_2017_flu_np.csv')

df <- deep_mut_data$ashenberg_2016_np %>%
  mutate(variants = str_c('p.', wt, site, mut),
         raw_score = diffsel) %>%
  select(variants, score = diffsel, raw_score)

formatted_deep_data$ashenberg_2016_np <- DeepMut(variant_data = df, gene_name = 'Nucleoprotein', species = species_options$flu,
                                                 aa_seq = str_c(raw_seqs$`H3N2_A-Aichi-2-1968-nucleoprotein`, collapse = ''), transform = 'None',
                                                 authour = 'Ashenberg et al.', year = 2016,
                                                 substrain = 'Human adapted strain A/Aichi/2/1968, H3N2',
                                                 title='Deep mutational scanning identifies sites in influenza nucleoprotein that affect viral inhibition by MxA',
                                                 doi='10.1371/journal.ppat.1006288', pmid='28346537',
                                                 url='https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1006288')

#### Weile 2017 ube2i ####
# deep_mut_data$weile_2017_ube2i <- read_csv('data/raw/processed/weile_2017_ube2i_score_comp.csv', na = c('NA','','None'))
# 
# aa_code <- structure(names(Biostrings::AMINO_ACID_CODE), names=Biostrings::AMINO_ACID_CODE)
# df <- deep_mut_data$weile_2017_ube2i %>%
#   mutate(variants = str_replace_all(hgvs_pro, aa_code),
#          score = log2(abs(pred.score))) %>%
#   select(variants, score, raw_score = pred.score)
# 
# formatted_deep_data$weile_2017_ube2i <- DeepMut(variant_data = df, gene_name = 'UBE2I', species = species_options$sapiens,
#                                                 aa_seq = str_c(raw_seqs$h_sapiens_ube2i), transform = 'log2(abs(x))', uniprot_id = 'P63279',
#                                                 authour = 'Weile et al.', year = 2017,
#                                                 misc = list(title='A framework for exhaustively mapping functional missense variants',
#                                                             doi='10.15252/msb.20177908', pmid='29269382',
#                                                             url='http://msb.embopress.org/content/13/12/957'))

deep_mut_data$weile_2017_ube2i <- read_tsv('data/raw/processed/weile_2017_raw_counts_ube2i_tileseq.tsv') %>%
  mutate_at(.vars = vars(nonselect1, nonselect2, select1, select2), funs(pseudocount = . + min(.[. > 0]))) %>%
  mutate(mean_nonselect = (nonselect1_pseudocount + nonselect2_pseudocount)/2,
         mean_select = (select1_pseudocount + select2_pseudocount)/2,
         er = mean_select / mean_nonselect)

wt_er <- median(filter(deep_mut_data$weile_2017_ube2i, annotation == 'SYN') %>% pull(er))

deep_mut_data$weile_2017_ube2i %<>% mutate(fitness = er/wt_er,
                                           log2 = log2(fitness))

df <- deep_mut_data$weile_2017_ube2i %>%
  mutate(variants = str_c('p.', wt_aa, pos, mut_aa)) %>%
  select(variants, score = log2, raw_score = fitness, wt_codon, mut_codon, annotation)

formatted_deep_data$weile_2017_ube2i <- DeepMut(variant_data = df, gene_name = 'UBE2I', species = species_options$sapiens,
                                                aa_seq = str_c(raw_seqs$h_sapiens_ube2i, collapse = ''), transform = 'log2(abs(x))',
                                                uniprot_id = 'P63279', authour = 'Weile et al.', year = 2017,
                                                title='A framework for exhaustively mapping functional missense variants',
                                                doi='10.15252/msb.20177908', pmid='29269382', pdb_id = c('1A3S:A', '1KPS:A'),
                                                url='http://msb.embopress.org/content/13/12/957')

#### Weile 2017 sumo1 ####
# deep_mut_data$weile_2017_sumo1 <- read_csv('data/raw/processed/weile_2017_sumo1_score_comp.csv', na = c('NA','','None'))
# 
# df <- deep_mut_data$weile_2017_sumo1 %>%
#   mutate(variants = str_replace_all(hgvs_pro, aa_code)) %>%
#   select(variants, score, raw_score = pred.score)

deep_mut_data$weile_2017_sumo1 <- read_tsv('data/raw/processed/weile_2017_raw_counts_sumo1_tileseq.tsv') %>%
  mutate_at(.vars = vars(nonselect1, nonselect2, select1, select2), funs(pseudocount = . + min(.[. > 0]))) %>%
  mutate(mean_nonselect = (nonselect1_pseudocount + nonselect2_pseudocount)/2,
         mean_select = (select1_pseudocount + select2_pseudocount)/2,
         er = mean_select / mean_nonselect)

wt_er <- median(filter(deep_mut_data$weile_2017_sumo1, annotation == 'SYN') %>% pull(er))

deep_mut_data$weile_2017_sumo1 %<>% mutate(fitness = er/wt_er,
                                           log2 = log2(fitness))

df <- deep_mut_data$weile_2017_sumo1 %>%
  mutate(variants = str_c('p.', wt_aa, pos, mut_aa)) %>%
  select(variants, score = log2, raw_score = fitness, wt_codon, mut_codon, annotation)

formatted_deep_data$weile_2017_sumo1 <- DeepMut(variant_data = df, gene_name = 'SUMO1', species = species_options$sapiens, transform = 'log2',
                                                aa_seq = str_c(raw_seqs$h_sapiens_sumo1, collapse = ''), uniprot_id = 'P63165',
                                                authour = 'Weile et al.', year = 2017, pdb_id = c('1WYW:B', '2PE6:B'),
                                                title='A framework for exhaustively mapping functional missense variants',
                                                doi='10.15252/msb.20177908', pmid='29269382',
                                                url='http://msb.embopress.org/content/13/12/957')

#### Weile 2017 tpk1 ####
## Does not have raw count data on lab website
deep_mut_data$weile_2017_tpk1 <- read_csv('data/raw/processed/weile_2017_tpk1_score_comp.csv', na = c('NA','','None'))

#### Weile 2017 calm ####
## Does not have raw count data on lab website
deep_mut_data$weile_2017_calm <- read_csv('data/raw/processed/weile_2017_calm_score_comp.csv', na = c('NA','','None'))

#### Findlay 2018 Brca1 ####
deep_mut_data$findlay_2018_brca1 <- read_xlsx('data/raw/processed/findlay_2018_brca1_ring_brct.xlsx', skip = 2, na = 'NA') %>%
  rename_all(funs(gsub('[\\/ \\(\\)]+', '_', .))) %>%
  rename(ref_nuc = reference,
         alt_nuc = alt,
         ref_aa = aa_ref,
         alt_aa = aa_alt)

# Currently only interested in protein variants
df <- drop_na(deep_mut_data$findlay_2018_brca1, aa_pos) %>%
  mutate(score = function.score.mean,
         raw_score = score,
         variants = protein_variant) %>%
  select(variants, score, raw_score, p.nonfunctional, transcript_variant, position_hg19_)

formatted_deep_data$findlay_2018_brca1 <- DeepMut(variant_data = df, gene_name = 'BRCA1', species = species_options$sapiens,
                                                  aa_seq = str_c(raw_seqs$h_sapiens_brca1, collapse = ''), transform = 'None',
                                                  uniprot_id = 'P38398', authour = 'Findlay et al.', year = 2018,
                                                  title='Accurate classification of BRCA1 variants with saturation genome editing',
                                                  doi='10.1038/s41586-018-0461-z', pmid='',
                                                  url='https://www.nature.com/articles/s41586-018-0461-z')

#### Lee 2018 Flu haemagglutinin ####
# Xlsx also has sheets with preferences for each repeat, but importing just the average result for now
# Reports per site AA preference, which is fairly dissimilar to other metrics
deep_mut_data$lee_2018_flu_ha <- read_xlsx('data/raw/processed/lee_2018_influenza_ha.xlsx', sheet = 'avg_prefs') %>%
  rename(position = site) %>%
  gather(key = 'alt_aa', value = 'aa_pref', -position, -entropy, -neffective)

lee_position <- function(x){
  if (grepl('\\-', x)){
    return(as.numeric(x) + 17)
  } else if (grepl('HA2', x)) {
    x <- gsub('\\(HA2\\)', '', x)
    return(as.numeric(x) + 345)
  } else {
    return(as.numeric(x) + 16)
  }
}

# TODO use a more objective transformation to turn AA prefs into an approximation of fitness
df <- deep_mut_data$lee_2018_flu_ha %>%
  mutate(score = log2(aa_pref * 5),
         pos = sapply(position, lee_position),
         ref_aa = raw_seqs$`H3N2_A-Perth-16-2009-Hemagglutinin`[pos],
         variants = str_c('p.', ref_aa, pos, alt_aa)) %>%
  select(variants, score, raw_score=aa_pref, entropy, neffective)

formatted_deep_data$lee_2018_flu_ha <- DeepMut(variant_data = df, gene_name = 'HA', species = species_options$flu,
                                               aa_seq = str_c(raw_seqs$`H3N2_A-Perth-16-2009-Hemagglutinin`, collapse = ''), transform = 'log2(x*5)',
                                               authour = 'Lee et al.', year = 2018,
                                               notes='Used a simple placeholder transformation to change AA preference (raw_score) into something resembling other studies',
                                               substrain='Human adapted A/Perth/16/2009, H3N2', pmid='30104379', doi='10.1073/pnas.1806133115',
                                               url='www.pnas.org/cgi/doi/10.1073/pnas.1806133115',
                                               title='Deep mutational scanning of hemagglutinin helps predict evolutionary fates of human H3N2 influenza variants')

#### Giacomelli 2018 tp53 ####
deep_mut_data$giacomelli_2018_tp53 <- read_xlsx('data/raw/processed/giacomelli_2018_tp53.xlsx', skip=1) %>%
  rename_all(funs(tolower(gsub('[\\/ \\(\\)\\-]+', '_', gsub('\\*', '_2', .))))) %>%
  rename(ref_aa = aa_wt,
         alt_aa = aa_variant)

df1 <- deep_mut_data$giacomelli_2018_tp53 %>%
  mutate(score = a549_p53wt_nutlin_3_z_score,
         raw_score = score,
         variants = str_c('p.', allele)) %>%
  select(variants, score, raw_score, a549_p53wt_nutlin_3_z_score, a549_p53null_nutlin_3_z_score, a549_p53null_etoposide_z_score)

df2 <- deep_mut_data$giacomelli_2018_tp53 %>% 
  mutate(score = a549_p53null_nutlin_3_z_score,
         raw_score = score,
         variants = str_c('p.', allele)) %>%
  select(variants, score, raw_score, a549_p53wt_nutlin_3_z_score, a549_p53null_nutlin_3_z_score, a549_p53null_etoposide_z_score)

df3 <- deep_mut_data$giacomelli_2018_tp53 %>%
  mutate(score = a549_p53null_etoposide_z_score,
         raw_score = score,
         variants = str_c('p.', allele)) %>%
  select(variants, score, raw_score, a549_p53wt_nutlin_3_z_score, a549_p53null_nutlin_3_z_score, a549_p53null_etoposide_z_score)

tp53_ref_seq <- deep_mut_data$giacomelli_2018_tp53 %>% group_by(position) %>% summarise(aa = first(ref_aa)) %>% pull(aa) %>% str_c(.,collapse = '')

formatted_deep_data$giacomelli_2018_tp53 <- DeepMutSet(
  list(p53_wt_nutlin3=DeepMut(variant_data = df1, gene_name = 'TP53', species = species_options$sapiens,
                              authour = 'Giacomelli et al.', year = 2018, transform = 'None',
                              uniprot_id = 'P04637', aa_seq = tp53_ref_seq,
                              notes='Score is in Z-score format, with wt p53 background and nutlin3 selecting for dominant negative variants (see paper)',
                              title='Mutational processes shape the landscape of TP53 mutations in human cancer',
                              doi='10.1038/s41588-018-0204-y', pmid='30224644',
                              url='https://www.nature.com/articles/s41588-018-0204-y'),
       p53_null_nutlin3=DeepMut(variant_data = df2,
                                gene_name = 'TP53', species = species_options$sapiens,
                                authour = 'Giacomelli et al.', year = 2018, transform = 'None',
                                uniprot_id = 'P04637', aa_seq = tp53_ref_seq,
                                notes='Score is in Z-score format, with null p53 background and nutlin3 selecting for LOF variants (see paper)',
                                title='Mutational processes shape the landscape of TP53 mutations in human cancer',
                                doi='10.1038/s41588-018-0204-y', pmid='30224644',
                                url='https://www.nature.com/articles/s41588-018-0204-y'),
       p53_null_etoposide=DeepMut(variant_data = df3,
                                  gene_name = 'TP53', species = species_options$sapiens,
                                  authour = 'Giacomelli et al.', year = 2018, transform = 'None',
                                  uniprot_id = 'P04637', aa_seq = tp53_ref_seq,
                                  notes='Score is in Z-score format, with null p53 background and etoposide selecting for WT like (or better) variants (see paper)',
                                  title='Mutational processes shape the landscape of TP53 mutations in human cancer',
                                  doi='10.1038/s41588-018-0204-y', pmid='30224644',
                                  url='https://www.nature.com/articles/s41588-018-0204-y')
       ))

#### Save processed output ####
saveRDS(deep_mut_data, 'data/raw_deep_mut_data.RDS')
saveRDS(formatted_deep_data, 'data/formatted_deep_mut_data.RDS')

# Files written in form /UniprotID_GeneName.set.dm
for (i in names(formatted_deep_data)){
  print(str_c('Writing DeepMut: ', i))
  if ('DeepMutSet' %in% class(formatted_deep_data[[i]])){
    write_deep_mut(formatted_deep_data[[i]], str_c('data/standardised/', i))
  } else if ('DeepMut' %in% class(formatted_deep_data[[i]])){
    dm_path <- gsub(' ', '',
                    str_c('data/standardised/', i, '/', str_replace_na(formatted_deep_data[[i]]$uniprot_id, replacement = ''), '_',
                          str_replace_na(formatted_deep_data[[i]]$gene_name, replacement = ''), '.dm'))
    write_deep_mut(formatted_deep_data[[i]], dm_path)
  } else {
    stop('Incorrect object added to formatted_deep_data, not of class DeepMut or DeepMutSet')
  }
}

# Save meta data about proteins
get_meta <- function(x, var){
  if ('DeepMutSet' %in% class(x)){
    return(unique(sapply(x, get_meta, var=var)))
  } else if ('DeepMut' %in% class(x)){
    return(x[[var]])
  }
}
meta <- data_frame(gene_name = sapply(formatted_deep_data, get_meta, var='gene_name'),
                   uniprot_id = sapply(formatted_deep_data, get_meta, var='uniprot_id'))
write_tsv(meta, 'meta/gene_meta_data.tsv')

dataset_size <- sapply(deep_mut_data, function(x){dim(x)[1]}) %>% unlist()
