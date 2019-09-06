#!/usr/bin/env Rscript 
# Parse and process data ready for analysis of mutation properties in DMS data

source('src/config.R')
source('src/analysis/data_processing.R')
source('src/analysis/secondary_structure.R')

dms_data <- readRDS('data/rdata/processed_variant_data.RDS')

# Meta data about each study, including chosen threshold for deleterousness
meta_df <- tibble(study = names(dms_data),
                  gene_type = sapply(dms_data, function(x){get_meta(x$dm, 'gene_type')}),
                  gene_name = sapply(dms_data, function(x){get_meta(x$dm, 'gene_name')}),
                  uniprot_id = sapply(dms_data, function(x){get_meta(x$dm, 'uniprot_id')}),
                  test_class = sapply(dms_data, function(x){get_meta(x$dm, 'test_class')}),
                  species = sapply(dms_data, function(x){get_meta(x$dm, 'species')}),
                  authour = sapply(dms_data, function(x){get_meta(x$dm, 'authour')}),
                  group = sapply(dms_data, function(x){get_meta(x$dm, 'group')}),
                  single = sapply(dms_data, function(x){"single_variant" %in% class(x)}),
                  thresh=sapply(dms_data, function(x){x$manual_threshold}),
                  factor=sapply(dms_data, function(x){x$norm_factor}),
                  norm_thresh=thresh/factor)
write_tsv(meta_df, 'meta/study_meta_data.tsv')
saveRDS(meta_df, 'data/rdata/study_meta_data.RDS')

# Helix propensity
a_helix_propensity <- read_tsv('meta/alpha_helix_energy', comment = '#') %>%
  rename_all(str_to_lower) %>%
  mutate(aa1 = aa_3_to_1(aa)) %>%
  select(aa1, aa3 = aa, exptl:luque)
saveRDS(a_helix_propensity, 'data/rdata/a_helix_propensity.RDS')

# Surface accessibility scores for all variants
surface_accesibility <- sapply(dms_data, function(x){if(!identical(NA, x$surface_accesibility)){x$surface_accesibility$combined}},
                               simplify = FALSE) %>%
  bind_rows(.id = 'study')
saveRDS(surface_accesibility, 'data/rdata/position_surface_accesibility.RDS')

# Secondary structure predictions for all positions
secondary_structure <- sapply(dms_data, function(x){if(!identical(NA, x$secondary_structure)){x$secondary_structure}},
                              simplify = FALSE) %>%
  bind_rows(.id = 'study') %>%
  mutate(ss_reduced = SS_REDUCED_HASH[ss])
saveRDS(secondary_structure, 'data/rdata/position_secondary_structure.RDS')

# Hydrophobicity
hydrophobicity <- read_tsv('meta/residue_hydrophobicity.tsv', comment = '#',
                           col_types = cols(AA = col_character(), .default = col_double())) %>%
  rename_all(str_to_lower)
saveRDS(hydrophobicity, 'data/rdata/hydrophobicity.RDS')

# Parse FoldX data for all structures
foldx_preds <- sapply(dms_data, function(x){if (!identical(x$foldx, NA)) bind_rows(x$foldx, .id='pdb_id') else NULL}, simplify = FALSE) %>%
  bind_rows(.id = 'study') %>%
  filter(sapply(variants, str_count, pattern=',') == 0) %>%
  mutate(wt = str_sub(variants, end = 1),
         mut = str_sub(variants, start = -1),
         position = as.integer(str_sub(variants, start = 2, end = -2))) %>% 
  select(-variants) %>%
  select(study, pdb_id, position, wt, mut, everything())
saveRDS(foldx_preds, 'data/rdata/all_foldx.RDS')

# Backbone angles for all structures
backbone_angles <- sapply(dms_data, function(x){if(!identical(NA, x$backbone_angles)){x$backbone_angles}}, simplify = FALSE) %>%
  bind_rows(.id = 'study')
saveRDS(backbone_angles, 'data/rdata/backbone_angles.RDS')

# Chemical environment profiles
chemical_environments <- sapply(dms_data,
                                function(x){if(!identical(NA, x$chem_env)){if(!identical(NA, x$chem_env$combine_long)){x$chem_env$combine_long}}},
                                simplify = TRUE) %>%
  bind_rows(.id = 'study')  %>%
  rename(struct_group = group) %>%
  left_join(., meta_df, by = 'study') %>%
  left_join(., select(secondary_structure, study, position = pos, aa, ss, ss_reduced), by = c('study', 'position', 'aa')) %>%
  left_join(., rename(surface_accesibility, aa=res1, position=pos), by = c('study', 'position', 'aa')) %>%
  mutate(aa_reduced=AA_REDUCED_HASH[aa]) %>%
  group_by(pdb_id) %>%
  mutate(relative_position = position/max(position)) %>%
  ungroup()
saveRDS(chemical_environments, 'data/rdata/position_chemical_environments.RDS')

# Read PTMs
ptms <- readRDS('data/mutfunc/human/ptm/psp_ptms_parsed.rds') %>%
  as_tibble() %>%
  transmute(modification = str_to_lower(modification), uniprot_id = entry, position) %>%
  bind_rows(., read_tsv('meta/human_phosphosites_ochoa_et_al.tsv') %>%
              transmute(modification = 'phosphorylation', uniprot_id = uniprot, position))
saveRDS(ptms, 'data/rdata/human_ptms.RDS')

# Load SIFT scores
# sift <- read_tsv('data/mutfunc/human/conservation/sift_parsed_all.tab') %>%
#   rename(wt=ref, mut=alt, uniprot_id = acc, position = pos) %>%
#   filter(wt %in% Biostrings::AA_STANDARD) %>%
#   spread(key = mut, value = score) %>%
#   left_join(ptms, by = c('uniprot_id', 'position')) %>%
#   mutate(modification = ifelse(is.na(modification), 'none', modification))
# saveRDS(sift, 'data/rdata/human_sift.RDS')
sift <- readRDS('data/rdata/human_sift.RDS') # Only store cached file locally normally to save disk space

sift_reduced <- filter(sift, uniprot_id %in% sample(sift$uniprot_id, 100))
saveRDS(sift_reduced, 'data/rdata/human_sift_reduced.RDS')

sift_ptms <- bind_rows(sample_n(filter(sift, modification == 'none'), 50000),
                       sample_n(filter(sift, !modification == 'none'), 50000))
saveRDS(sift_ptms, 'data/rdata/human_sift_ptms.RDS')

# Load FoldX Scores - uncomment to regenerate
# foldx <- read_tsv('data/mutfunc/human/structure/exp_full.tab') %>%
#   rename_all(.funs = ~str_replace_all(., ' ', '_')) %>%
#   filter(pos > 0) %>%
#   arrange(uniprot_id, pdb_id, pos) %>%
#   distinct(uniprot_id, pos, mut, .keep_all = TRUE) %>% # Choose one pdb_id per position (randomly TODO - change to better filter)
#   rename(position = pos) %>%
#   left_join(., ptms, by = c('uniprot_id', 'position')) %>%
#   mutate(modification = ifelse(is.na(modification), 'none', modification), m = TRUE) %>%
#   distinct() %>%
#   spread(key = modification, value = m, fill=FALSE) %>%
#   rename(modified = none,
#          o_glcnac_glycosilation = `o-linked glycosilation (o-glcnac)`,
#          o_galnac_glycosilation = `o-linked glycosilation (o-galnac)`) %>%
#   select(uniprot_id:mut, modified, sd:entropy_complex, acetylation:ubiquitination) %>%
#   mutate(modified = !modified)
# saveRDS(foldx, 'data/rdata/human_foldx.RDS')
foldx <- readRDS('data/rdata/human_foldx.RDS')

foldx_tiny <- filter(foldx, uniprot_id %in% sample(foldx$uniprot_id, 10))
saveRDS(foldx_tiny, 'data/rdata/human_foldx_tiny.RDS')

foldx_reduced <- filter(foldx, uniprot_id %in% sample(foldx$uniprot_id, 50))
saveRDS(foldx_reduced, 'data/rdata/human_foldx_reduced.RDS')

# Sample 1000 positions with/without ptms
foldx_ptms <- distinct(foldx, uniprot_id, position, modified) %>%
  group_by(modified) %>%
  sample_n(1000) %>%
  ungroup() %>%
  semi_join(foldx, ., by = c('uniprot_id', 'position', 'modified'))
saveRDS(foldx_ptms, 'data/rdata/human_foldx_ptms.RDS')

# Dataframe of all individually scored variant/sets of variants in all studies
all_variants <- bind_rows(lapply(dms_data, function(x){x$dm$variant_data}), .id = 'study') %>%
  select(study, variants, score, raw_score, norm_score) %>%
  mutate(variants = str_replace_all(variants, 'p.', '')) %>%
  left_join(., meta_df, by='study')
saveRDS(all_variants, 'data/rdata/all_study_variants.RDS')

# Dataframe with matrix of mutational profiles for all positions and all studies
variant_matrices <- list()
variant_matrices$all_variants <- bind_rows(lapply(dms_data, make_var_matrix, score='score'), .id = 'study') %>%
  filter(!wt %in% c('Z', 'B')) %>%
  drop_na(wt, pos) %>%
  left_join(., meta_df, by='study') %>%
  left_join(., rename(surface_accesibility, wt=res1), by=c('study', 'wt', 'pos')) %>%
  left_join(., rename(secondary_structure, wt=aa) %>%
               rename_at(vars(-study, -pos, -wt, -ss, -ss_reduced), .funs= ~ str_c('ss_prob_', .)),
            by=c('study', 'wt', 'pos')) %>%
  left_join(., drop_na(backbone_angles, phi, psi) %>%
              group_by(study, aa, position) %>% 
              summarise(phi = first(phi), psi = first(psi)) %>%
              rename(pos = position, wt = aa),
            by = c("study", "pos", "wt")) %>%
  select(-factor, -norm_thresh) %>%
  mutate(sig_count = mutate_at(., .vars = vars(A:Y), .funs = list(~ . < thresh)) %>%
           select(A:Y) %>%
           rowSums(na.rm = TRUE))

variant_matrices$norm_all_variants <- bind_rows(lapply(dms_data, make_var_matrix, score='norm_score'), .id = 'study') %>%
  filter(!wt %in% c('Z', 'B')) %>%
  drop_na(wt, pos) %>%
  left_join(., meta_df, by='study') %>%
  left_join(., rename(surface_accesibility, wt=res1), by=c('study', 'wt', 'pos')) %>%
  left_join(., rename(secondary_structure, wt=aa) %>%
               rename_at(vars(-study, -pos, -wt, -ss, -ss_reduced), .funs= ~ str_c('ss_prob_', .)),
            by=c('study', 'wt', 'pos')) %>%
  left_join(., drop_na(backbone_angles, phi, psi) %>%
              group_by(study, aa, position) %>% 
              summarise(phi = first(phi), psi = first(psi)) %>%
              rename(pos = position, wt = aa),
            by = c("study", "pos", "wt")) %>%
  select(-factor, -thresh, thresh = norm_thresh) %>%
  mutate(sig_count = mutate_at(., .vars = vars(A:Y), .funs = list(~ . < thresh)) %>%
           select(A:Y) %>%
           rowSums(na.rm = TRUE))

# Label secondary structure before filtering data
variant_matrices <- sapply(variant_matrices, label_secondary_structure, simplify = FALSE)

# Filter interesting subsets
variant_matrices$sig_positions <- filter(variant_matrices$all_variants, sig_count > 0)
variant_matrices$norm_sig_positions <- filter(variant_matrices$norm_all_variants, sig_count > 0)
variant_matrices$single_variants <- filter(variant_matrices$all_variants, single)
variant_matrices$sig_single_variants <- filter(variant_matrices$sig_positions, single)
variant_matrices$norm_single_variants <- filter(variant_matrices$norm_all_variants, single)
variant_matrices$norm_sig_single_variants <- filter(variant_matrices$norm_sig_positions, single)

imputed_matrices <- mapply(impute_variant_profiles, variant_matrices,
                           list(variant_matrices$all_variants, variant_matrices$norm_all_variants, variant_matrices$all_variants,
                                variant_matrices$norm_all_variants, variant_matrices$all_variants, variant_matrices$all_variants,
                                variant_matrices$norm_all_variants, variant_matrices$norm_all_variants), SIMPLIFY = FALSE)

saveRDS(variant_matrices, 'data/rdata/all_study_position_matrices.RDS')
saveRDS(imputed_matrices, 'data/rdata/all_study_imputed_position_matrices.RDS')
