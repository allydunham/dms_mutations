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

imputed_matrices <- sapply(variant_matrices, impute_variant_profiles, background_matrix=variant_matrices$all_variants, simplify=FALSE)

saveRDS(variant_matrices, 'data/rdata/all_study_position_matrices.RDS')
saveRDS(imputed_matrices, 'data/rdata/all_study_imputed_position_matrices.RDS')
