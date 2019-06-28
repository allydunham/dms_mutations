#!/usr/bin/env Rscript 
# Script analysing chemical environment in proteins, in the context of our deep mutagenesis data

source('src/config.R')
source('src/analysis/chemical_environment.R')
source('src/analysis/position_profile_clustering.R')

variant_matrices <- readRDS('data/rdata/all_study_position_matrices.RDS')

chemical_environments <- readRDS('data/rdata/position_chemical_environments.RDS') %>%
  left_join(., select(variant_matrices$norm_all_variants, study, position=pos, aa=wt, sig_count, A:Y), by = c("study", "position", "aa")) %>%
  mutate_at(.vars = vars(one_of(c('nearest_10', 'within_10.0'))),
            .funs = list(reduced=~sapply(., reduce_aa_profile, simplify = FALSE)))

chem_env_profiles <- c('nearest_10', 'nearest_10_reduced', 'within_10.0', 'within_10.0_reduced')
chem_env_profile_col_names <- list(str_c('prof_', SORTED_AA_1_CODE),
                                   str_c('prof_', names(AA_REDUCED_CLASSES)),
                                   str_c('prof_', SORTED_AA_1_CODE),
                                   str_c('prof_', names(AA_REDUCED_CLASSES)))

chem_env_anlyses <- mapply(function(prof_col, prof_names){analyse_chem_env_profile(chemical_environments, prof_col = !!prof_col,
                                                                                   prof_col_names = prof_names)},
                           SIMPLIFY =  FALSE, chem_env_profiles, chem_env_profile_col_names)
saveRDS(chem_env_anlyses, 'data/rdata/chemical_environment_analyses.RDS')

# Uncomment and comment above to switch to cached version
#chem_env_anlyses <- readRDS('data/rdata/chemical_environment_analyses.RDS')

plots <- sapply(chem_env_anlyses, function(x){x$plots}, simplify = FALSE)

# Save plots
save_plot_list(plots, root='figures/6_chemical_environment/')
