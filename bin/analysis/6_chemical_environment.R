#!/usr/bin/env Rscript 
# Script analysing chemical environment in proteins, in the context of our deep mutagenesis data

source('src/config.R')
source('src/analysis/chemical_environment.R')
source('src/analysis/position_profile_clustering.R')

variant_matrices <- readRDS('data/rdata/all_study_position_matrices.RDS')
chemical_environments <- readRDS('data/rdata/position_chemical_environments.RDS') %>%
  left_join(., select(variant_matrices$all_variants, study, position=pos, aa=wt, sig_count), by = c("study", "position", "aa"))

chem_env_metric_names <- c('nearest_10', 'within_10.0')

plots <- list()

## PCA on profiles
chem_env_deduped <- distinct(chemical_environments, pdb_id, chain, position, aa, .keep_all = TRUE)
chem_env_pcas <- sapply(chem_env_metric_names, function(x){chem_env_pca(chem_env_deduped, var=x)}, simplify = FALSE)

plots <- list_modify(plots, !!!sapply(chem_env_pcas, function(x){list(
  pca=plot_chem_env_basic_pca_plots(x, cont_factors=c('all_atom_rel', 'relative_position', 'sig_count'),
                                    discrete_factors=c('ss', 'ss_reduced', 'aa', 'aa_reduced', 'pdb_id', 'gene_name', 'species'))
  )}, simplify = FALSE))

## Quantify relationship of factors with PCs
chem_env_pca_factor_cors <- sapply(chem_env_pcas, pca_factor_cor, .vars = vars(all_atom_abs:polar_rel, relative_position, sig_count),
                                   simplify = FALSE)
plots <- list_modify(plots, !!!sapply(chem_env_pca_factor_cors, function(x){list(pca=list(pca_factor_cor_heatmap=pca_factor_heatmap(x)))},
                                      simplify = FALSE))

# Test factors against PCs
p <- ggplot(gather(chem_env_pcas$within_10.0$profiles, key='pc', value = 'value', starts_with('PC')), aes(x=ss_reduced, y=value)) +
  facet_wrap(~pc) +
  geom_boxplot()

# Save plots
save_plot_list(plots, root='figures/6_chemical_environment/')
