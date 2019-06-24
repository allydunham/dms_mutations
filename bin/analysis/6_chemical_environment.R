#!/usr/bin/env Rscript 
# Script analysing chemical environment in proteins, in the context of our deep mutagenesis data

source('src/config.R')
source('src/analysis/chemical_environment.R')
source('src/analysis/position_profile_clustering.R')

chemical_environments <- readRDS('data/rdata/position_chemical_environments.RDS')

chem_env_metric_names <- c('nearest_10', 'within_10.0')
plots <- sapply(chem_env_metric_names, function(x){return(list())}, simplify = FALSE)

# Cluster profiles
chem_env_deduped <- distinct(chemical_environments, pdb_id, chain, position, aa, .keep_all = TRUE)
chem_env_pcas <- sapply(chem_env_metric_names, function(x){chem_env_pca(chem_env_deduped, var=x)}, simplify = FALSE)

factors <- c('ss', 'ss_reduced', 'all_atom_rel')
for (name in chem_env_metric_names){
  plots[[name]] <- list_modify(plots[[name]], pca=sapply(factors, function(x){plot_all_pcs(chem_env_pcas[[name]]$profiles, colour_var = x)},
                                                         simplify = FALSE))
}

# Save plots
save_plot_list(plots, root='figures/6_chemical_environment/')
