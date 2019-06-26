#!/usr/bin/env Rscript 
# Script analysing chemical environment in proteins, in the context of our deep mutagenesis data

source('src/config.R')
source('src/analysis/chemical_environment.R')
source('src/analysis/position_profile_clustering.R')

variant_matrices <- readRDS('data/rdata/all_study_position_matrices.RDS')
chemical_environments <- readRDS('data/rdata/position_chemical_environments.RDS') %>%
  left_join(., select(variant_matrices$all_variants, study, position=pos, aa=wt, sig_count), by = c("study", "position", "aa"))
chem_env_deduped <- distinct(chemical_environments, pdb_id, chain, position, aa, .keep_all = TRUE)

chem_env_metric_names <- c('nearest_10', 'within_10.0')

plots <- list()

#### Basic profile analysis ####
plots <- list_modify(plots, !!!sapply(chem_env_metric_names, function(x){
  plot_basic_profile_analysis(chem_env_deduped, !!x, names = sorted_aa_1_code) %>%
    set_names(str_c('aa_profile_', names(.)))
}, simplify = FALSE))

## PCA on profiles
chem_env_pcas <- sapply(chem_env_metric_names, function(x){chem_env_pca(chem_env_deduped, var=x, names=sort(Biostrings::AA_STANDARD))},
                        simplify = FALSE)

plots <- list_modify(plots, !!!sapply(chem_env_pcas, function(x){list(
  pca=plot_chem_env_basic_pca_plots(x, cont_factors=c('all_atom_rel', 'relative_position', 'sig_count'),
                                    discrete_factors=c('ss', 'ss_reduced', 'aa', 'aa_reduced', 'pdb_id', 'gene_name', 'species'))
  )}, simplify = FALSE))

## Quantify relationship of factors with PCs
chem_env_pca_factor_cors <- sapply(chem_env_pcas, pca_factor_cor, .vars = vars(all_atom_abs:polar_rel, relative_position, sig_count),
                                   simplify = FALSE)
plots <- list_modify(plots, !!!sapply(chem_env_pca_factor_cors, function(x){list(pca=list(pca_factor_cor_heatmap=pca_factor_heatmap(x)))},
                                      simplify = FALSE))

# # Test factors against PCs
# p <- ggplot(gather(chem_env_pcas$within_10.0$profiles, key='pc', value = 'value', starts_with('PC')), aes(x=ss_reduced, y=value)) +
#   facet_wrap(~pc) +
#   geom_boxplot()

## Profile tSNE
chem_env_tsnes <- sapply(chem_env_metric_names, function(x){chem_env_tsne(chem_env_deduped, var=x, names=sort(Biostrings::AA_STANDARD))},
                         simplify = FALSE)
plots <- list_modify(plots, !!!sapply(chem_env_tsnes, function(x){
  list(tSNE=plot_factors(x$profiles, tSNE1, tSNE2, quos(ss_reduced=ss_reduced, aa_reduced=aa_reduced, gene_name=gene_name,
                                                        sig_count=sig_count, sqrt_suf_acc=sqrt(all_atom_rel))))
}, simplify = FALSE))

# tSNE plot tests
p <- ggplot(chem_env_tsnes$within_10.0$profiles, aes(x=tSNE1, y=tSNE2, colour=gene_name)) + geom_point()

########

#### Profiles reduced to AA groups ####
chem_env_reduced <- mutate_at(chem_env_deduped, .vars = vars(one_of(chem_env_metric_names)),
                              .funs = ~sapply(., reduce_aa_profile, simplify = FALSE))

# Basic profiles
plots <- list_modify(plots, !!!sapply(chem_env_metric_names, function(x){
  plot_basic_profile_analysis(chem_env_reduced, !!x) %>%
    set_names(str_c('aa_groups_profile_', names(.)))
}, simplify = FALSE))

ch <- expand_profile_column(chem_env_reduced, within_10.0)

# PCA
chem_env_reduced_pcas <- sapply(chem_env_metric_names, function(x){chem_env_pca(chem_env_reduced, var=x)}, simplify = FALSE)

plots <- list_modify(plots, !!!sapply(chem_env_reduced_pcas, function(x){list(
  pca=plot_chem_env_basic_pca_plots(x, cont_factors=c('all_atom_rel', 'relative_position', 'sig_count'),
                                    discrete_factors=c('ss', 'ss_reduced', 'aa', 'aa_reduced', 'pdb_id',
                                                       'gene_name', 'species'),
                                    max_pc=6, nrow=1, ncol=3) %>% # Are actually 7 classes, but doesn't plot well with current function
    set_names(str_c(names(.), '_grouped_aa_prof')))}, simplify = FALSE))

# Quantify relationship of factors with PCs
chem_env_reduced_pca_factor_cors <- sapply(chem_env_reduced_pcas, pca_factor_cor, .vars = vars(all_atom_abs:polar_rel, relative_position, sig_count),
                                           simplify = FALSE)
plots <- list_modify(plots, !!!sapply(chem_env_reduced_pca_factor_cors,
                                      function(x){list(pca=list(pca_factor_cor_heatmap_grouped_aa_prof=pca_factor_heatmap(x)))},
                                      simplify = FALSE))

# Profile tSNE
chem_env_tsnes <- sapply(chem_env_metric_names, function(x){chem_env_tsne(chem_env_reduced, var=x)}, simplify = FALSE)
plots <- list_modify(plots, !!!sapply(chem_env_tsnes, function(x){
  list(tSNE=plot_factors(x$profiles, tSNE1, tSNE2, quos(ss_reduced=ss_reduced, aa_reduced=aa_reduced, gene_name=gene_name,
                                                        sig_count=sig_count, sqrt_suf_acc=sqrt(all_atom_rel))) %>%
         set_names(str_c(names(.), '_grouped_aa_prof')))
}, simplify = FALSE))

########

# Save plots
save_plot_list(plots, root='figures/6_chemical_environment/')
