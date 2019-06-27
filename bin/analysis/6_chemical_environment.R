#!/usr/bin/env Rscript 
# Script analysing chemical environment in proteins, in the context of our deep mutagenesis data

source('src/config.R')
source('src/analysis/chemical_environment.R')
source('src/analysis/position_profile_clustering.R')

variant_matrices <- readRDS('data/rdata/all_study_position_matrices.RDS')
chemical_environments <- readRDS('data/rdata/position_chemical_environments.RDS') %>%
  left_join(., select(variant_matrices$norm_all_variants, study, position=pos, aa=wt, sig_count, A:Y), by = c("study", "position", "aa"))
chem_env_deduped <- distinct(chemical_environments, pdb_id, chain, position, aa, .keep_all = TRUE)

chem_env_metric_names <- c('nearest_10', 'within_10.0')

plots <- list()

#### Basic profile analysis ####
plots <- list_modify(plots, !!!sapply(chem_env_metric_names, function(x){
  plot_basic_profile_analysis(chem_env_deduped, !!x, names = str_c('prof_', sorted_aa_1_code)) %>%
    set_names(str_c('aa_profile_', names(.)))
}, simplify = FALSE))

## Linear model from profiles (currently just with within_10.0 profile, should rework whole script to better cope with multiples)
chem_env_raw_prof_lm <- calc_all_profile_lms(chemical_environments, within_10.0, A:Y, prof_col_names = str_c('prof_', sorted_aa_1_code))

plots$within_10.0$lm_summary <- labeled_ggplot(
  ggplot(chem_env_raw_prof_lm, aes(x=mut_aa, y=r.squared, fill=-log10(p.value), label=n)) + 
    facet_wrap(~study) + 
    geom_col() + 
    geom_text(colour='red', check_overlap = FALSE, angle=90, hjust=0, vjust=0.5, nudge_y = 0.05),
  width=15, height=10)

tmp <- filter(chem_env_raw_prof_lm, study=='ALL', mut_aa=='P') %>%
  pull(model)
tmp <- augment(tmp[[1]])
plots$within_10.0$lm_example1 <- ggplot(tmp, aes(x=er, y=.fitted, colour=.se.fit)) + geom_point() + ggtitle('Proline, ALL')

tmp <- filter(chem_env_raw_prof_lm, study=='weile_2017_ube2i', mut_aa=='K') %>%
  pull(model)
tmp <- augment(tmp[[1]])
plots$within_10.0$lm_example2 <- ggplot(tmp, aes(x=er, y=.fitted, colour=.se.fit)) + geom_point() + ggtitle('Lysine, Weile 2017 ube2i')

tmp <- filter(chem_env_raw_prof_lm, study=='brenan_2016_erk2', mut_aa=='H') %>%
  pull(model)
tmp <- augment(tmp[[1]])
plots$within_10.0$lm_example3 <- ggplot(tmp, aes(x=er, y=.fitted, colour=.se.fit)) + geom_point() + ggtitle('Histadine, Brenan 2016 Erk2')

# Linear model with significance of position included
chem_env_prof_sig_count_lm <- calc_all_profile_lms(chemical_environments, within_10.0, A:Y,
                                                   prof_col_names = str_c('prof_', sorted_aa_1_code), include_sig_count = TRUE)

plots$within_10.0$lm_sig_count_summary <- labeled_ggplot(
  ggplot(chem_env_prof_sig_count_lm, aes(x=mut_aa, y=r.squared, fill=-log10(p.value), label=n)) + 
    facet_wrap(~study) + 
    geom_col() + 
    geom_text(colour='red', check_overlap = FALSE, angle=90, hjust=0, vjust=0.5, nudge_y = 0.05),
  width=15, height=10)

tmp <- filter(chem_env_prof_sig_count_lm, study=='ALL', mut_aa=='P') %>%
  pull(model)
tmp <- augment(tmp[[1]])
plots$within_10.0$lm_sig_count_example1 <- ggplot(tmp, aes(x=er, y=.fitted, colour=.se.fit)) + geom_point() + ggtitle('Proline, ALL')

tmp <- filter(chem_env_prof_sig_count_lm, study=='weile_2017_ube2i', mut_aa=='K') %>%
  pull(model)
tmp <- augment(tmp[[1]])
plots$within_10.0$lm_sig_count_example2 <- ggplot(tmp, aes(x=er, y=.fitted, colour=.se.fit)) + geom_point() + ggtitle('Lysine, Weile 2017 ube2i')

tmp <- filter(chem_env_prof_sig_count_lm, study=='brenan_2016_erk2', mut_aa=='H') %>%
  pull(model)
tmp <- augment(tmp[[1]])
plots$within_10.0$lm_sig_count_example3 <- ggplot(tmp, aes(x=er, y=.fitted, colour=.se.fit)) + geom_point() + ggtitle('Histadine, Brenan 2016 Erk2')

chem_env_prof_sig_count_lm_loadings <- filter(chem_env_prof_sig_count_lm, study == 'ALL') %>%
  select(mut_aa, model, n, r.squared) %>%
  mutate(coef_df = lapply(model, tidy)) %>%
  unnest(coef_df)

plots$within_10.0$lm_sig_count_loadings <- ggplot(chem_env_prof_sig_count_lm_loadings,
                                                  aes(x=mut_aa, y=term, fill=estimate)) + 
  geom_tile() +
  scale_fill_gradient2() +
  theme(axis.ticks = element_blank(), panel.background = element_blank()) +
  xlab('Substituted AA') +
  ylab('LM Term') +
  geom_point(data = filter(chem_env_prof_sig_count_lm_loadings, p.value < 0.0001), aes(shape='p < 0.0001')) + 
  geom_point(data = filter(chem_env_prof_sig_count_lm_loadings, p.value < 0.001, p.value > 0.0001), aes(shape='p < 0.001')) + 
  geom_point(data = filter(chem_env_prof_sig_count_lm_loadings, p.value < 0.01, p.value > 0.001), aes(shape='p < 0.01')) +
  scale_shape_manual(values = c('p < 0.0001'=8, 'p < 0.001'=3, 'p < 0.01'=20))

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

chem_env_reduced_full <- mutate_at(chemical_environments, .vars = vars(one_of(chem_env_metric_names)),
                                   .funs = ~sapply(., reduce_aa_profile, simplify = FALSE))

# Basic profiles
plots <- list_modify(plots, !!!sapply(chem_env_metric_names, function(x){
  plot_basic_profile_analysis(chem_env_reduced, !!x) %>%
    set_names(str_c('aa_groups_profile_', names(.)))
}, simplify = FALSE))

ch <- expand_profile_column(chem_env_reduced, within_10.0)

# LM of ER vs profile
chem_env_reduced_prof_lm <- calc_all_profile_lms(chem_env_reduced_full, within_10.0, A:Y, prof_col_names = NULL)
plots$within_10.0$lm_reduced_summary <- labeled_ggplot(
  ggplot(chem_env_reduced_prof_lm, aes(x=mut_aa, y=r.squared, fill=-log10(p.value), label=n)) + 
    facet_wrap(~study) + 
    geom_col() + 
    geom_text(colour='red', check_overlap = FALSE, angle=90, hjust=0, vjust=0.5, nudge_y = 0.05),
  width=15, height=10)

tmp <- filter(chem_env_reduced_prof_lm, study=='ALL', mut_aa=='P') %>%
  pull(model)
tmp <- augment(tmp[[1]])
plots$within_10.0$lm_reduced_example1 <- ggplot(tmp, aes(x=er, y=.fitted, colour=.se.fit)) + geom_point() + ggtitle('Proline, ALL')

tmp <- filter(chem_env_reduced_prof_lm, study=='weile_2017_ube2i', mut_aa=='K') %>%
  pull(model)
tmp <- augment(tmp[[1]])
plots$within_10.0$lm_reduced_example2 <- ggplot(tmp, aes(x=er, y=.fitted, colour=.se.fit)) + geom_point() + ggtitle('Lysine, Weile 2017 ube2i')

tmp <- filter(chem_env_reduced_prof_lm, study=='brenan_2016_erk2', mut_aa=='H') %>%
  pull(model)
tmp <- augment(tmp[[1]])
plots$within_10.0$lm_reduced_example3 <- ggplot(tmp, aes(x=er, y=.fitted, colour=.se.fit)) + geom_point() + ggtitle('Histadine, Brenan 2016 Erk2')


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
