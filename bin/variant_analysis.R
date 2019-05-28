#!/usr/bin/env Rscript 
# Script to analyse deep mutagenesis data itself

source('bin/config.R')

#### Import and process data ####
dms_data <- readRDS('data/variant_data.RDS')

deep_variant_plots <- list(score_distributions=list())

# Meta data about each study, including chosen threshold for deleterousness
meta_df <- tibble(study = names(dms_data),
                  gene_type = sapply(dms_data, function(x){get_meta(x$dm, 'gene_type')}),
                  gene_name = sapply(dms_data, function(x){get_meta(x$dm, 'gene_name')}),
                  test_class = sapply(dms_data, function(x){get_meta(x$dm, 'test_class')}),
                  species = sapply(dms_data, function(x){get_meta(x$dm, 'species')}),
                  authour = sapply(dms_data, function(x){get_meta(x$dm, 'authour')}),
                  thresh=sapply(dms_data, function(x){x$manual_threshold}),
                  factor=sapply(dms_data, function(x){x$norm_factor}),
                  norm_thresh=thresh/factor)

# Dataframe of all individually scored variant/sets of variants in all studies
all_variants <- bind_rows(lapply(dms_data, function(x){x$dm$variant_data}), .id='study') %>%
  select(study, variants, score, raw_score, norm_score) %>%
  mutate(variants = str_replace_all(variants, 'p.', '')) %>%
  left_join(., meta_df, by='study')

# Dataframe with matrix of mutational profiles for all positions and all studies
variant_matrix <- bind_rows(lapply(dms_data, make_var_matrix, score='score'), .id = 'study') %>%
  separate(study, into = c('authour', 'year', 'gene'), sep='_', extra = 'merge', remove = FALSE) %>%
  select(-`*`, -`_`, -B, -Z, -`<NA>`) %>%
  drop_na(wt, pos)

# Per study/per AA median mutational profiles
per_study_mean_profiles <- group_by(variant_matrix, study, authour, year, gene, wt) %>%
  summarise_at(.vars = vars(-pos), .funs = median, na.rm=TRUE) %>%
  replace(is.na(.), 0)

# Matrix of mutational profiles with missing values imputed from the per study/AA medians
imputed_matrix <- gather(variant_matrix, key='mut', value='score', A:Y) %>%
  left_join(., gather(per_study_mean_profiles, key='mut', value='imp_score', A:Y),
            by = c("study", "authour", "year", "gene", "wt", "mut"))

imputed_matrix$score[is.na(imputed_matrix$score)] <- imputed_matrix$imp_score[is.na(imputed_matrix$score)]

imputed_matrix <- imputed_matrix %>%
  select(-imp_score) %>%
  spread(key=mut, value=score)

#### Enrichment Score Distributions ####
deep_variant_plots$score_distributions$dm_hists <- labeled_ggplot(p = plot_study_histogram(all_variants, meta_df),
                                                                  width = 14, height = 9)

deep_variant_plots$score_distributions$dm_norm_hists <- labeled_ggplot(p = plot_study_histogram(all_variants, meta_df,
                                                                                                x='norm_score', thresh = 'norm_thresh'),
                                                                       width = 14, height = 9)

deep_variant_plots$score_distributions$dm_gene_type_density <- labeled_ggplot(p = plot_factor_density(all_variants, facet = '~gene_type',
                                                                                                      x='norm_score'),
                                                                              width=14, height=9)

deep_variant_plots$score_distributions$dm_test_class_density <- labeled_ggplot(p = plot_factor_density(all_variants, facet = '~test_class',
                                                                                                       x='norm_score'),
                                                                               width=14, height=9)

deep_variant_plots$score_distributions$dm_species_density <- labeled_ggplot(p = plot_factor_density(all_variants, facet = '~species',
                                                                                                    x='norm_score'),
                                                                            width=14, height=9)

deep_variant_plots$score_distributions$dm_gene_density <- labeled_ggplot(p = plot_factor_density(all_variants, facet = '~gene_name',
                                                                                                 x='norm_score'),
                                                                         width=14, height=9)



per_position_boxplots <- gather(variant_matrix, key='mut', value='score', A:Y) %>%
  select(study, pos, wt, mut, score) %>%
  group_by(study) %>%
  do(plots = plot_score_distribution_per_position(.)) 

deep_variant_plots$score_distributions$per_position_boxplots <- per_position_boxplots$plots
names(deep_variant_plots$score_distributions$per_position_boxplots) <- per_position_boxplots$study
  
#### Variant Score PCAs ####
pca <- prcomp(as.matrix(select(imputed_matrix, -authour, -year, -gene, -wt, -pos)), center = FALSE, scale. = FALSE)

pca_variants <- bind_cols(select(imputed_matrix, authour, year, gene, wt, pos), as_tibble(pca$x)) %>%
  unite(col = 'study', authour, year, gene, remove = FALSE)

deep_variant_plots$pcas <- list()

# Plot PCAs againt basic possible underlying variables
deep_variant_plots$pcas$by_aa <- ggplot(pca_variants, aes(x=PC1, y=PC2, colour=study)) + geom_point()

deep_variant_plots$pcas$by_authour <- ggplot(pca_variants, aes(x=PC1, y=PC2, colour=study)) + geom_point()

deep_variant_plots$pcas$by_authour_fields <- ggplot(filter(pca_variants,
                                                           authour %in% c('araya', 'melamed', 'starita', 'kitzman', 'weile')),
                                                    aes(x=PC1, y=PC2, colour=study)) + geom_point()

#### Per AA PCAs ####

per_aa_pca <- function()


#### Save plots #####
save_plot_list(deep_variant_plots, root='figures/variant_analysis/')
