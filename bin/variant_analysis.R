#!/usr/bin/env Rscript 
# Script to analyse deep mutagenesis data itself

source('bin/config.R')

# Import data
dms_data <- readRDS('data/variant_data.RDS')

deep_variant_plots <- list(score_distributions=list())

#### Score Distribution Plots ####
meta_df <- data_frame(study = names(deep_variant_data),
                      gene_type = sapply(deep_variant_data, function(x){get_meta(x$dm, 'gene_type')}),
                      gene_name = sapply(deep_variant_data, function(x){get_meta(x$dm, 'gene_name')}),
                      test_class = sapply(deep_variant_data, function(x){get_meta(x$dm, 'test_class')}),
                      species = sapply(deep_variant_data, function(x){get_meta(x$dm, 'species')}),
                      authour = sapply(deep_variant_data, function(x){get_meta(x$dm, 'authour')}))

all_dm <- bind_rows(lapply(deep_variant_data, function(x){x$dm$variant_data}), .id='study') %>%
  select(study, variants, score, raw_score, norm_score) %>%
  mutate(variants = str_replace_all(variants, 'p.', '')) %>%
  left_join(., meta_df, by='study')

thresh_df <- data_frame(study=names(deep_variant_data),
                        thresh=sapply(deep_variant_data, function(x){x$manual_threshold}),
                        factor=sapply(deep_variant_data, function(x){x$norm_factor}),
                        norm_thresh=thresh/factor) %>%
  left_join(., meta_df, by='study')

deep_variant_plots$score_distributions$dm_hists <- labeled_ggplot(p = plot_study_histogram(all_dm, thresh_df),
                                                                  width = 14, height = 9)

deep_variant_plots$score_distributions$dm_norm_hists <- labeled_ggplot(p = plot_study_histogram(all_dm, thresh_df,
                                                                                                x='norm_score', thresh = 'norm_thresh'),
                                                                       width = 14, height = 9)

deep_variant_plots$score_distributions$dm_gene_type_density <- labeled_ggplot(p = plot_factor_density(all_dm, facet = '~gene_type',
                                                                                                      x='norm_score'),
                                                                              width=14, height=9)

deep_variant_plots$score_distributions$dm_test_class_density <- labeled_ggplot(p = plot_factor_density(all_dm, facet = '~test_class',
                                                                                                       x='norm_score'),
                                                                               width=14, height=9)

deep_variant_plots$score_distributions$dm_species_density <- labeled_ggplot(p = plot_factor_density(all_dm, facet = '~species',
                                                                                                    x='norm_score'),
                                                                            width=14, height=9)

deep_variant_plots$score_distributions$dm_gene_density <- labeled_ggplot(p = plot_factor_density(all_dm, facet = '~gene_name',
                                                                                                 x='norm_score'),
                                                                         width=14, height=9)


#### Cluster Variants ####
make_var_matrix <- function(x, score='score'){
  variants <- select(x$single_variants, variants, score=!!score) %>%
    mutate(wt = str_sub(variants, end=1),
           mut = str_sub(variants, start=-1),
           pos = str_sub(variants, start=2, end=-2)) %>%
    select(-variants)
    
  if ('=' %in% variants$mut){
    variants$mut[variants$mut=='='] <- variants$wt[variants$mut=='=']
  }
  
  variants <- spread(variants, key = 'mut', value = 'score') %>%
    arrange(pos)
  
  return(variants)
}

combined_variants <- bind_rows(lapply(dms_data, make_var_matrix, score='norm_score'), .id = 'study') %>%
  separate(study, into = c('authour', 'year', 'gene'), sep='_', extra = 'merge') %>%
  select(-`*`, -`_`, -B, -Z, -`<NA>`) %>%
  drop_na(wt, pos)

# Use medians for missing values
per_study_mean_profiles <- group_by(combined_variants, authour, year, gene, wt) %>%
  summarise_at(.vars = vars(-pos), .funs = median, na.rm=TRUE) %>%
  replace(is.na(.), 0)

imputed_variants <- gather(combined_variants, key='mut', value='score', A:Y) %>%
  left_join(., gather(per_study_mean_profiles, key='mut', value='imp_score', A:Y), by = c("authour", "year", "gene", "wt", "mut"))
imputed_variants$score[is.na(imputed_variants$score)] <- imputed_variants$imp_score[is.na(imputed_variants$score)]

imputed_variants <- imputed_variants %>%
  select(-imp_score) %>%
  spread(key=mut, value=score)

pca <- prcomp(as.matrix(select(imputed_variants, -authour, -year, -gene, -wt, -pos)), center = FALSE, scale. = FALSE)

pca_variants <- bind_cols(select(imputed_variants, authour, year, gene, wt, pos), as_tibble(pca$x)) %>%
  unite(col = 'study', authour, year, gene, remove = FALSE)


p_pca <- ggplot(filter(pca_variants, authour %in% c('araya', 'melamed', 'starita', 'kitzman', 'weile')),
                aes(x=PC1, y=PC2, colour=wt)) + geom_point()

#### Save plots #####
save_plot_list(deep_variant_plots, root='figures/variant_analysis/')
