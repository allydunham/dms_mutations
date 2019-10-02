#!/usr/bin/env Rscript 
# Script to explore distribution of enrichment scores across deep mutagenesis studies

source('src/config.R')
source('src/analysis/enrichment_score_summary.R')

all_variants <- readRDS('data/rdata/all_study_variants.RDS')
variant_matrices <- readRDS('data/rdata/all_study_position_matrices.RDS')
meta_df <- readRDS('data/rdata/study_meta_data.RDS')

plots <- list()

## Plot overall distribution of enrichment scores
plots$dm_hists <- labeled_ggplot(p = plot_study_histogram(all_variants, meta_df),
                                              width = 14, height = 9)

plots$dm_norm_hists <- labeled_ggplot(p = plot_study_histogram(all_variants, meta_df,
                                                                            x='norm_score', thresh = 'norm_thresh'),
                                                   width = 14, height = 9)

plots$dm_norm_hists_tidy <- ggplot(all_variants, aes(x=norm_score, fill=gene_name)) + 
  geom_histogram() +
  geom_vline(aes(xintercept=norm_thresh), colour='red') +
  facet_wrap(~study, scales = 'free_y', ncol = 3, 
             labeller = as_labeller(sapply(unique(all_variants$study), format_study, max_width = 25))) +
  theme_pubclean(base_size = 10) +
  theme(strip.background = element_blank()) +
  guides(fill=FALSE) +
  labs(x = 'Normalised ER', y = 'Count')
plots$dm_norm_hists_tidy <- labeled_ggplot(p = plots$dm_norm_hists_tidy, units = 'cm', height = 26, width = 18)
  
## Plot ER distribution split by various factors
factor_density_config <- list(gene_type=list(facet='~gene_type'), test_class=list(facet='~test_class'),
                              species=list(facet='~species'), gene=list(facet='~gene_name'))
plots$score_factors <- sapply(factor_density_config, function(x){
  labeled_ggplot(
    do.call(plot_factor_density, c(list(tbl=all_variants), x)),
    width=14, height=9
  )
}, simplify = FALSE)

## Plot ER distribution along length of genes
per_position_boxplots <- gather(variant_matrices$all_variants, key='mut', value='score', A:Y) %>%
  select(study, pos, wt, mut, score) %>%
  group_by(study) %>%
  do(plots = plot_score_distribution_per_position(.)) 

plots$per_position_boxplots <- per_position_boxplots$plots
names(plots$per_position_boxplots) <- per_position_boxplots$study

## Save plots
save_plot_list(plots, root='figures/1_study_summaries/')
