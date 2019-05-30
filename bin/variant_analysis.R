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
                  group = sapply(dms_data, function(x){get_meta(x$dm, 'group')}),
                  thresh=sapply(dms_data, function(x){x$manual_threshold}),
                  factor=sapply(dms_data, function(x){x$norm_factor}),
                  norm_thresh=thresh/factor)

# Surface accessibility scores for all variants
surface_accesibility <- sapply(dms_data, function(x){if(!identical(NA, x$surface_accesibility)){x$surface_accesibility$combined}},
                                   simplify = FALSE) %>%
  bind_rows(.id = 'study')

# Dataframe of all individually scored variant/sets of variants in all studies
all_variants <- bind_rows(lapply(dms_data, function(x){x$dm$variant_data}), .id = 'study') %>%
  select(study, variants, score, raw_score, norm_score) %>%
  mutate(variants = str_replace_all(variants, 'p.', '')) %>%
  left_join(., meta_df, by='study')

# Dataframe with matrix of mutational profiles for all positions and all studies
variant_matrices <- list()
variant_matrices$all_variants <- bind_rows(lapply(dms_data, make_var_matrix, score='score'), .id = 'study') %>%
  filter(!wt %in% c('Z', 'B')) %>%
  drop_na(wt, pos) %>%
  left_join(., meta_df, by='study') %>%
  left_join(., rename(surface_accesibility, wt=res1), by=c('study', 'wt', 'pos')) %>%
  select(-factor, -norm_thresh) # Need to switch thresh if using norm_score

variant_matrices$all_variants <- mutate(variant_matrices$all_variants,
                                        sig_count = variant_matrices$all_variants %>%
                                          mutate_at(.vars = vars(A:Y), .funs = list(~ . < thresh)) %>%
                                          select(A:Y) %>%
                                          rowSums(na.rm = TRUE)
)

variant_matrices$sig_positions <- filter(variant_matrices$all_variants, sig_count > 0)

imputed_matrices <- sapply(variant_matrices, impute_variant_profiles, background_matrix=variant_matrices$all_variants, simplify=FALSE)

#### Enrichment Score Distributions ####
deep_variant_plots$dm_hists <- labeled_ggplot(p = plot_study_histogram(all_variants, meta_df),
                                              width = 14, height = 9)

deep_variant_plots$dm_norm_hists <- labeled_ggplot(p = plot_study_histogram(all_variants, meta_df,
                                                                            x='norm_score', thresh = 'norm_thresh'),
                                                   width = 14, height = 9)

factor_density_config <- list(gene_type=list(facet='~gene_type'), test_class=list(facet='~test_class'),
                              species=list(facet='~species'), gene=list(facet='~gene_name'))
deep_variant_plots$score_factor_densities <- sapply(factor_density_config, function(x){
  labeled_ggplot(
    do.call(plot_factor_density, c(list(tbl=all_variants), x)),
    width=14, height=9
    )
}, simplify = FALSE)

per_position_boxplots <- gather(variant_matrices$all_variants, key='mut', value='score', A:Y) %>%
  select(study, pos, wt, mut, score) %>%
  group_by(study) %>%
  do(plots = plot_score_distribution_per_position(.)) 

deep_variant_plots$score_distributions$per_position_boxplots <- per_position_boxplots$plots
names(deep_variant_plots$score_distributions$per_position_boxplots) <- per_position_boxplots$study
########

#### Average severity vs Blosum ####
data(BLOSUM62)
sub_scores <- variant_matrices$all_variants %>%
  select(wt, A:Y) %>%
  group_by(wt) %>%
  summarise_all(.funs = mean, na.rm=TRUE) %>%
  gather(key = 'mut', value = 'er', -wt) %>%
  left_join(., BLOSUM62 %>%
              as_tibble(rownames = 'wt') %>%
              gather(key='mut', value = 'blosum62', -wt),
            by=c('wt', 'mut'))

deep_variant_plots$avg_scores <- list()
deep_variant_plots$avg_scores$er_vs_blosum <- ggplot(sub_scores, aes(x=er, y=blosum62)) + geom_point()
  
########

#### Variant Score PCAs ####
variant_pcas <- sapply(imputed_matrices, positional_profile_PCA, simplify = FALSE)
  
deep_variant_plots$pcas <- list()

# plot PCAs against basic covariates
deep_variant_plots$pcas <- sapply(variant_pcas, basic_pca_plots, simplify = FALSE)

# Average PCA loading for each AA
avg_AA_pca_profiles <- sapply(variant_pcas, get_avg_aa_pca_profile, simplify = FALSE)

deep_variant_plots <- list_modify(deep_variant_plots, pcas=sapply(avg_AA_pca_profiles, aa_avg_profile_plot, simplify = FALSE))

deep_variant_plots <- list_modify(deep_variant_plots, pcas=sapply(avg_AA_pca_profiles, aa_profile_heatmap, simplify = FALSE))

# Corelation between PCs and surface accesibility
pca_surf_cor <- sapply(variant_pcas, pca_surf_acc_cor, simplify = FALSE)

deep_variant_plots <- list_modify(deep_variant_plots, pcas=sapply(pca_surf_cor, pca_surface_heatmap, simplify = FALSE))

# Per AA PCAs #
AAs <- sort(unique(variant_matrices$all_variants$wt))
deep_variant_plots <- list_modify(deep_variant_plots, pcas=sapply(imputed_matrices, function(x){
  list(per_aa_pcas=sapply(AAs, per_aa_pcas, variant_matrix=x, simplify = FALSE))
}, simplify = FALSE))

#### Profile vs SIFT/FoldX ####
prediction_martices <- list()
prediction_martices$sift <- bind_rows(lapply(dms_data, make_var_matrix, score='sift_score'), .id = 'study') %>%
  mutate(sift = TRUE)
prediction_martices$foldx <- bind_rows(lapply(dms_data, make_foldx_var_matrix), .id = 'study')  %>%
  mutate(foldx = TRUE)


var_preds_mat <- variant_matrices$all_variants %>%
  left_join(., prediction_martices$foldx, suffix = c('', '_foldx'), by=c('study', 'pos', 'wt')) %>%
  left_join(., prediction_martices$sift, suffix = c('', '_sift'), by=c('study', 'pos', 'wt')) %>%
  mutate(foldx = !is.na(foldx),
         sift = !is.na(sift))

# FoldX
foldx_dist <- pdist(var_preds_mat %>%
                      filter(foldx) %>%
                      select(A:Y) %>%
                      as.matrix(),
                    var_preds_mat %>%
                      filter(foldx) %>%
                      select(A_foldx:Y_foldx) %>%
                      as.matrix())

# FoldX
sift_dist <- pdist(var_preds_mat %>%
                      filter(sift) %>%
                      select(A:Y) %>%
                      as.matrix(),
                    var_preds_mat %>%
                      filter(sift) %>%
                      select(A_sift:Y_sift) %>%
                      as.matrix())


#### Save plots #####
save_plot_list(deep_variant_plots, root='figures/variant_analysis/')
