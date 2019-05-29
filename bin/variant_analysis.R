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
  bind_rows(.id = 'study') %>%
  select(-chain)

# Dataframe of all individually scored variant/sets of variants in all studies
all_variants <- bind_rows(lapply(dms_data, function(x){x$dm$variant_data}), .id = 'study') %>%
  select(study, variants, score, raw_score, norm_score) %>%
  mutate(variants = str_replace_all(variants, 'p.', '')) %>%
  left_join(., meta_df, by='study')

# Dataframe with matrix of mutational profiles for all positions and all studies
variant_matrix <- bind_rows(lapply(dms_data, make_var_matrix, score='score'), .id = 'study') %>%
  filter(!wt %in% c('Z', 'B')) %>%
  drop_na(wt, pos) %>%
  left_join(., meta_df, by='study') %>%
  select(-factor, -norm_thresh) # Need to switch thresh if using norm_score

imputed_variants <- impute_variant_profiles(variant_matrix)

sig_pos_variant_matrix_filter <- variant_matrix %>%
  mutate_at(.vars = vars(A:Y), .funs = list(~ . < thresh)) %>%
  select(A:Y) %>%
  rowSums(na.rm = TRUE) %>%
  is_greater_than(., 0)

sig_pos_variant_matrix <- filter(variant_matrix, sig_pos_variant_matrix_filter)
sig_pos_imputed_variants <- impute_variant_profiles(sig_pos_variant_matrix, background_matrix = variant_matrix)

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
pca <- prcomp(as.matrix(select(imputed_matrix, A:Y)), center = TRUE, scale. = TRUE)
pca_variants <- bind_cols(select(imputed_matrix, -(A:Y)), as_tibble(pca$x)) %>%
  add_column(sig_pos = sig_pos_variant_matrix_filter)

sig_pca <- prcomp(as.matrix(select(sig_pos_imputed_variants, A:Y)), center = TRUE, scale. = TRUE)
sig_pca_variants <- bind_cols(select(sig_pos_imputed_variants, -(A:Y)), as_tibble(sig_pca$x))
  
deep_variant_plots$pcas <- list()

deep_variant_plots$pcas$tt <- ggplot(pca_variants, aes(x=PC1, y=PC2, colour=sig_pos)) + geom_point()

# PCAs against batchs
deep_variant_plots$pcas$by_authour <- ggplot(pca_variants, aes(x=PC1, y=PC2, colour=study)) + geom_point()

deep_variant_plots$pcas$by_authour_fields <- ggplot(filter(pca_variants,
                                                           authour %in% c('araya', 'melamed', 'starita', 'kitzman', 'weile')),
                                                    aes(x=PC1, y=PC2, colour=study)) + geom_point()

# PCAs against AA types
deep_variant_plots$pcas$by_aa <- ggplot(pca_variants, aes(x=PC1, y=PC2, colour=wt)) + geom_point()

deep_variant_plots$pcas$by_aa_facet <- ggplot(pca_variants, aes(x=PC1, y=PC2, colour=gene_name)) + 
  facet_wrap(~wt) +
  geom_point()

deep_variant_plots$pcas$sig_by_aa_facet <- ggplot(sig_pca_variants, aes(x=PC1, y=PC3, colour=gene_name)) + 
  facet_wrap(~wt) +
  geom_point()

# Average AA PCA loadings
avg_sig_pca_variants <- sig_pca_variants %>%
  group_by(wt) %>%
  summarise_at(.vars = vars(starts_with('PC')), .funs = list(~ mean(.)))

deep_variant_plots$pcas$sig_aa_avg <- ggplot(avg_sig_pca_variants, aes(x=PC2, y=PC3, label=wt)) + geom_text()

avg_sig_pca_variant_cor_mat <- select(avg_sig_pca_variants, -wt) %>% 
  t() %>% 
  set_colnames(avg_sig_pca_variants$wt) %>%
  cor()

aa_order <- rownames(avg_sig_pca_variant_cor_mat)[hclust(dist(avg_sig_pca_variant_cor_mat))$order]

avg_sig_pca_variant_cors <- avg_sig_pca_variant_cor_mat %>%
  as_tibble(rownames = 'AA1') %>%
  gather(key = 'AA2', value = 'cor', -AA1) %>%
  mutate(AA1 = factor(AA1, levels = aa_order),
         AA2 = factor(AA2, levels = aa_order))

deep_variant_plots$pcas$sig_aa_cor_heatmap <- ggplot(avg_sig_pca_variant_cors, aes(x=AA1, y=AA2, fill=cor)) + 
  geom_tile(colour='white') + 
  scale_fill_gradient2() +
  theme(axis.ticks = element_blank(), panel.background = element_blank())

# PCA against surface accesibility
pca_surface <- left_join(pca_variants, rename(all_surface_accesibility, wt=res1), by=c('study', 'wt', 'pos')) %>%
  drop_na(all_atom_abs:polar_rel)

sig_pca_surface <- left_join(sig_pca_variants, rename(all_surface_accesibility, wt=res1), by=c('study', 'wt', 'pos')) %>%
  drop_na(all_atom_abs:polar_rel)

deep_variant_plots$pcas$surface_acc <- ggplot(pca_surface, aes(x=PC1, y=PC2, colour=all_atom_rel)) +
  geom_point() +
  scale_colour_gradientn(colours = c('blue', 'green', 'yellow', 'orange', 'red'))

deep_variant_plots$pcas$surface_acc <- ggplot(pca_surface, aes(x=PC1, y=PC2, colour=all_atom_rel)) +
  geom_point() +
  scale_colour_gradientn(colours = c('blue', 'green', 'yellow', 'orange', 'red'))

# Cor against surface accesibility
sig_pca_variant_mat <- select(sig_pca_surface, starts_with('PC')) %>% 
  set_colnames(str_c('PC', 1:dim(.)[2])) %>%
  as.matrix()

surface_acc_mat <- select(sig_pca_surface, all_atom_abs:polar_rel) %>% 
  set_colnames(colnames(all_surface_accesibility)[-(1:3)]) %>%
  as.matrix()

pca_surface_cor_mat <- cor(sig_pca_variant_mat, surface_acc_mat)

pca_surface_cors <- pca_surface_cor_mat %>%
  as_tibble(rownames = 'PC') %>%
  gather(key = 'Surf', value = 'cor', -PC)

deep_variant_plots$pcas$sig_pc_surf_acc_cor_heatmap <- ggplot(pca_surface_cors, aes(x=PC, y=Surf, fill=cor)) + 
  geom_tile(colour='white') + 
  scale_fill_gradient2() +
  theme(axis.ticks = element_blank(), panel.background = element_blank())


#### Per AA PCAs ####
per_aa_pca <- function(aa, matrix)


#### Save plots #####
save_plot_list(deep_variant_plots, root='figures/variant_analysis/')
