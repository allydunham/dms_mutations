#!/usr/bin/env Rscript 
# Script to analyse clustering of deep mutagenesis position profiles

source('src/config.R')
source('src/analysis/position_profile_clustering.R')

variant_matrices <- readRDS('data/rdata/all_study_position_matrices.RDS')
imputed_matrices <- readRDS('data/rdata/all_study_imputed_position_matrices.RDS')

plots <- list()

#### PCA Analysis ####
variant_pcas <- sapply(imputed_matrices, positional_profile_PCA, simplify = FALSE)

plots$pca <- list()

# plot PCAs against basic covariates
plots$pca <- sapply(variant_pcas, basic_pca_plots, simplify = FALSE)

# Average PCA loading for each AA
avg_AA_pca_profiles <- sapply(variant_pcas, get_avg_aa_pca_profile, simplify = FALSE)

plots <- list_modify(plots, pca=sapply(avg_AA_pca_profiles, aa_avg_profile_plot, simplify = FALSE))

plots <- list_modify(plots, pca=sapply(avg_AA_pca_profiles, aa_profile_heatmap, simplify = FALSE))

plots <- list_modify(plots, pca=sapply(variant_pcas, function(x){list(aa_pca_profile_avg_cor=plot_aa_pca_profile_average_cor(x))},
                                       simplify = FALSE))

# Corelation between PCs and surface accesibility
pca_surf_cor <- sapply(variant_pcas, pca_factor_cor, simplify = FALSE, .vars = vars(all_atom_abs:polar_rel))

plots <- list_modify(plots, pca=sapply(pca_surf_cor, function(x){list(pca_surf_acc_cor_heatmap=pca_factor_heatmap(x))}, simplify = FALSE))

# # Per AA PCAs #
# AAs <- sort(unique(variant_matrices$all_variants$wt))
# plots <- list_modify(plots, pca=sapply(imputed_matrices, function(x){
#   list(per_aa_pcas=sapply(AAs, per_aa_pcas, variant_matrix=x, simplify = FALSE))
# }, simplify = FALSE))
########

#### Tidy plots ####
plots$pca$sig_positions$poster_surf_acc <- ggplot(drop_na(variant_pcas$sig_positions$profiles, all_atom_abs),
                                                  aes(x=PC2, y=PC4, colour=all_atom_abs)) +
  scale_colour_gradientn(colours = c('blue', 'green', 'yellow', 'orange', 'red')) +
  geom_point() +
  lims(x=c(-5,5), y=c(-4,4)) +
  theme_pubclean() +
  theme(legend.position = 'right', panel.grid.major = element_line(colour = 'lightgrey', linetype = 'dotted')) +
  guides(colour=guide_colorbar(title='Surface\nAccesibility'))

plots$pca$sig_positions$poster_sig_count <- ggplot(variant_pcas$sig_positions$profiles,
                                                   aes(x=PC1, y=PC2, colour=sig_count)) +
  scale_colour_gradientn(colours = c('blue', 'red')) +
  geom_point() +
  theme_pubclean() +
  theme(legend.position = 'right', panel.grid.major = element_line(colour = 'lightgrey', linetype = 'dotted')) +
  guides(colour=guide_colorbar(title='# Significant\nSubstitutions'))

# Comparison between avg AA PCA profile and blosum62
data("BLOSUM62")

aa_prof_blosum <- as_tibble(BLOSUM62, rownames='AA1') %>%
  gather(key = 'AA2', value = 'BLOSUM62', -AA1) %>%
  filter(AA1 %in% Biostrings::AA_STANDARD, AA2 %in% Biostrings::AA_STANDARD) %>%
  left_join(., avg_AA_pca_profiles$sig_positions$cor_tbl, by=c('AA1', 'AA2'))

plots$pca$sig_positions$avg_aa_profile_blosum_cor <- ggplot(aa_prof_blosum, aes(x=cor, y=BLOSUM62)) +
  geom_point()

########

save_plot_list(plots, root='figures/4_position_profile_clustering/')
