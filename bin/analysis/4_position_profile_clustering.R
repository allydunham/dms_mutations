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

# Corelation between PCs and surface accesibility
pca_surf_cor <- sapply(variant_pcas, pca_surf_acc_cor, simplify = FALSE)

plots <- list_modify(plots, pca=sapply(pca_surf_cor, pca_surface_heatmap, simplify = FALSE))

# Per AA PCAs #
AAs <- sort(unique(variant_matrices$all_variants$wt))
plots <- list_modify(plots, pca=sapply(imputed_matrices, function(x){
  list(per_aa_pcas=sapply(AAs, per_aa_pcas, variant_matrix=x, simplify = FALSE))
}, simplify = FALSE))
########

save_plot_list(plots, root='figures/4_position_profile_clustering/')
