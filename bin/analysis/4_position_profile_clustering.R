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

# Average correlation between PCA profiles for each AA (rather than cor of average profile)
# Long and not particularly exciting, don't regenerate until needed
#plots <- list_modify(plots, pca=sapply(variant_pcas, function(x){list(aa_pca_profile_avg_cor=plot_aa_pca_profile_average_cor(x))},
#                                       simplify = FALSE))

# Corelation between PCs and surface accesibility
pca_surf_cor <- sapply(variant_pcas, pca_factor_cor, simplify = FALSE, .vars = vars(all_atom_abs:polar_rel))

plots <- list_modify(plots, pca=sapply(pca_surf_cor, function(x){list(pca_surf_acc_cor_heatmap=pca_factor_heatmap(x))}, simplify = FALSE))

# # Per AA PCAs #
# AAs <- sort(unique(variant_matrices$all_variants$wt))
# plots <- list_modify(plots, pca=sapply(imputed_matrices, function(x){
#   list(per_aa_pcas=sapply(AAs, per_aa_pcas, variant_matrix=x, simplify = FALSE))
# }, simplify = FALSE))
########

#### K means clustering ####
n <- 3
kmean_clusters <- group_by(imputed_matrices$norm_sig_positions, wt) %>%
  do(kmean = make_kmeans_clusters(., A:Y, n=n))

kmean_tbl <- map_dfr(kmean_clusters$kmean, .f = ~ .[[1]]) %>%
  mutate(cluster = str_c(wt, '_', cluster))

kmean_cluster_summaries <- group_by(kmean_tbl, cluster) %>%
  summarise_at(., .vars = vars(A:Y), .funs = mean)

kmean_cluster_cors <- transpose_tibble(kmean_cluster_summaries, cluster, name_col = 'aa') %>%
  tibble_correlation(-aa) %>%
  rename(cluster1 = cat1, cluster2 = cat2) %>%
  mutate(wt1 = str_sub(cluster1, end = 1),
         wt2 = str_sub(cluster2, end = 1))

plots$kmeans$kmeans_centre_cor <- labeled_ggplot(
  p = ggplot(kmean_cluster_cors, aes(x=cluster1, y=cluster2, fill=cor)) +
    geom_tile() +
    scale_fill_gradient2() +
    ggtitle(str_c('Correlation of Kmeans cluster centroids (n = ', n, ') for clusters of AA position based on ER')) + 
    coord_fixed() +
    theme(axis.ticks = element_blank(),
          panel.background = element_blank(),
          axis.title = element_blank(),
          axis.text.x = element_text(colour = AA_COLOURS[str_sub(levels(kmean_cluster_cors$cluster1), end = 1)],
                                     angle = 90, vjust = 0.5),
          axis.text.y = element_text(colour = AA_COLOURS[str_sub(levels(kmean_cluster_cors$cluster2), end = 1)])),
  width = 10, height = 10)
########

#### hclust ####
h <- 18
hclust_clusters <- group_by(imputed_matrices$norm_sig_positions, wt) %>%
  do(hclust = make_hclust_clusters(., A:Y, h = h))

hclust_tbl <- map_dfr(hclust_clusters$hclust, .f = ~ .[[1]]) %>%
  mutate(cluster = str_c(wt, '_', cluster))

hclust_cluster_summaries <- group_by(hclust_tbl, cluster) %>%
  summarise_at(., .vars = vars(A:Y), .funs = mean)

hclust_cluster_cors <- transpose_tibble(hclust_cluster_summaries, cluster, name_col = 'aa') %>%
  tibble_correlation(-aa) %>%
  rename(cluster1 = cat1, cluster2 = cat2) %>%
  mutate(wt1 = str_sub(cluster1, end = 1),
         wt2 = str_sub(cluster2, end = 1))

plots$kmeans$hclust_centre_cor <- labeled_ggplot(
  p = ggplot(hclust_cluster_cors, aes(x=cluster1, y=cluster2, fill=cor)) +
    geom_tile() +
    scale_fill_gradient2() +
    ggtitle(str_c('Correlation of hclust cluster centroids (h = ', h, ') for clusters of AA position based on ER')) + 
    coord_fixed() +
    theme(axis.ticks = element_blank(),
          panel.background = element_blank(),
          axis.title = element_blank(),
          axis.text.x = element_text(colour = AA_COLOURS[str_sub(levels(hclust_cluster_cors$cluster1), end = 1)],
                                     angle = 90, vjust = 0.5),
          axis.text.y = element_text(colour = AA_COLOURS[str_sub(levels(hclust_cluster_cors$cluster2), end = 1)])),
  width = 10, height = 10)

########

#### Tidy plots ####
plots$pca$sig_positions$poster_surf_acc <- labeled_ggplot(
  p=ggplot(drop_na(variant_pcas$sig_positions$profiles, all_atom_abs),
           aes(x=PC2, y=PC4, colour=all_atom_abs)) +
    scale_colour_gradientn(colours = c('blue', 'green', 'yellow', 'orange', 'red')) +
    geom_point() +
    lims(x=c(-5,5), y=c(-4,4)) +
    guides(colour=guide_colorbar(title='Surface\nAccesibility')) +
    theme_pubclean() +
    theme(legend.position = 'right',
          panel.grid.major = element_line(colour = 'lightgrey', linetype = 'dotted'),
          axis.title = element_text(size=30),
          axis.text = element_text(size=20),
          legend.text = element_text(size=20),
          legend.title = element_text(size=30)),
  width = 9.82, height = 7.09)


plots$pca$sig_positions$poster_sig_count <- labeled_ggplot(
  p=ggplot(variant_pcas$sig_positions$profiles,
                                                   aes(x=PC1, y=PC2, colour=sig_count)) +
  scale_colour_gradientn(colours = c('blue', 'red')) +
  geom_point() +
  guides(colour=guide_colorbar(title='# Significant\nSubstitutions')) +  
  theme_pubclean() +
  theme(legend.position = 'right', panel.grid.major = element_line(colour = 'lightgrey', linetype = 'dotted'),
        axis.title = element_text(size=30),
        axis.text = element_text(size=20),
        legend.text = element_text(size=20),
        legend.title = element_text(size=30)),
width = 9.82, height = 7.09)

plots$pca$sig_positions$poster_aa_profile_heatmap <- labeled_ggplot(
  p=ggplot(avg_AA_pca_profiles$sig_positions$cor_tbl, aes(x=AA1, y=AA2, fill=cor)) + 
    geom_tile(colour='white') + 
    scale_fill_gradient2() +
    coord_fixed() + 
    guides(fill=guide_colourbar(title='Correlation')) +
    theme(axis.ticks = element_blank(),
          panel.background = element_blank(),
          axis.text = element_text(size=20),
          axis.title = element_blank(),
          legend.text = element_text(size=15),
          legend.title = element_text(size=20)),
  width = 8, height = 6)

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
