#!/usr/bin/env Rscript 
# Script to analyse clustering of deep mutagenesis position profiles

source('src/config.R')
source('src/analysis/position_profile_clustering.R')

meta <- readRDS('data/rdata/study_meta_data.RDS')
variant_matrices <- readRDS('data/rdata/all_study_position_matrices.RDS')
imputed_matrices <- readRDS('data/rdata/all_study_imputed_position_matrices.RDS')
foldx <- readRDS('data/rdata/all_foldx.RDS')
sift <- readRDS('data/rdata/human_sift.RDS')
backbone_angles <- readRDS('data/rdata/backbone_angles.RDS')
data("BLOSUM62")

sift_long <- select(sift, uniprot_id, pos=position, wt, A:Y) %>%
  filter(uniprot_id %in% unique(variant_matrices$all_variants$uniprot_id)) %>%
  gather(key = 'mut', value = 'sift', A:Y) %>%
  mutate(sift = log10(sift + min(sift[sift > 0], na.rm=TRUE)))

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
pca_factor_cors <- sapply(variant_pcas, pca_factor_cor, simplify = FALSE, .vars = vars(all_atom_abs:polar_rel, psi, phi))

plots <- list_modify(plots, pca=sapply(pca_factor_cors, function(x){list(pca_factor_cor_heatmap=pca_factor_heatmap(x))}, simplify = FALSE))

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

kmean_analysis <- cluster_analysis(kmean_tbl, backbone_angles = backbone_angles, foldx = rename(foldx, pos=position),
                                   er_str = 'Norm ER', cluster_str = str_c('Kmean, n = ', n), pos_col = pos)

plots$kmean <- kmean_analysis$plots

kmean_cluster_mean_sift <- select(kmean_tbl, uniprot_id, cluster, pos, wt, A:Y) %>%
  filter(!is.na(uniprot_id)) %>%
  gather(key = 'mut', value = 'er', A:Y) %>%
  left_join(., sift_long,
            by = c('uniprot_id', 'pos', 'wt', 'mut')) %>%
  group_by(cluster, wt, mut) %>%
  summarise(er = mean(er, na.rm=TRUE),
            sift = mean(sift, na.rm=TRUE)) %>%
  ungroup() %>%
  mutate(cluster = factor(cluster, levels = kmean_analysis$cluster_cor_order))

plots$kmean$average_sift <- ggplot(kmean_cluster_mean_sift, aes(x = mut, y = cluster, fill = sift)) +
  geom_tile() +
  scale_fill_gradient2() +
  coord_fixed() +
  ggtitle(str_c('Mean Log10(SIFT) score for kmeans (n = ', n, ') clusters')) +
  guides(fill=guide_colourbar(title = 'log10(SIFT)')) +
  theme(axis.ticks = element_blank(),
        panel.background = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(colour = AA_COLOURS[unique(kmean_cluster_mean_sift$mut)]),
        axis.text.y = element_text(colour = AA_COLOURS[str_sub(levels(kmean_cluster_mean_sift$cluster), end = 1)]))
plots$kmean$average_sift <- labeled_ggplot(p = plots$kmean$average_sift, units = 'cm', width = 20,
                                           height = length(levels(kmean_cluster_mean_sift$cluster)) * 0.5 + 3)

########

#### hclust ####
h <- 8
hclust_clusters <- group_by(imputed_matrices$norm_sig_positions, wt) %>%
  do(hclust = make_hclust_clusters(., A:Y, h = h))

hclust_tbl <- map_dfr(hclust_clusters$hclust, .f = ~ .[[1]]) %>%
  mutate(cluster = str_c(wt, '_', cluster))

hclust_analysis <- cluster_analysis(hclust_tbl, backbone_angles = backbone_angles, foldx = rename(foldx, pos=position),
                                    er_str = 'Norm ER', cluster_str = str_c('Hclust, h = ', h), pos_col = pos)
plots$hclust <- hclust_analysis$plots

hclust_cluster_mean_sift <- select(hclust_tbl, uniprot_id, cluster, pos, wt, A:Y) %>%
  filter(!is.na(uniprot_id)) %>%
  gather(key = 'mut', value = 'er', A:Y) %>%
  left_join(., sift_long,
            by = c('uniprot_id', 'pos', 'wt', 'mut')) %>%
  group_by(cluster, wt, mut) %>%
  summarise(er = mean(er, na.rm=TRUE),
            sift = mean(sift, na.rm=TRUE)) %>%
  ungroup() %>%
  mutate(cluster = factor(cluster, levels = hclust_analysis$cluster_cor_order))

plots$hclust$average_sift <- ggplot(hclust_cluster_mean_sift, aes(x = mut, y = cluster, fill = sift)) +
  geom_tile() +
  scale_fill_gradient2() +
  coord_fixed() +
  ggtitle(str_c('Mean Log10(SIFT) score for hclust (h = ', h, ') clusters')) +
  guides(fill=guide_colourbar(title = 'log10(SIFT)')) +
  theme(axis.ticks = element_blank(),
        panel.background = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(colour = AA_COLOURS[unique(hclust_cluster_mean_sift$mut)]),
        axis.text.y = element_text(colour = AA_COLOURS[str_sub(levels(hclust_cluster_mean_sift$cluster), end = 1)]))
plots$hclust$average_sift <- labeled_ggplot(p = plots$hclust$average_sift, units = 'cm', width = 20,
                                            height = length(levels(hclust_cluster_mean_sift$cluster)) * 0.5 + 3)
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
aa_prof_blosum <- as_tibble(BLOSUM62, rownames='AA1') %>%
  gather(key = 'AA2', value = 'BLOSUM62', -AA1) %>%
  filter(AA1 %in% Biostrings::AA_STANDARD, AA2 %in% Biostrings::AA_STANDARD) %>%
  left_join(., mutate(avg_AA_pca_profiles$sig_positions$cor_tbl, AA1 = as.character(AA1), AA2 = as.character(AA2)),
            by=c('AA1', 'AA2'))

plots$pca$sig_positions$avg_aa_profile_blosum_cor <- ggplot(aa_prof_blosum, aes(x=cor, y=BLOSUM62)) +
  geom_point()

########

# Save plots
save_plot_list(plots, root='figures/4_position_profile_clustering/')

# Save analyses
saveRDS(hclust_analysis, str_c('data/rdata/deep_mut_hclust_clusters.RDS'))
