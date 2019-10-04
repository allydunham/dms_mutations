#!/usr/bin/env Rscript 
# Script to analyse clustering of deep mutagenesis position profiles

source('src/config.R')
source('src/analysis/position_profile_clustering.R')

library(dbscan)

meta <- readRDS('data/rdata/study_meta_data.RDS')
variant_matrices <- readRDS('data/rdata/all_study_position_matrices.RDS')
imputed_matrices <- readRDS('data/rdata/all_study_imputed_position_matrices.RDS')
foldx <- readRDS('data/rdata/all_foldx.RDS')
sift_human <- readRDS('data/rdata/human_sift.RDS')
sift_dms <- readRDS('data/rdata/sift_dms_log10.RDS')
backbone_angles <- readRDS('data/rdata/backbone_angles.RDS')
propensity <- readRDS('data/rdata/a_helix_propensity.RDS')
hydrophobicity <- readRDS('data/rdata/hydrophobicity.RDS')
chemical_environments <- readRDS('data/rdata/position_chemical_environments.RDS') %>%
  select(study, position, wt=aa, within_10=within_10.0) %>%
  distinct(study, position, wt, .keep_all = TRUE) %>%
  unnest(cols = within_10) %>%
  mutate(mut = rep(sort(AA_STANDARD), n() / 20)) %>%
  spread(key = mut, value = within_10) %>%
  rename_at(vars(A:Y), ~str_c('chem_env_', .))
data("BLOSUM62")

sift_human_long <- select(sift_human, uniprot_id, pos=position, wt, A:Y) %>%
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

# Conservation
pos_conervation <- filter(sift_human, uniprot_id %in% unique(variant_pcas$norm_all_variants$profiles$uniprot_id)) %>%
  gather('mut', 'sift_score', A:Y) %>%
  filter(!wt == mut) %>%
  group_by(uniprot_id, position, wt) %>%
  summarise(mean_sift = mean(sift_score, na.rm = TRUE),
            mean_log10_sift = mean(log10(sift_score + 0.0001), na.rm = TRUE))

# Cor = -0.36 between sig_count and conservation
plots$pca$norm_all_variants$conservation <- left_join(variant_pcas$norm_all_variants$profiles, rename(pos_conervation, pos=position),
                                                      by=c('uniprot_id', 'pos', 'wt')) %>%
  filter(!is.na(mean_sift)) %>%
  ggplot(aes(x=PC1, y=PC2, colour=mean_log10_sift)) +
  geom_point() +
  theme_pubclean() +
  theme(legend.position = 'right', panel.grid.major  = element_line(linetype = 'dotted', colour = 'grey'))

# Hydrophobicity
plots$pca$norm_all_variants$clean_hydrophobicity <- select(hydrophobicity, wt = aa, hydrophobicity = white4) %>%
  left_join(variant_pcas$norm_all_variants$profiles, ., by='wt') %>%
  ggplot(aes(x=PC1, y=PC2, colour=hydrophobicity)) +
  geom_point() +
  theme_pubclean() +
  theme(legend.position = 'right', panel.grid.major  = element_line(linetype = 'dotted', colour = 'grey'))

# Chemical environment
chem_env <- inner_join(variant_pcas$norm_all_variants$profiles, rename(chemical_environments, pos=position), by=c('study', 'pos', 'wt'))

chem_env_pc_cor <- cor(tibble_to_matrix(chem_env, PC1:PC20), tibble_to_matrix(chem_env, starts_with('chem_env_'))) %>%
  as_tibble(rownames = 'pc') %>%
  gather(key = 'chem_env', value = 'raw_cor', -pc) %>%
  mutate(chem_env = str_sub(chem_env, start = -1),
         pc = factor(pc, levels = str_c('PC', 1:20)))

chem_env_pc_cor <- cor(tibble_to_matrix(chem_env, PC1:PC20),
                       tibble_pca(chem_env, starts_with('chem_env_'))$x) %>%
  as_tibble(rownames = 'pc') %>%
  gather(key = 'chem_env_pc', value = 'pc_cor', -pc) %>%
  mutate(chem_env_pc = factor(chem_env_pc, levels = str_c('PC', 1:20))) %>%
  select(chem_env_pc, pc_cor) %>%
  bind_cols(chem_env_pc_cor, .)

plots$pca$norm_all_variants$chemical_environment <- labeled_ggplot(
  p = ggplot(chem_env_pc_cor, aes(x = pc, y = chem_env, fill = raw_cor)) +
    geom_tile(colour='white') + 
    scale_fill_gradient2() +
    coord_fixed() +
    theme(axis.ticks = element_blank(), panel.background = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    guides(fill = guide_colourbar(title = 'Pearson\nCorrelation')) +
    labs(x = 'Mutational Landscape', y = 'Chemical Environment'),
  units = 'cm', width = 18, height = 14)

plots$pca$norm_all_variants$chemical_environment_pcs <- labeled_ggplot(
  p = ggplot(chem_env_pc_cor, aes(x = pc, y = chem_env_pc, fill = pc_cor)) +
    geom_tile(colour='white') + 
    scale_fill_gradient2() +
    coord_fixed() +
    theme(axis.ticks = element_blank(), panel.background = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    guides(fill = guide_colourbar(title = 'Pearson\nCorrelation')) +
    labs(x = 'Mutational Landscape', y = 'Chemical Environment'),
  units = 'cm', width = 18, height = 14)

# FoldX properties
pcs_foldx <- distinct(foldx, study, position, wt, mut, .keep_all = TRUE) %>%
  group_by(study, pos=position, wt) %>%
  summarise_at(vars(total_energy:entropy_complex), mean) %>%
  inner_join(variant_pcas$norm_all_variants$profiles, ., by = c('study', 'pos', 'wt'))

pcs_foldx_cor <- cor(tibble_to_matrix(pcs_foldx, PC1:PC20),
                     select(pcs_foldx, -water_bridge, -sloop_entropy, -mloop_entropy, -entropy_complex) %>%
                       tibble_to_matrix(total_energy:energy_ionisation)) %>%
  as_tibble(rownames = 'pc') %>%
  gather(key = 'foldx_term', value = 'cor', -pc) %>%
  add_factor_order(pc, foldx_term, cor)

plots$pca$norm_all_variants$foldx_terms <- ggplot(pcs_foldx_cor, aes(x = pc, y = foldx_term, fill = cor)) +
  geom_tile(colour='white') + 
  scale_fill_gradient2() +
  theme(axis.ticks = element_blank(), panel.background = element_blank())

plots$pca$norm_all_variants$foldx_solv_hydrophobic <- ggplot(pcs_foldx, aes(x=PC2, y=solvation_hydrophobic)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  theme_pubclean() +
  theme(legend.position = 'right', panel.grid.major  = element_line(linetype = 'dotted', colour = 'grey'))

plots$pca$norm_all_variants$foldx_solv_hydrophobic <- ggplot(pcs_foldx, aes(x=PC4, y=PC5, colour=entropy_sidechain)) +
  geom_point() +
  scale_colour_gradient2(mid = 'lightgrey', low = 'red', high = 'blue') +
  theme_pubclean() +
  theme(legend.position = 'right', panel.grid.major  = element_line(linetype = 'dotted', colour = 'grey'))

# Phi/Psi
phi_psi_cor <- cor(tibble_to_matrix(variant_pcas$norm_all_variants$profiles, PC1:PC20),
                   tibble_to_matrix(variant_pcas$norm_all_variants$profiles, phi, psi),
                   use = 'pairwise.complete.obs') %>%
  as_tibble(rownames = 'pc') %>%
  gather(key = 'angle', value = 'cor', -pc) %>%
  mutate(pc = factor(pc, levels = str_c('PC', 1:20)))

plots$pca$norm_all_variants$psi_phi <- ggplot(phi_psi_cor, aes(x = pc, y = angle, fill = cor)) +
  geom_tile(colour='white') + 
  scale_fill_gradient2() +
  theme(axis.ticks = element_blank(), panel.background = element_blank())


########

#### tSNE ####
tsne <- tibble_tsne(imputed_matrices$norm_all_variants, A:Y)
  
tsne$tbl <- left_join(tsne$tbl, select(propensity, wt = aa1, a_helix_propensity = exptl), by = 'wt')

plots$tsne <- list()
plots$tsne$gene_confounding <- labeled_ggplot(
  p = tsne_plot(tsne$tbl, gene_name) + guides(colour=guide_legend(title = 'Gene')),
  units = 'cm', height = 14, width = 20
)

plots$tsne$study_confounding <- labeled_ggplot(
  p = tsne_plot(tsne$tbl, study) + guides(colour=guide_legend(title = 'Study')),
  units = 'cm', height = 14, width = 30
)

plots$tsne$helix_propensity <- labeled_ggplot(
  p = tsne_plot(tsne$tbl, a_helix_propensity) + guides(colour=guide_colourbar(title = 'alpha helix propensity')),
  units = 'cm', height = 14, width = 20
)

plots$tsne$sec_struct <- labeled_ggplot(
  p = tsne_plot(tsne$tbl, ss) + guides(colour=guide_legend(title = 'Sec. Struct.')),
  units = 'cm', height = 14, width = 20
)
#######

#### K means clustering ####
n <- 3
kmean_clusters <- group_by(imputed_matrices$norm_sig_positions, wt) %>%
  do(kmean = make_kmeans_clusters(., A:Y, n=n))

kmean_tbl <- map_dfr(kmean_clusters$kmean, .f = ~ .[[1]]) %>%
  mutate(cluster = str_c(wt, '_', cluster))

kmean_analysis <- cluster_analysis(kmean_tbl, backbone_angles = backbone_angles, foldx = rename(foldx, pos=position),
                                   er_str = 'Norm ER', cluster_str = str_c('Kmean, n = ', n), pos_col = pos)

plots$kmean <- kmean_analysis$plots

# Cluster mean sift
kmean_cluster_mean_sift <- select(kmean_tbl, uniprot_id, cluster, pos, wt, A:Y) %>%
  filter(!is.na(uniprot_id)) %>%
  gather(key = 'mut', value = 'er', A:Y) %>%
  left_join(., sift_human_long,
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

# Cluster mean chem env
kmean_cluster_mean_chem_env <- select(kmean_tbl, study, uniprot_id, cluster, position=pos, wt) %>%
  inner_join(chemical_environments, by = c('study', 'position', 'wt')) %>%
  group_by(cluster) %>%
  summarise_at(vars(chem_env_A:chem_env_Y), list(mean = ~mean(.), n = ~n())) %>%
  rename(n=chem_env_A_n) %>%
  select(-(chem_env_C_n:chem_env_Y_n)) %>%
  rename_at(vars(-cluster, -n), ~str_sub(., start = -6, end = -6)) %>%
  gather(key = 'aa', value = 'mean_count', -cluster, -n) %>%
  group_by(aa) %>%
  mutate(norm_count = log2((mean_count + min(mean_count[mean_count > 0])/2)/mean(mean_count))) %>%
  ungroup() %>%
  add_factor_order(cluster, aa, norm_count, sym = FALSE)

plots$kmean$average_chem_env <- ggplot(kmean_cluster_mean_chem_env, aes(x = aa, y = cluster, fill = norm_count)) +
  geom_tile() +
  scale_fill_gradient2() +
  coord_fixed() +
  ggtitle(str_c('Log2(Cluster Mean) score for kmeans (n = ', n, ') clusters')) +
  guides(fill=guide_colourbar(title = 'Norm. Log2(Cluster Mean)')) +
  theme(axis.ticks = element_blank(),
        panel.background = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(colour = AA_COLOURS[unique(kmean_cluster_mean_chem_env$aa)]),
        axis.text.y = element_text(colour = AA_COLOURS[str_sub(levels(kmean_cluster_mean_chem_env$cluster), end = 1)]))
plots$kmean$average_chem_env <- labeled_ggplot(p = plots$kmean$average_chem_env, units = 'cm', width = 18,
                                           height = length(levels(kmean_cluster_mean_chem_env$cluster)) * 0.5 + 3)

########

#### K means no PC1 ####
n <- 3
kmean_no_sig_clusters <- bind_cols(imputed_matrices$norm_all_variants,
                            select(variant_pcas$norm_all_variants$profiles, PC1:PC20)) %>%
  group_by(wt) %>%
  do(kmean = make_kmeans_clusters(., PC2:PC20, n=n))

kmean_no_sig_tbl <- map_dfr(kmean_no_sig_clusters$kmean, .f = ~ .[[1]]) %>%
  mutate(cluster = str_c(wt, '_', cluster))

kmean_no_sig_analysis <- cluster_analysis(kmean_no_sig_tbl, backbone_angles = backbone_angles,
                                          foldx = rename(foldx, pos=position),
                                          er_str = 'Norm ER', cluster_str = str_c('Kmean PC2:PC20, n = ', n), pos_col = pos)

plots$kmean_no_sig <- kmean_no_sig_analysis$plots

kmean_no_sig_er_prof <- group_by(kmean_no_sig_tbl, cluster) %>%
  summarise_at(vars(A:Y), mean) %>%
  gather(key = 'mut', value = 'norm_er', -cluster) %>%
  add_factor_order(cluster, mut, norm_er) %>%
  mutate(cluster = factor(cluster, levels = sort(levels(cluster))))

plots$kmean_no_sig$mean_er_profiles <- labeled_ggplot(
  p=ggplot(kmean_no_sig_er_prof, aes(x=mut, y=cluster, fill=norm_er)) +
    geom_tile() +
    scale_fill_gradient2() +
    coord_fixed() +
    ggtitle('Mean norm er profile for Kmean PC2:PC20, n = 3') +
    guides(fill=guide_colourbar(title = 'Norm ER')) +
    theme(axis.ticks = element_blank(),
          panel.background = element_blank(),
          axis.title = element_blank(),
          axis.text.x = element_text(colour = AA_COLOURS[levels(kmean_no_sig_er_prof$mut)]),
          axis.text.y = element_text(colour = AA_COLOURS[str_sub(levels(kmean_no_sig_er_prof$cluster), end = 1)])),
  units = 'cm', width = 0.5*length(unique(kmean_no_sig_er_prof$mut)) + 4,
  height = 0.5*length(unique(kmean_no_sig_er_prof$cluster)) + 2, limitsize=FALSE)
  
########

#### hclust ####
hclust_settings <- list(h=10, k=NULL, max_k=6, min_k=3)
hclust_clusters <- group_by(imputed_matrices$norm_sig_positions, wt) %>%
  do(hclust = make_hclust_clusters(., A:Y, conf = hclust_settings))
#sapply(hclust_clusters$hclust, function(x){plot(x$hclust); abline(h = hclust_settings$h)})

hclust_tbl <- map_dfr(hclust_clusters$hclust, .f = ~ .[[1]]) %>%
  mutate(cluster = str_c(wt, '_', cluster))

hclust_analysis <- cluster_analysis(hclust_tbl, backbone_angles = backbone_angles, foldx = rename(foldx, pos=position),
                                    er_str = 'Norm ER', cluster_str = make_hclust_cluster_str(hclust_settings), pos_col = pos)
plots$hclust <- hclust_analysis$plots
hclust_analysis$tbl <- hclust_tbl

hclust_cluster_mean_sift <- select(hclust_tbl, uniprot_id, cluster, pos, wt, A:Y) %>%
  filter(!is.na(uniprot_id)) %>%
  gather(key = 'mut', value = 'er', A:Y) %>%
  left_join(., sift_human_long,
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
  ggtitle(str_c('Mean Log10(SIFT) score for ', make_hclust_cluster_str(hclust_settings),' clusters')) +
  guides(fill=guide_colourbar(title = 'log10(SIFT)')) +
  theme(axis.ticks = element_blank(),
        panel.background = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(colour = AA_COLOURS[unique(hclust_cluster_mean_sift$mut)]),
        axis.text.y = element_text(colour = AA_COLOURS[str_sub(levels(hclust_cluster_mean_sift$cluster), end = 1)]))
plots$hclust$average_sift <- labeled_ggplot(p = plots$hclust$average_sift, units = 'cm', width = 20,
                                            height = length(levels(hclust_cluster_mean_sift$cluster)) * 0.5 + 3)
########

#### Hclust all residues ####
hclust_all_settings <- list(h=8, k=NULL, max_k=Inf, min_k=0)
hclust_all <- make_hclust_clusters(imputed_matrices$norm_sig_positions, A:Y, conf = hclust_all_settings)
hclust_all$tbl <- mutate(hclust_all$tbl, cluster = as.character(cluster))

hclust_all_str <- make_hclust_cluster_str(hclust_all_settings)
plots$hclust_all <- list()

# Cluster mean profiles
hclust_all_mean_profiles <- group_by(hclust_all$tbl, cluster) %>%
  summarise_at(.vars = vars(A:Y), .funs = mean)

# Cluster mean profile correlation
hclust_all_cluster_cors <- transpose_tibble(hclust_all_mean_profiles, cluster, name_col = 'aa') %>%
  tibble_correlation(-aa) %>%
  rename(cluster1 = cat1, cluster2 = cat2)

hclust_all_mean_prof_long <- gather(hclust_all_mean_profiles, key='mut', value = 'norm_er', -cluster) %>%
  add_factor_order(cluster, mut, norm_er, sym=FALSE) %>%
  mutate(cluster = factor(cluster, levels = levels(hclust_all_cluster_cors$cluster1)))

# Cluster Sizes
hclust_all_cluster_sizes <- count(hclust_all$tbl, cluster) %>%
  mutate(cluster = factor(cluster, levels = levels(hclust_all_cluster_cors$cluster1)))

hclust_all_plot_height <- 0.5*n_distinct(hclust_all_mean_prof_long$cluster) + 2

# Plot mean profile
plots$hclust_all$mean_prof <- labeled_ggplot(
  p=ggplot(hclust_all_mean_prof_long, aes(x=mut, y=cluster, fill=norm_er)) +
    geom_tile() +
    scale_fill_gradient2() +
    coord_fixed() +
    ggtitle(str_c('Cluster centroid', er_str, 'for', cluster_str, 'clusters', sep=' ')) +
    guides(fill=guide_colourbar(title = er_str)) +
    theme(axis.ticks = element_blank(),
          panel.background = element_blank(),
          axis.title = element_blank()),
  units = 'cm', width = 20, height = hclust_all_plot_height, limitsize=FALSE)

# Plot sizes
plots$hclust_all$cluster_size <- labeled_ggplot(
  p = ggplot(hclust_all_cluster_sizes, aes(x=cluster, y=n)) +
    geom_col(fill = 'cornflowerblue') +
    xlab('Cluster') +
    ylab('Size') +
    scale_y_log10() +
    coord_flip() +
    theme_pubclean() +
    guides(fill=FALSE),
  units = 'cm', height = 15, width = hclust_all_plot_height)

# Plot cors
plots$hclust_all$centre_cor <- labeled_ggplot(
  p = ggplot(hclust_all_cluster_cors, aes(x=cluster1, y=cluster2, fill=cor)) +
    geom_tile() +
    scale_fill_gradient2() +
    ggtitle(str_c('Correlation of', hclust_all_str, 'centroids for clusters based on Norm ER')) +
    coord_fixed() +
    theme(axis.ticks = element_blank(),
          panel.background = element_blank(),
          axis.title = element_blank()),
  units = 'cm', width = hclust_all_plot_height, height = hclust_all_plot_height, limitsize=FALSE)

# FoldX params
hclust_all_fx <- group_by(foldx, study, position, wt) %>%
  summarise_at(.vars = vars(-mut, -pdb_id, -sd), .funs = mean, na.rm=TRUE) %>%
  rename(pos=position) %>%
  inner_join(hclust_all$tbl, ., by=c('study', 'pos', 'wt')) %>%
  select(cluster, study, pos, wt, total_energy:entropy_complex, everything())

hclust_all_fx_cluster_mean <- gather(hclust_all_fx, key = 'foldx_term', value = 'ddG', total_energy:entropy_complex) %>%
  select(cluster, study, pos, wt, foldx_term, ddG, everything()) %>%
  group_by(cluster, foldx_term) %>%
  summarise(ddG = mean(ddG)) %>%
  group_by(foldx_term) %>%
  mutate(max_ddG = max(abs(ddG))) %>%
  filter(max_ddG != 0) %>% # Filter any terms that are all 0
  ungroup() %>%
  mutate(rel_ddG = ddG/max_ddG) %>%
  add_factor_order(cluster, foldx_term, rel_ddG, sym = FALSE) %>%
  mutate(cluster = factor(cluster, levels = levels(hclust_all_cluster_cors$cluster1)))

plots$hclust_all$cluster_avg_foldx_profile <- labeled_ggplot(
  p=ggplot(hclust_all_fx_cluster_mean,
           aes(x=foldx_term, y=cluster, fill=rel_ddG)) +
    geom_tile() +
    scale_fill_gradient2() +
    ggtitle(str_c('Mean FoldX energy terms for each ', hclust_all_str, ' cluster (norm ER)')) +
    coord_fixed() +
    theme(plot.title = element_text(hjust = 0.5, size=8),
          axis.ticks = element_blank(),
          panel.background = element_blank(),
          axis.title = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)),
  units='cm', height=hclust_all_plot_height, width=13, limitsize=FALSE)

# Count of AAs in each cluster
hclust_all_cluster_aa_count <- count(hclust_all$tbl, wt, cluster) %>%
  ungroup() %>%
  spread(key = wt, value = n) %>%
  gather(key = 'wt', value = 'n', -cluster) %>%
  replace_na(list(n=0)) %>%
  add_factor_order(cluster, wt, n) %>%
  group_by(cluster) %>%
  mutate(n = n/sum(n)) %>%
  ungroup() %>%
  mutate(cluster = factor(cluster, levels = hclust_all_analysis$cluster_cor_order))

plots$hclust_all$aa_breakdown <- labeled_ggplot(
  p = ggplot(hclust_all_cluster_aa_count, aes(y=cluster, x=wt, fill=n)) +
    geom_raster() +
    scale_fill_gradientn(colours = c('white', 'blue', 'green', 'yellow', 'red'), values = c(0, 0.25, 0.5, 0.75, 1)) +
    labs(x='Cluster', y='WT',
         title = str_c('Count of WT positions in each ', make_hclust_cluster_str(hclust_all_settings), ' cluster')) +
    theme(axis.ticks = element_blank(),
          panel.background = element_blank(),
          strip.background = element_blank(),
          axis.title = element_blank(),
          plot.title = element_text(hjust = 0.5)),
  units = 'cm', width = 20, height = hclust_all_plot_height)

########

#### HDBSCAN clustering ####
minPts <- 5
hdbscan_clusters <- group_by(imputed_matrices$norm_sig_positions, wt) %>%
  do(hdbscan = make_hdbscan_clusters(., A:Y, minPts = minPts))

hdbscan_tbl <- map_dfr(hdbscan_clusters$hdbscan, .f = ~ .[[1]]) %>%
  mutate(cluster = str_c(wt, '_', cluster))

hdbscan_analysis <- cluster_analysis(hdbscan_tbl, backbone_angles = backbone_angles, foldx = rename(foldx, pos=position),
                                    er_str = 'Norm ER', cluster_str = str_c('HDBSCAN (min = ', minPts, ')'), pos_col = pos)
plots$hdbscan <- hdbscan_analysis$plots
hdbscan_analysis$tbl <- hdbscan_tbl

hdbscan_cluster_mean_sift <- select(hdbscan_tbl, uniprot_id, cluster, pos, wt, A:Y) %>%
  filter(!is.na(uniprot_id)) %>%
  gather(key = 'mut', value = 'er', A:Y) %>%
  left_join(., sift_human_long,
            by = c('uniprot_id', 'pos', 'wt', 'mut')) %>%
  group_by(cluster, wt, mut) %>%
  summarise(er = mean(er, na.rm=TRUE),
            sift = mean(sift, na.rm=TRUE)) %>%
  ungroup() %>%
  mutate(cluster = factor(cluster, levels = hdbscan_analysis$cluster_cor_order))

plots$hdbscan$average_sift <- ggplot(hdbscan_cluster_mean_sift, aes(x = mut, y = cluster, fill = sift)) +
  geom_tile() +
  scale_fill_gradient2() +
  coord_fixed() +
  ggtitle(str_c('Mean Log10(SIFT) score for ', str_c('HDBSCAN (min = ', minPts, ')'),' clusters')) +
  guides(fill=guide_colourbar(title = 'log10(SIFT)')) +
  theme(axis.ticks = element_blank(),
        panel.background = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(colour = AA_COLOURS[unique(hdbscan_cluster_mean_sift$mut)]),
        axis.text.y = element_text(colour = AA_COLOURS[str_sub(levels(hdbscan_cluster_mean_sift$cluster), end = 1)]))
plots$hdbscan$average_sift <- labeled_ggplot(p = plots$hdbscan$average_sift, units = 'cm', width = 20,
                                            height = length(levels(hdbscan_cluster_mean_sift$cluster)) * 0.5 + 3)
########

#### Classifying hdbscan ####
# Version of hdbscan limited to positions with sift data
minPts_sift <- 4
hdbscan_sift_clusters <- semi_join(imputed_matrices$norm_sig_positions, rename(sift_dms, pos=position),
                                   by = c('uniprot_id', 'pos', 'wt')) %>%
  group_by(wt) %>%
  do(hdbscan = make_hdbscan_clusters(., A:Y, minPts = minPts_sift))

hdbscan_sift_tbl <- map_dfr(hdbscan_sift_clusters$hdbscan, .f = ~ .[[1]]) %>%
  mutate(cluster = str_c(wt, '_', cluster))

hdbscan_sift_analysis <- cluster_analysis(hdbscan_sift_tbl, backbone_angles = backbone_angles, foldx = rename(foldx, pos=position),
                                          er_str = 'Norm ER', cluster_str = str_c('HDBSCAN (min = ', minPts, ')'), pos_col = pos)
plots$hdbscan_sift <- hdbscan_sift_analysis$plots
hdbscan_sift_analysis$tbl <- hdbscan_sift_tbl

hdbscan_sift_cluster_mean_sift <- select(hdbscan_sift_tbl, uniprot_id, cluster, pos, wt, A:Y) %>%
  gather(key = 'mut', value = 'er', A:Y) %>%
  left_join(., select(sift_dms, uniprot_id, pos=position, wt, A:Y) %>% gather(key = 'mut', value = 'log10_sift', A:Y),
            by = c('uniprot_id', 'pos', 'wt', 'mut')) %>%
  group_by(cluster, wt, mut) %>%
  summarise(er = mean(er, na.rm=TRUE),
            log10_sift = mean(log10_sift, na.rm=TRUE)) %>%
  ungroup() %>%
  mutate(cluster = factor(cluster, levels = hdbscan_sift_analysis$cluster_cor_order))

plots$hdbscan_sift$average_sift <- ggplot(hdbscan_sift_cluster_mean_sift, aes(x = mut, y = cluster, fill = log10_sift)) +
  geom_tile() +
  scale_fill_gradient2() +
  coord_fixed() +
  ggtitle(str_c('Mean Log10(SIFT) score for ', str_c('HDBSCAN (min = ', minPts, ')'),' clusters')) +
  guides(fill=guide_colourbar(title = 'log10(SIFT)')) +
  theme(axis.ticks = element_blank(),
        panel.background = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(colour = AA_COLOURS[unique(hdbscan_sift_cluster_mean_sift$mut)]),
        axis.text.y = element_text(colour = AA_COLOURS[str_sub(levels(hdbscan_sift_cluster_mean_sift$cluster), end = 1)]))
plots$hdbscan_sift$average_sift <- labeled_ggplot(p = plots$hdbscan_sift$average_sift, units = 'cm', width = 20,
                                                  height = length(levels(hdbscan_sift_cluster_mean_sift$cluster)) * 0.5 + 3)

# Save tsv for classifying
gene_lengths <- structure(meta$gene_length, names=meta$uniprot_id) %>%
  extract(unique(names(.))) %>%
  extract(!is.na(names(.)))

hdbscan_sift_class_tbl <- select(hdbscan_sift_tbl, cluster, uniprot_id, gene_name, position=pos, wt, ss, A:Y) %>%
  filter(uniprot_id %in% unique(sift_dms$uniprot_id)) %>%
  rename_at(vars(A:Y), ~str_c('er_', .)) %>%
  left_join(sift_dms, by = c('uniprot_id', 'position', 'wt')) %>%
  filter(!is.na(species)) %>%
  rename_at(vars(A:Y), ~str_c('sift_', .)) %>%
  mutate(relative_position = position / unname(gene_lengths[uniprot_id]),
         var_sp = 1, var_ss = 1) %>%
  spread(species, var_sp, fill = 0) %>%
  spread(ss, var_ss, fill = 0)

########

#### Tidy plots ####
# Poster Surface Accesibility
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

# Norm Surface Accesibility
plots$pca$norm_all_variants$clean_surf_acc <- labeled_ggplot(
  p=ggplot(drop_na(variant_pcas$norm_all_variants$profiles, all_atom_abs),
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

# Poster Sig Count
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

# Norm Sig Count
plots$pca$norm_all_variants$clean_sig_count <- labeled_ggplot(
  p=ggplot(variant_pcas$norm_all_variants$profiles,
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

# Poster AA profile heatmap
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

# PCA AA vectors
plots$pca$norm_all_variants$aa_vectors <- labeled_ggplot(
  p=ggplot(variant_pcas$norm_all_variants$profiles, aes(x=PC1, y=PC2, colour=wt)) +
  geom_smooth(method = 'lm') + 
  theme_pubclean() +
  scale_colour_manual(values = AA_COLOURS) +
  theme(legend.position = 'right') +
  labs(title = 'Differing trends in slope between positions in PC space based on wt AA'),
  units = 'cm', height = 16, width = 24)

# confounding by study
plots$pca$norm_all_variants$study_confounding <- labeled_ggplot(
  p = select(variant_pcas$norm_all_variants$profiles, study, PC1:PC20) %>%
    gather(key = 'PC', value = 'value', -study) %>%
    mutate(PC = factor(PC, levels = unique(PC)[unique(PC) %>% str_sub(start = 3) %>% as.integer() %>% order()])) %>%
    ggplot(aes(x = value, y = ..ncount.., colour = study)) +
    facet_wrap(~PC, scales = 'free_x') +
    geom_freqpoly() +
    guides(colour=FALSE) +
    labs(title='Distribution of each PC of Norm ER across studies', x = 'PC Value', y = 'Normalised Count') + 
    theme_pubclean() +
    theme(strip.background = element_blank(),
          plot.title = element_text(hjust = 0.5)),
  units = 'cm', height = 20, width = 30
)

########

# Save plots
save_plot_list(plots, root='figures/4_position_profile_clustering/')

# Save analyses
saveRDS(hclust_analysis, 'data/rdata/deep_mut_hclust_clusters.RDS')
