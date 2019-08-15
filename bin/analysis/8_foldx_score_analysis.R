#!/usr/bin/env Rscript 
# Perform complimentary analysis on FoldX scores to that done on DMS data

source('src/config.R')
source('src/analysis/foldx_score_analysis.R')
source('src/analysis/position_profile_clustering.R')

foldx <- readRDS('data/rdata/human_foldx_tiny.RDS')

foldx_ptms <- readRDS('data/rdata/human_foldx_ptms.RDS')

deep_mut_hclust <- readRDS('data/rdata/deep_mut_hclust_clusters.RDS')

plots <- list()

#### Analyse by average at position ####
# foldx_pos_avg <- group_by(foldx, uniprot_id, position, wt) %>%
#   summarise_at(vars(total_energy:entropy_complex), mean, na.rm=TRUE) %>%
#   ungroup()
# saveRDS(foldx_pos_avg, 'data/rdata/human_foldx_tiny_pos_avg.RDS')
foldx_pos_avg <- readRDS('data/rdata/human_foldx_tiny_pos_avg.RDS')

# drop unused factors then center and scale (pssobly add water_bridge/kon back if structures where they're meaningful are used later)
foldx_pos_avg_scaled <- select(foldx_pos_avg, -total_energy, -sloop_entropy, -mloop_entropy,
                               -entropy_complex, -electrostatic_kon, -water_bridge) %>% 
  tibble_to_matrix(., backbone_hbond:energy_ionisation) %>% 
  scale(center = FALSE, scale = TRUE) %>% 
  as_tibble() %>%
  bind_cols(select(foldx_pos_avg, uniprot_id:wt), .)

pos_avg_settings <- list(h=12.5, k=NULL, max_k=6, min_k=3)
pos_avg_hclust <- group_by(foldx_pos_avg_scaled, wt) %>%
  do(hclust = make_hclust_clusters(., backbone_hbond:energy_ionisation, conf = pos_avg_settings))
#sapply(pos_avg_hclust$hclust, function(x){plot(x$hclust); abline(h = pos_avg_settings$h)})

pos_avg_hclust_tbl <- map_dfr(pos_avg_hclust$hclust, .f = ~ .[[1]]) %>%
  mutate(cluster = str_c(wt, '_', cluster))

pos_avg_hclust_analysis <- analyse_clusters(pos_avg_hclust_tbl, make_hclust_cluster_str(pos_avg_settings), 'mean FoldX terms',
                                            backbone_hbond:energy_ionisation,
                                            transform_ddg = function(x){x / max(abs(x))})

plots$pos_avg_clusters <- pos_avg_hclust_analysis$plots
########

#### Analyse by total energy of each substitution ####
foldx_total_energy <- select(foldx, uniprot_id:mut, total_energy) %>%
  spread(key = mut, value = total_energy)

tot_eng_settings <- list(h=75, k=NULL, max_k=6)
tot_eng_hclust <- group_by(foldx_total_energy, wt) %>%
  do(hclust = make_hclust_clusters(., A:Y, conf = tot_eng_settings))
#sapply(tot_eng_hclust$hclust, function(x){plot(x$hclust); abline(h = tot_eng_settings$h)})

tot_eng_hclust_tbl <- map_dfr(tot_eng_hclust$hclust, .f = ~ .[[1]]) %>%
  mutate(cluster = str_c(wt, '_', cluster))

tot_eng_hclust_analysis <- analyse_clusters(tot_eng_hclust_tbl, make_hclust_cluster_str(tot_eng_settings),
                                            'per substitution total_energy', A:Y,
                                            transform_ddg = function(x){atan(0.5 * x)})
tot_eng_hclust_analysis$plots$mean_profile$plot <- tot_eng_hclust_analysis$plots$mean_profile$plot + 
  theme(axis.text.x = element_text(angle = 0))

plots$total_energy <- tot_eng_hclust_analysis$plots
########

#### Analyse by total energy of each substitution use PTMs dataset ####
foldx_ptms_total_energy <- select(foldx_ptms, uniprot_id:mut, modified, acetylation:ubiquitination, total_energy) %>%
  spread(key = mut, value = total_energy)

tot_eng_ptms_settings <- list(h=75, k=NULL, max_k=6)
tot_eng_ptms_hclust <- group_by(foldx_ptms_total_energy, wt) %>%
  do(hclust = make_hclust_clusters(., A:Y, conf = tot_eng_settings))
#sapply(tot_eng_ptms_hclust$hclust, function(x){plot(x$hclust); abline(h = tot_eng_ptms_settings$h)})

tot_eng_ptms_hclust_tbl <- map_dfr(tot_eng_ptms_hclust$hclust, .f = ~ .[[1]]) %>%
  mutate(cluster = str_c(wt, '_', cluster))

tot_eng_ptms_hclust_analysis <- analyse_clusters(tot_eng_ptms_hclust_tbl, make_hclust_cluster_str(tot_eng_ptms_settings),
                                            'per substitution total_energy', A:Y,
                                            transform_ddg = function(x){atan(0.5 * x)})
tot_eng_ptms_hclust_analysis$plots$mean_profile$plot <- tot_eng_ptms_hclust_analysis$plots$mean_profile$plot + 
  theme(axis.text.x = element_text(angle = 0))

plots$total_energy_ptms <- tot_eng_ptms_hclust_analysis$plots

# Mean profiles for each cluster under each modification
tot_eng_ptm_mean_profs <- mutate(tot_eng_ptms_hclust_tbl, none = !modified) %>%
  gather(key = 'modification', value = 'present', none, acetylation:ubiquitination) %>%
  filter(present) %>%
  group_by(cluster, modification) %>%
  summarise_at(vars(A:Y), mean) %>%
  gather(key = 'aa', value = 'ddg', A:Y) %>%
  spread(modification, ddg) %>%
  gather(key = 'ptm', value = 'ddg', -cluster, -aa, -none) %>%
  drop_na() %>%
  mutate(cluster_aa = str_sub(cluster, end=1),
         cluster_num = str_sub(cluster, start=3)) %>%
  filter(!ptm %in% c('o_galnac_glycosilation', 'o_glcnac_glycosilation'))

plots$total_energy_ptms$ptm_profile_cor <- ggplot(tot_eng_ptm_mean_profs, aes(x=none, y=ddg, shape=ptm, colour=aa)) +
  geom_point() +
  scale_colour_manual(values = AA_COLOURS) +
  facet_wrap(~cluster, scales = 'free') +
  theme_pubclean()
########

#### Analyse all FoldX terms ####
foldx_wide <- select(foldx, -sd, -total_energy, -sloop_entropy, -mloop_entropy, -entropy_complex, -electrostatic_kon, -water_bridge) %>%
  gather(key = 'term', value = 'ddg', backbone_hbond:energy_ionisation) %>%
  unite(term, mut, term) %>%
  spread(term, ddg)

foldx_wide_scaled <- tibble_to_matrix(foldx_wide, A_backbone_clash:Y_van_der_waals_clashes) %>% 
  scale(center = FALSE, scale = TRUE) %>% 
  as_tibble() %>%
  bind_cols(select(foldx_wide, uniprot_id:sumoylation), .)

all_terms_settings <- list(h=300, k=NULL, max_k=6, min_k=3)
all_terms_hclust <- group_by(foldx_wide_scaled, wt) %>%
  do(hclust = make_hclust_clusters(., A_backbone_clash:Y_van_der_waals_clashes, conf = all_terms_settings))
#sapply(all_terms_hclust$hclust, function(x){plot(x$hclust); abline(h = all_terms_settings$h)})

all_terms_hclust_tbl <- map_dfr(all_terms_hclust$hclust, .f = ~ .[[1]]) %>%
  mutate(cluster = str_c(wt, '_', cluster))

all_terms_hclust_analysis <- analyse_clusters(all_terms_hclust_tbl, make_hclust_cluster_str(all_terms_settings),
                                              'all FoldX terms for all substitutions',
                                              A_backbone_clash:Y_van_der_waals_clashes,
                                              transform_ddg = function(x){x / max(abs(x))})
all_terms_hclust_analysis$plots$mean_profile$width <- 100

plots$all_terms <- all_terms_hclust_analysis$plots

# Correlation of terms between amino acid
term_correlation <- select(all_terms_hclust_analysis$mean_profiles_long, cluster, term, ddG) %>%
  spread(term, ddG) %>%
  tibble_correlation(-cluster) %>%
  rename(term1 = cat1, term2 = cat2)

plots$all_terms$term_correlations <- labeled_ggplot(
  p = ggplot(term_correlation, aes(x=term1, y=term2, fill=cor)) +
    geom_raster() +
    scale_fill_gradient2() +
    labs(x='', y='', title = 'Correlation of each FoldX term accross clusters') +
    coord_fixed() +
    theme(axis.ticks = element_blank(),
          panel.background = element_blank(),
          axis.title = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          plot.title = element_text(hjust = 0.5)),
  units='cm',
  height = 100,
  width = 100)

term_correlation_reduced <- mutate(term_correlation,
                                   aa1 = str_sub(term1, end = 1), term1 = str_sub(term1, start = 3),
                                   aa2 = str_sub(term2, end = 1), term2 = str_sub(term2, start = 3)) %>%
  filter(term1 == term2) %>%
  select(term = term1, cor, aa1, aa2)

plots$all_terms$term_correlation_reduced <- labeled_ggplot(
  p = ggplot(term_correlation_reduced, aes(x=aa1, y=aa2, fill=cor)) +
    facet_wrap(~term) +
    geom_raster() +
    scale_fill_gradient2() +
    labs(x='', y='', title = 'Correlation of each AA within FoldX terms') +
    coord_fixed() +
    theme(axis.ticks = element_blank(),
          panel.background = element_blank(),
          axis.title = element_blank(),
          plot.title = element_text(hjust = 0.5)),
  units='cm',
  height = 100,
  width = 100)
########

#### Compare to Deep Mut Clustering ####
deep_mut_cor <- cor(select(deep_mut_hclust$foldx_cluster_mean_energy, cluster, term=foldx_term, ddg=ddG) %>%
                      filter(term %in% colnames(pos_avg_hclust_analysis$mean_profiles)) %>%
                      mutate(cluster = str_c('dm_', cluster)) %>%
                      spread(key = term, value = ddg) %>%
                      tibble_to_matrix(-cluster, row_names = 'cluster') %>%
                      t() %>%
                      alphabetise_matrix(by='rows'),
                    mutate(pos_avg_hclust_analysis$mean_profiles, cluster = str_c('fx_', cluster)) %>%
                      tibble_to_matrix(-cluster, row_names = 'cluster') %>%
                      t() %>%
                      alphabetise_matrix(by='rows')) %>%
  as_tibble(rownames = 'dm_cluster') %>%
  gather(key = 'fx_cluster', value = 'cor', -dm_cluster) %>%
  mutate(dm_aa = str_sub(dm_cluster, start=4, end=4),
         fx_aa = str_sub(fx_cluster, start=4, end=4)) %>%
  add_factor_order(dm_cluster, fx_cluster, cor, sym=FALSE)

plots$pos_avg_clusters$dm_comparison_cor <- labeled_ggplot(
  p = filter(deep_mut_cor, dm_aa == fx_aa) %>%
  ggplot(aes(x=dm_cluster, y=fx_cluster, fill=cor)) +
  geom_raster() +
  facet_wrap(~dm_aa, nrow = 4, ncol = 5, scales = 'free') +
  scale_fill_gradient2() +
  labs(x='', y='', title = 'Comparison between Deep Mut and FoldX based hclust using pearsons correlation') +
  theme(axis.ticks = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5)),
  units = 'cm', height = 20, width = 30) 

deep_mut_dist <- bind_rows(select(deep_mut_hclust$foldx_cluster_mean_energy, cluster, term=foldx_term, ddg=ddG) %>%
                             filter(term %in% colnames(pos_avg_hclust_analysis$mean_profiles)) %>%
                             mutate(cluster = str_c('dm_', cluster)) %>%
                             spread(key = term, value = ddg),
                           mutate(pos_avg_hclust_analysis$mean_profiles, cluster = str_c('fx_', cluster))) %>%
  tibble_to_matrix(-cluster, row_names = 'cluster') %>%
  dist() %>%
  as.matrix() %>%
  as_tibble(rownames = 'dm_cluster') %>%
  gather(key = 'fx_cluster', value = 'dist', -dm_cluster) %>%
  filter(startsWith(dm_cluster, 'dm_'), startsWith(fx_cluster, 'fx_')) %>%
  mutate(dm_aa = str_sub(dm_cluster, start=4, end=4),
         fx_aa = str_sub(fx_cluster, start=4, end=4)) %>%
  add_factor_order(dm_cluster, fx_cluster, dist, sym=FALSE)

plots$pos_avg_clusters$dm_comparison_dist <- labeled_ggplot(
  p = filter(deep_mut_dist, dm_aa == fx_aa) %>%
    ggplot(aes(x=dm_cluster, y=fx_cluster, fill=atan(0.5*dist) )) +
    geom_raster() +
    facet_wrap(~dm_aa, nrow = 4, ncol = 5, scales = 'free') +
    scale_fill_viridis_c() +
    labs(x='', y='', title = 'Comparison between Deep Mut and FoldX based hclust using euclidean distance') +
    theme(axis.ticks = element_blank(),
          panel.background = element_blank(),
          strip.background = element_blank(),
          axis.title = element_blank(),
          plot.title = element_text(hjust = 0.5)),
  units = 'cm', height = 20, width = 30) 

########

#### Assign proteome FoldX results to Deep Mut clusters ####
# NOTE: not all deep mut clusters include positions in structures so we don't have average FoldX profiles for all of them
plots$dm_cluster_assignment <- list()

# Generate matrix of foldx position profiles and scaling factors for each term
foldx_mat <- select(foldx, total_energy:entropy_complex) %>%
  select(-total_energy, -sloop_entropy, -mloop_entropy, -electrostatic_kon, -entropy_complex, -water_bridge) %>% # drop unused/all 0 terms
  as.matrix() %>%
  alphabetise_matrix(by = 'columns') %>%
  scale(., center = FALSE, scale = TRUE)

term_scaling_factors <- attr(foldx_mat, "scaled:scale")
foldx_scaled_pos_avg_mat <- bind_cols(select(foldx, uniprot_id, position, wt), as_tibble(foldx_mat)) %>%
  group_by(uniprot_id, position, wt) %>%
  summarise_all(mean) %>%
  ungroup()

# Generate scaled profile matrix for DM derived clusters
dm_cluster_foldx_scaled_profiles <- select(deep_mut_hclust$foldx, cluster, total_energy:entropy_complex) %>%
  select(-total_energy, -sloop_entropy, -mloop_entropy, -electrostatic_kon, -entropy_complex, -water_bridge) %>% # drop unused/all 0 terms
  tibble_to_matrix(-cluster, row_names = 'cluster') %>%
  alphabetise_matrix(by = 'both') %>%
  scale(., center = FALSE, scale = term_scaling_factors) %>% # Scale by factors calculated for large FoldX set
  as_tibble(rownames = 'cluster') %>%
  group_by(cluster) %>%
  summarise_all(mean) %>%
  tibble_to_matrix(-cluster, row_names = 'cluster') %>%
  t() # matrix ops done by column

# return a vector of distances between query vector and columns of ref
col_distances <- function(query, ref){
  sqrt(colSums((ref - query)^2))
}

foldx_cluster_dists <- select(foldx_scaled_pos_avg_mat, -uniprot_id, -position, -wt) %>%
  apply(1, col_distances, ref=dm_cluster_foldx_scaled_profiles) %>%
  t()

wt_aas <- foldx_scaled_pos_avg_mat$wt

# Assign cluster limited to those for wt AA
assign_cluster <- function(x, wt, dists){
  d <- dists[x, startsWith(colnames(dists), wt[x])]
  return(d[which.min(d)])
}
wt_cluster_distance <- sapply(1:length(wt_aas), assign_cluster, wt=wt_aas, dists=foldx_cluster_dists)
wt_cluster_assignment <- names(wt_cluster_distance)

# Assign cluster from all options
cluster_assignment <- colnames(foldx_cluster_dists)[apply(foldx_cluster_dists, 1, which.min)]
cluster_distance <- foldx_cluster_dists[1:nrow(foldx_cluster_dists) + (apply(foldx_cluster_dists, 1, which.min)-1)*nrow(foldx_cluster_dists)]
message('Proportion of time wt cluster == cluster: ',
        sum(cluster_assignment == wt_cluster_assignment) / length(wt_cluster_assignment))

foldx_scaled_pos_avg_clustered <- bind_cols(foldx_scaled_pos_avg_mat, 
                                            as_tibble(foldx_cluster_dists)) %>%
  mutate(wt_cluster = wt_cluster_assignment,
         wt_cluster_distance = wt_cluster_distance,
         free_cluster = cluster_assignment,
         free_cluster_distance = cluster_distance)

foldx_dm_cluster_avg <- group_by(foldx_scaled_pos_avg_clustered, wt_cluster) %>%
  summarise_at(vars(backbone_clash:van_der_waals_clashes), mean) %>%
  gather(key = 'term', value = 'foldx_ddG', -wt_cluster) %>%
  rename(cluster=wt_cluster) %>%
  full_join(.,
            t(dm_cluster_foldx_scaled_profiles) %>%
              as_tibble(rownames='cluster') %>%
              gather(key = 'term', value = 'dm_ddG', -cluster),
            by = c('cluster', 'term'))

plots$dm_cluster_assignment$foldx_dm_cluster_term_cor <- labeled_ggplot(
  p=ggplot(foldx_dm_cluster_avg, aes(x=dm_ddG, y=foldx_ddG, label=cluster)) +
    facet_wrap(~term, scales = 'free') +
    geom_point(colour='black') +
    geom_abline(slope = 1, linetype='dotted', colour='red') +
    geom_smooth(method = 'lm', colour='black') +
    labs(x='Deep Mut Cluster Mean ddG', y='FoldX Cluster Matches Mean ddG',
         title = 'Comparison of average FoldX profiles for deep mutagenesis derived clusters between the originals and proteome wide results'),
  units='cm', height=20, width=20)

plots$dm_cluster_assignment$foldx_dm_cluster_term_rel <- labeled_ggplot(
  p=ggplot(foldx_dm_cluster_avg, aes(x=term, y=cluster, fill=foldx_ddG - dm_ddG)) +
    geom_raster() +
    scale_fill_gradient2() +
    labs(x='', y='', title = 'Difference between proteome derived ddG and Deep Mut ddG for deep mut clusters') +
    coord_fixed() +
    theme(axis.ticks = element_blank(),
          panel.background = element_blank(),
          axis.title = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          axis.text.y = element_text(colour = AA_COLOURS[str_sub(sort(unique(foldx_dm_cluster_avg$cluster)), end = 1)]),
          plot.title = element_text(hjust = 0.5)),
  units='cm', height=length(unique(foldx_dm_cluster_avg$cluster))*0.5, width=20)

########

#### Cluster deep mut results based on FoldX terms ####
dm_foldx_scaled <- select(deep_mut_hclust$foldx, gene_name, position=pos, wt, total_energy:entropy_complex) %>%
  select(-total_energy, -sloop_entropy, -mloop_entropy, -electrostatic_kon, -entropy_complex, -water_bridge) %>%
  distinct() %>%
  bind_cols(select(., gene_name, position, wt),
            tibble_to_matrix(., backbone_hbond:energy_ionisation) %>%
              scale(center = FALSE, scale=term_scaling_factors) %>%
              as_tibble())

dm_foldx_settings <- list(h=8, k=NULL, max_k=6, min_k=3)
dm_foldx_hclust <- group_by(dm_foldx_scaled, wt) %>%
  do(hclust = make_hclust_clusters(., backbone_hbond:energy_ionisation, conf = dm_foldx_settings))
#sapply(dm_foldx_hclust$hclust, function(x){plot(x$hclust); abline(h = dm_foldx_settings$h)})

dm_foldx_hclust_tbl <- map_dfr(dm_foldx_hclust$hclust, .f = ~ .[[1]]) %>%
  mutate(cluster = str_c(wt, '_', cluster))

dm_foldx_hclust_analysis <- analyse_clusters(dm_foldx_hclust_tbl, make_hclust_cluster_str(dm_foldx_settings), 'mean FoldX terms',
                                            backbone_hbond:energy_ionisation,
                                            transform_ddg = function(x){x / max(abs(x))})

plots$dm_foldx_clusters <- dm_foldx_hclust_analysis$plots

# Compare to Deep Mut ER based clusters
dm_foldx_profiles <- select(dm_foldx_hclust_tbl, cluster, backbone_hbond:energy_ionisation) %>%
  group_by(cluster) %>%
  summarise_all(mean) %>%
  tibble_to_matrix(-cluster, row_names = 'cluster') %>%
  t() %>%
  set_colnames(str_c('fx_', colnames(.)))

# Correlation
dm_foldx_cluster_cor <- cor(dm_foldx_profiles, dm_cluster_foldx_scaled_profiles %>% set_colnames(str_c('er_', colnames(.)))) %>%
  as_tibble(rownames = 'fx_cluster') %>%
  gather(key = 'er_cluster', value = 'cor', -fx_cluster) %>%
  add_factor_order(er_cluster, fx_cluster, cor, sym = FALSE) %>%
  mutate(er_aa = str_sub(er_cluster, start=4, end=4),
         fx_aa = str_sub(fx_cluster, start=4, end=4))

plots$dm_foldx_clusters$er_cluster_cor <- labeled_ggplot(
  p = ggplot(dm_foldx_cluster_cor, aes(x=er_cluster, y=fx_cluster, fill=cor)) +
  geom_raster() +
  scale_fill_gradient2() +
  labs(x='', y='', title = 'Correlation between mean FoldX profiles of clusters derived from ER and FoldX scores') +
  coord_fixed() +
  theme(axis.ticks = element_blank(),
        panel.background = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(colour = AA_COLOURS[str_sub(levels(dm_foldx_cluster_cor$er_cluster), start=4, end=4)],
                                   angle = 90, vjust = 0.5),
        axis.text.y = element_text(colour = AA_COLOURS[str_sub(levels(dm_foldx_cluster_cor$fx_cluster), start=4, end=4)]),
        plot.title = element_text(hjust = 0.5)),
  units='cm', height=length(levels(dm_foldx_cluster_cor$fx_cluster)) * 0.5, width=length(levels(dm_foldx_cluster_cor$er_cluster)) * 0.5)

plots$dm_foldx_clusters$er_cluster_cor_per_aa <- labeled_ggplot(
  p = filter(dm_foldx_cluster_cor, fx_aa == er_aa) %>%
    ggplot(aes(x=er_cluster, y=fx_cluster, fill=cor)) +
    facet_wrap(~er_aa, scales = 'free') +
    geom_raster() +
    scale_fill_gradient2() +
    labs(x='', y='', title = 'Correlation between mean FoldX profiles of clusters derived from ER and FoldX scores') +
    theme(axis.ticks = element_blank(),
          panel.background = element_blank(),
          axis.title = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5),
          plot.title = element_text(hjust = 0.5)),
  units='cm', height=length(levels(dm_foldx_cluster_cor$fx_cluster)) * 0.5, width=length(levels(dm_foldx_cluster_cor$er_cluster)) * 0.5)

plots$dm_foldx_clusters$er_cluster_cor_vs_size <- filter(dm_foldx_cluster_cor, fx_aa == er_aa) %>%
  mutate_at(vars(fx_cluster, er_cluster), as.character) %>%
  left_join(select(dm_foldx_hclust_analysis$sizes, fx_cluster=cluster, FoldX=n) %>% 
              mutate(fx_cluster = str_c('fx_', fx_cluster)),
            by = "fx_cluster") %>%
  left_join(select(deep_mut_hclust$cluster_sizes, er_cluster=cluster, ER=n) %>% 
              mutate(er_cluster = str_c('er_', er_cluster)),
            by = 'er_cluster') %>%
  gather(key = 'type', value = 'size', FoldX, ER) %>%
  ggplot(aes(x=cor, y=size, colour=er_aa)) +
  facet_wrap(~type) +
  geom_point() +
  scale_colour_manual(values = AA_COLOURS)

# Distance
dm_foldx_cluster_dist <- cbind(dm_foldx_profiles,
                               dm_cluster_foldx_scaled_profiles %>% set_colnames(str_c('er_', colnames(.)))) %>%
  t() %>%
  dist() %>%
  as.matrix() %>%
  as_tibble(rownames = 'fx_cluster') %>%
  gather(key = 'er_cluster', value = 'dist', -fx_cluster) %>%
  filter(startsWith(fx_cluster, 'fx'), startsWith(er_cluster, 'er')) %>%
  add_factor_order(er_cluster, fx_cluster, dist, sym = FALSE) %>%
  mutate(er_aa = str_sub(er_cluster, start=4, end=4),
         fx_aa = str_sub(fx_cluster, start=4, end=4))

plots$dm_foldx_clusters$er_cluster_dist <- labeled_ggplot(
  p = ggplot(dm_foldx_cluster_dist, aes(x=er_cluster, y=fx_cluster, fill=atan(0.5*dist))) +
    geom_raster() +
    scale_fill_viridis_c() +
    labs(x='', y='', title = 'Distance between mean FoldX profiles of clusters derived from ER and FoldX scores') +
    coord_fixed() +
    theme(axis.ticks = element_blank(),
          panel.background = element_blank(),
          axis.title = element_blank(),
          axis.text.x = element_text(colour = AA_COLOURS[str_sub(levels(dm_foldx_cluster_dist$er_cluster), start=4, end=4)],
                                     angle = 90, vjust = 0.5),
          axis.text.y = element_text(colour = AA_COLOURS[str_sub(levels(dm_foldx_cluster_dist$fx_cluster), start=4, end=4)]),
          plot.title = element_text(hjust = 0.5)),
  units='cm', height=length(levels(dm_foldx_cluster_dist$fx_cluster)) * 0.5, width=length(levels(dm_foldx_cluster_dist$er_cluster)) * 0.5)

plots$dm_foldx_clusters$er_cluster_dist_per_aa <- labeled_ggplot(
  p = filter(dm_foldx_cluster_dist, fx_aa == er_aa) %>%
    ggplot(aes(x=er_cluster, y=fx_cluster, fill=atan(0.5*dist))) +
    facet_wrap(~er_aa, scales = 'free') +
    geom_raster() +
    scale_fill_viridis_c() +
    labs(x='', y='', title = 'Distance between mean FoldX profiles of clusters derived from ER and FoldX scores') +
    theme(axis.ticks = element_blank(),
          panel.background = element_blank(),
          axis.title = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5),
          plot.title = element_text(hjust = 0.5)),
  units='cm', height=length(levels(dm_foldx_cluster_dist$fx_cluster)) * 0.5, width=length(levels(dm_foldx_cluster_dist$er_cluster)) * 0.5)
########

# Save Plots
save_plot_list(plots, root='figures/8_foldx_score_analysis/')
