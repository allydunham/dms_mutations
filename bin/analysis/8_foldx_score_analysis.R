#!/usr/bin/env Rscript 
# Perform complimentary analysis on FoldX scores to that done on DMS data

source('src/config.R')
source('src/analysis/foldx_score_analysis.R')
source('src/analysis/position_profile_clustering.R')

foldx <- readRDS('data/rdata/human_foldx_tiny.RDS')

plots <- list()

#### Analyse by average at position ####
plots$pos_avg <- list()

# foldx_pos_avg <- group_by(foldx, uniprot_id, position, wt, modification) %>%
#   summarise_at(vars(total_energy:entropy_complex), mean, na.rm=TRUE) %>%
#   ungroup()
# saveRDS(foldx_pos_avg, 'data/rdata/human_foldx_tiny_pos_avg.RDS')
foldx_pos_avg <- readRDS('data/rdata/human_foldx_tiny_pos_avg.RDS')

# drop unused factors then center and scale (pssobly add water_bridge/kon back if structures where they're meaningful are used later)
foldx_pos_avg_scaled <- select(foldx_pos_avg, -sloop_entropy, -mloop_entropy, -entropy_complex, -electrostatic_kon, -water_bridge) %>% 
  tibble_to_matrix(., total_energy:energy_ionisation) %>% 
  scale(center = FALSE, scale = TRUE) %>% 
  as_tibble() %>%
  bind_cols(select(foldx_pos_avg, uniprot_id:modification), .)

h <- 2.5
k <- NULL
max_k <- 6
cluster_str <- str_c('hclust (h = ', h, ', k = ', k, ', max_k = ', max_k, ')')
hclust_clusters <- group_by(foldx_pos_avg_scaled, wt) %>%
  do(hclust = make_hclust_clusters(., total_energy:energy_ionisation, h = h, max_k = max_k))
#sapply(hclust_clusters$hclust, function(x){plot(x$hclust); abline(h = 12)})

hclust_cluster_tbl <- map_dfr(hclust_clusters$hclust, .f = ~ .[[1]]) %>%
  mutate(cluster = str_c(wt, '_', cluster))

hclust_cluster_profiles <- group_by(hclust_cluster_tbl, cluster) %>%
  summarise_at(.vars = vars(total_energy:energy_ionisation), .funs = mean)

hclust_cluster_profiles_long <- gather(hclust_cluster_profiles, key = 'term', value = 'ddG', -cluster) %>%
  add_factor_order(cluster, term, ddG, sym = FALSE) %>%
  group_by(term) %>%
  mutate(trans_ddg = ddG / max(abs(ddG))) %>%
  ungroup()

hclust_cluster_sizes <- group_by(hclust_cluster_tbl, cluster) %>%
  summarise(n = n()) %>%
  mutate(aa = str_sub(cluster, end = 1), cluster = factor(cluster, levels = levels(hclust_cluster_profiles_long$cluster)))

plots$pos_avg$cluster_sizes <- ggplot(hclust_cluster_sizes, aes(x=cluster, y=n, fill=aa)) +
  geom_col() +
  scale_fill_manual(values = AA_COLOURS) +
  labs(x='Cluster', y='Size', title = str_c('Size of ', cluster_str,  'clusters from mean position FoldX terms')) +
  theme_pubclean() +
  scale_y_log10() +
  guides(fill=FALSE) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

plots$pos_avg$cluster_mean_profile <- labeled_ggplot(
  p = ggplot(hclust_cluster_profiles_long, aes(x=term, y=cluster, fill=trans_ddg)) +
    geom_raster() +
    scale_fill_gradient2() +
    coord_fixed() +
    theme(axis.ticks = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5),
          axis.text.y = element_text(colour = AA_COLOURS[str_sub(levels(hclust_cluster_profiles_long$cluster), end = 1)]),
          plot.title = element_text(hjust = 0.5)) +
    labs(y='Cluster', x='FoldX Term', title = str_c('Mean profile of ', cluster_str, 'clusters from mean position FoldX terms')),
  units='cm', height = length(levels(hclust_cluster_profiles_long$cluster)) * 0.5, width = 15)

hclust_cluster_cors <- transpose_tibble(hclust_cluster_profiles, cluster, name_col = 'term') %>%
  tibble_correlation(-term) %>%
  rename(cluster1 = cat1, cluster2 = cat2) %>%
  mutate(wt1 = str_sub(cluster1, end = 1),
         wt2 = str_sub(cluster2, end = 1)) %>%
  mutate(pair = mapply(function(x, y){str_c(str_sort(c(x, y)), collapse = '')}, wt1, wt2))

plots$pos_avg$cluster_correlation <- ggplot(hclust_cluster_cors, aes(x=cluster1, y=cluster2, fill=cor)) +
  geom_raster() +
  scale_fill_gradient2() +
  labs(x='', y='', title = str_c('Correlation of ', cluster_str, ' cluster centroids based on mean FoldX profile')) +
  coord_fixed() +
  theme(axis.ticks = element_blank(),
        panel.background = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(colour = AA_COLOURS[str_sub(levels(hclust_cluster_cors$cluster1), end = 1)],
                                   angle = 90, vjust = 0.5),
        axis.text.y = element_text(colour = AA_COLOURS[str_sub(levels(hclust_cluster_cors$cluster2), end = 1)]),
        plot.title = element_text(hjust = 0.5))
plots$pos_avg$cluster_correlation <- labeled_ggplot(p = plots$pos_avg$cluster_correlation, units='cm',
                                                    height = length(levels(hclust_cluster_cors$cluster1)) * 0.5,
                                                    width = length(levels(hclust_cluster_cors$cluster2)) * 0.5)

########

# Save Plots
save_plot_list(plots, root='figures/8_foldx_score_analysis/')
