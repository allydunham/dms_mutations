#!/usr/bin/env Rscript 
# Perform complimentary analysis on SIFT scores as DMS data

source('src/config.R')
source('src/analysis/sift_score_analysis.R')
source('src/analysis/position_profile_clustering.R')

deep_mut_hclust <- readRDS('data/rdata/deep_mut_hclust_clusters.RDS')

foldx <- readRDS('data/rdata/human_foldx_reduced.RDS') %>%
  rename(position = pos)

sift_std <- readRDS('data/rdata/human_sift_reduced.RDS') %>%
  rename(uniprot_id = acc)

invariant_sift_position_inds <- tibble_to_matrix(sift_std, A:Y) == 0
invariant_sift_position_inds <- rowSums(invariant_sift_position_inds, na.rm = TRUE) == 19

permissive_sift_positon_inds <- tibble_to_matrix(sift_std, A:Y) < 0.05
permissive_sift_positon_inds <- rowSums(permissive_sift_positon_inds) == 0

sift_log10 <- readRDS('data/rdata/human_sift_reduced_log10.RDS') %>%
  rename(uniprot_id = acc)

# Normalised Sift Score, by wt
si_mat <- tibble_to_matrix(sift_std, A:Y)
wt_si <- sapply(1:nrow(si_mat), function(x){si_mat[x, sift_std$wt[x]]}) %>% unname()
si_mat_norm <- si_mat / wt_si

min_non_zero <- min(si_mat_norm[si_mat_norm > 0])
sift_wt_norm <- select(sift_std, uniprot_id, pos, wt, median_ic, n_aa, n_seq) %>%
  bind_cols(., as_tibble(si_mat_norm)) %>%
  mutate_at(vars(A:Y), .funs = ~ log10(. + min_non_zero))

# Select sift version to use?
sift_ver_name <- 'wt_norm'
sift <- sift_wt_norm

# filter invariant/permissive rows?
sift <- filter(sift, !invariant_sift_position_inds, !permissive_sift_positon_inds)

sum_missing <- tibble_to_matrix(sift, A:Y) %>% is.na() %>% sum()
total_scores <- dim(sift)[1] * 20
prop_missing <- sum_missing / total_scores

message(str_c(signif(prop_missing, 4), '% of SIFT scores missing (', sum_missing, ' / ', total_scores,')'))

plots <- list()
#### Summary statistics for sift ####
plots$summary <- plot_sift_score_summary(sift, score_labeler = function(x){x})

########

#### Cluster AA profiles ####
# kmeans
n <- 4
kmeans_clusters <- group_by(sift, wt) %>%
  do(kmeans = make_kmeans_clusters(., A:Y, n = n))

kmeans_tbl <- map_dfr(kmeans_clusters$kmeans, .f = ~ .[[1]]) %>%
  mutate(cluster = str_c(wt, '_', cluster))

kmeans_analysis <- cluster_analysis(kmeans_tbl, cluster_str = str_c('kmeans (n = ', n, ')'), er_str = 'log10(SIFT + e)',
                                    foldx = select(foldx, -chain), id_col = uniprot_id)

plots$kmeans <- kmeans_analysis$plots
plots$kmeans <- plots$kmeans[!sapply(plots$kmeans, is.null)]

# TODO move this to func
# Plot clusters againsts statistics
kmeans_cluster_metrics <- mutate(kmeans_tbl, cluster = factor(cluster, levels = kmeans_analysis$cluster_cor_order)) %>%
  gather('metric', 'value', median_ic, n_aa, n_seq)

plots$kmeans$cluster_metrics <- labeled_ggplot(
  p = ggplot(kmeans_cluster_metrics, aes(x=cluster, y=value, fill=wt)) +
    facet_wrap(~metric, ncol = 1, scales = 'free_y', strip.position = 'left') +
    ylab(NULL) +
    geom_boxplot() +
    scale_fill_manual(values = AA_COLOURS) +
    guides(fill=FALSE) +
    theme_pubclean() +
    theme(axis.text.x = element_text(colour = AA_COLOURS[str_sub(levels(kmeans_cluster_metrics$cluster), end = 1)],
                                     angle = 90, hjust = 1, vjust = 0.5),
          strip.background = element_blank(),
          strip.placement = 'outside'),
  units = 'cm', height = 30, width = length(unique(kmeans_cluster_metrics$cluster)) * 0.5 + 2)
  
# Hclust
h <- 30
hclust_clusters <- group_by(sift, wt) %>%
  do(hclust = make_hclust_clusters(., A:Y, h = h))
# sapply(hclust_clusters$hclust, function(x){plot(x$hclust); abline(h = h)})

hclust_tbl <- map_dfr(hclust_clusters$hclust, .f = ~ .[[1]]) %>%
  mutate(cluster = str_c(wt, '_', cluster))

hclust_analysis <- cluster_analysis(hclust_tbl, cluster_str = str_c('hclust (h = ', h, ')'), er_str = 'log10(SIFT + e)',
                                    foldx = select(foldx, -chain), id_col = uniprot_id)

plots$hclust <- hclust_analysis$plots
plots$hclust <- plots$hclust[!sapply(plots$hclust, is.null)]

# Plot clusters againsts statistics
hclust_cluster_metrics <- mutate(hclust_tbl, cluster = factor(cluster, levels = hclust_analysis$cluster_cor_order)) %>%
  gather('metric', 'value', median_ic, n_aa, n_seq)

plots$hclust$cluster_metrics <- labeled_ggplot(
  p = ggplot(hclust_cluster_metrics, aes(x=cluster, y=value, fill=wt)) +
    facet_wrap(~metric, ncol = 1, scales = 'free_y', strip.position = 'left') +
    ylab(NULL) +
    geom_boxplot() +
    scale_fill_manual(values = AA_COLOURS) +
    guides(fill=FALSE) +
    theme_pubclean() +
    theme(axis.text.x = element_text(colour = AA_COLOURS[str_sub(levels(hclust_cluster_metrics$cluster), end = 1)],
                                     angle = 90, hjust = 1, vjust = 0.5),
          strip.background = element_blank(),
          strip.placement = 'outside'),
  units = 'cm', height = 30, width = length(unique(hclust_cluster_metrics$cluster)) * 0.5 + 2)

########

#### Compare to deep mut clusters ####
cluster_cors <- cor(tibble_to_matrix(hclust_analysis$mean_profiles, A:Y, row_names = 'cluster') %>% t(),
                    tibble_to_matrix(deep_mut_hclust$mean_profiles, A:Y, row_names = 'cluster') %>% t(),
                    method = 'spearman') %>%
  as_tibble(rownames = 'sift_cluster') %>%
  gather(key = 'dm_cluster', value = 'cor', -sift_cluster) %>%
  add_factor_order(sift_cluster, dm_cluster, cor, sym = FALSE)

plots$hclust$deep_mut_cluster_comparison <- labeled_ggplot(
  p = ggplot(cluster_cors, aes(x=sift_cluster, y=dm_cluster, fill=cor)) +
    geom_tile() +
    scale_fill_gradient2() +
    coord_fixed() +
    ggtitle('Correlation between cluster profiles determined from SIFT and deep mutational scanning scores') +
    guides(fill=guide_colourbar(title = 'Spearmans Rank Correlation')) +
    xlab('SIFT') +
    ylab('Deep Mutational Scanning') +
    theme(axis.ticks = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(colour = AA_COLOURS[str_sub(levels(cluster_cors$sift_cluster), end = 1)],
                                     angle=90, hjust = 1, vjust = 0.5),
          axis.text.y = element_text(colour = AA_COLOURS[str_sub(levels(cluster_cors$dm_cluster), end = 1)])),
  units = 'cm', width = 0.5*length(unique(cluster_cors$sift_cluster)) + 4,
  height = 0.5*length(unique(cluster_cors$dm_cluster)) + 2, limitsize=FALSE)

top_cors <- arrange(cluster_cors, desc(cor)) %>%
  group_by(dm_cluster) %>%
  summarise(first = nth(sift_cluster, n = 1), second = nth(sift_cluster, n = 2),
            third = nth(sift_cluster, n = 3), fourth = nth(sift_cluster, n = 4),
            fifth = nth(sift_cluster, n = 5), bot_first = nth(sift_cluster, n = -1),
            bot_second = nth(sift_cluster, n = -2), bot_third = nth(sift_cluster, n = -3),
            bot_fourth = nth(sift_cluster, n = -4), bot_fifth = nth(sift_cluster, n = -5)) %>%
  arrange(as.character(dm_cluster))
########

# Save plots
save_plot_list(plots, root=str_c('figures/7_sift_score_analysis/', sift_ver_name))

# Save analyses
saveRDS(hclust_analysis, str_c('data/rdata/', sift_ver_name, '_hclust_clusters.RDS'))
