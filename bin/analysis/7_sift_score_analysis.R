#!/usr/bin/env Rscript 
# Perform complimentary analysis on SIFT scores as DMS data

source('src/config.R')
source('src/analysis/sift_score_analysis.R')
source('src/analysis/position_profile_clustering.R')

sift <- readRDS('data/rdata/human_sift_reduced_log10.RDS') %>%
  rename(uniprot_id = acc)
foldx <- readRDS('data/rdata/human_foldx_reduced.RDS') %>%
  rename(position = pos)

sum_missing <- tibble_to_matrix(sift, A:Y) %>% is.na() %>% sum()
total_scores <- dim(sift)[1] * 20
prop_missing <- sum_missing / total_scores

message(str_c(signif(prop_missing, 4), '% of SIFT scores missing (', sum_missing, ' / ', total_scores,')'))

plots <- list()

#### Summary statistics for sift ####
plots$summary <- plot_sift_score_summary(sift)

########

#### Cluster AA profiles ####
# kmeans
n <- 4
kmeans_clusters <- group_by(sift, wt) %>%
  do(hclust = make_kmeans_clusters(., A:Y, n = n))

kmeans_tbl <- map_dfr(kmeans_clusters$hclust, .f = ~ .[[1]]) %>%
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

# Save plots
save_plot_list(plots, root='figures/7_sift_score_analysis')
