#!/usr/bin/env Rscript 
# Perform complimentary analysis on SIFT scores as DMS data

source('src/config.R')
source('src/analysis/sift_score_analysis.R')
source('src/analysis/position_profile_clustering.R')

sift <- readRDS('data/rdata/human_sift_reduced_log10.RDS')
sift_no_missing <- readRDS('data/rdata/human_sift_no_missing_reduced_log10.RDS')

sum_missing <- tibble_to_matrix(sift, A:Y) %>% is.na() %>% sum()
total_scores <- dim(sift)[1] * 20
prop_missing <- sum_missing / total_scores

message(str_c(signif(prop_missing, 4), '% of SIFT scores missing (', sum_missing, ' / ', total_scores,')'))

plots <- list()

#### Summary statistics for sift ####
plots$summary <- plot_sift_score_summary(sift)
p <- plot_sift_score_summary(sift_no_missing)
names(p) <- str_c(names(p), '_no_missing')
plots$summary <- c(plots$summary, p)

# Compare positions where information is missing or not
plots$summary$missing_comparison <- list()
comb <- bind_rows(yes=gather(sift, key = 'mut', value = 'score', A:Y),
                  no=gather(sift_no_missing, key = 'mut', value = 'score', A:Y),
                  .id = 'missing')

plots$summary$missing_comparison$per_sub_type_score_distribution <- labeled_ggplot(
  p = ggplot(comb, aes(x=score, fill=missing)) +
    geom_histogram(alpha = 0.5) +
    facet_grid(rows = vars(mut), cols = vars(wt)) + 
    scale_fill_manual(values = c(yes='red', no='black')) +
    scale_x_continuous(labels = make_log_labeler(base = 10, force = 'exp')) +
    scale_y_log10() +
    theme_pubclean() +
    theme(strip.background = element_blank(),
          legend.position = 'right'),
  units = 'cm', height = 44, width = 44)

plots$summary$missing_comparison$score_distribution <- labeled_ggplot(
  p = ggplot(gather(comb, key = 'type', value = 'aa', wt, mut) %>% mutate(type = factor(type, levels = c('wt', 'mut'))),
             aes(x=score, fill=missing)) +
    geom_histogram(alpha = 0.5) +
    facet_grid(rows = vars(type), cols = vars(aa)) + 
    scale_fill_manual(values = c(yes='red', no='black')) +
    scale_x_continuous(labels = make_log_labeler(base = 10, force = 'exp')) +
    scale_y_log10() +
    theme_pubclean() +
    theme(strip.background = element_blank(),
          legend.position = 'right',
          axis.text.x = element_text(angle = 45)),
  units = 'cm', height = 10, width = 44)

plots$summary$missing_comparison$metric_distribution <- labeled_ggplot(
  p = ggplot(gather(comb, key = 'metric', value = 'value', median_ic, n_aa, n_seq),
             aes(x=value, fill=missing)) +
    geom_histogram(alpha = 0.5) +
    facet_wrap(~metric, nrow = 3, scales = 'free') + 
    scale_fill_manual(values = c(yes='red', no='black')) +
    scale_y_log10() +
    theme_pubclean() +
    theme(strip.background = element_blank(),
          legend.position = 'right'),
  units = 'cm', height = 30, width = 15)

plots$summary$missing_comparison$aa_freq_distribution <- labeled_ggplot(
  p = ggplot(gather(comb, key = 'type', value = 'aa', wt, mut) %>%
               drop_na(score) %>%
               mutate(type = factor(type, levels = c('wt', 'mut'))),
             aes(x=aa, fill=missing, y=..prop.., group=missing)) +
    geom_bar(position = 'dodge') + 
    facet_wrap(~type, nrow = 2) +
    scale_fill_manual(values = c(yes='red', no='black')) +
    theme_pubclean() +
    theme(strip.background = element_blank(),
          legend.position = 'right'),
  units = 'cm', height = 20, width = 20)
########

#### Cluster AA profiles ####
# kmeans
n <- 4
kmeans_clusters <- group_by(sift_no_missing, wt) %>%
  do(hclust = make_kmeans_clusters(., A:Y, n = n))

kmeans_tbl <- map_dfr(kmeans_clusters$hclust, .f = ~ .[[1]]) %>%
  mutate(cluster = str_c(wt, '_', cluster))

kmeans_analysis <- cluster_analysis(kmeans_tbl, cluster_str = str_c('kmeans (n = ', n, ')'), er_str = 'log10(SIFT + e)')

plots$kmeans <- kmeans_analysis$plots
plots$kmeans <- plots$kmeans[!is.null(plots$kmeans)]

# TODO move this to func
# Plot clusters againsts statistics
kmeans_cluster_metrics <- mutate(kmeans_tbl, cluster = factor(cluster, levels = kmeans_analysis$cluster_mean_order)) %>%
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
h <- 18
hclust_clusters <- group_by(sift_no_missing, wt) %>%
  do(hclust = make_hclust_clusters(., A:Y, h = h))

hclust_tbl <- map_dfr(hclust_clusters$hclust, .f = ~ .[[1]]) %>%
  mutate(cluster = str_c(wt, '_', cluster))

hclust_analysis <- cluster_analysis(hclust_tbl, cluster_str = str_c('hclust (h = ', h, ')'), er_str = 'log10(SIFT + e)')

plots$hclust <- hclust_analysis$plots
plots$hclust <- plots$hclust[!is.null(plots$hclust)]

# Plot clusters againsts statistics
hclust_cluster_metrics <- mutate(hclust_tbl, cluster = factor(cluster, levels = hclust_analysis$cluster_mean_order)) %>%
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
