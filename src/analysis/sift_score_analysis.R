#!/usr/bin/env Rscript 
# Functions for analysis of SIFT Scores in the manner of dms data

#### Overall analysis ####
# import sift scores in standard/log10/wt normalised forms, plus masks for invariant/permisive positions 
import_sift_scores <- function(tbl){
  invariant <- tibble_to_matrix(tbl, A:Y) == 0
  invariant <- rowSums(invariant, na.rm = TRUE) == 19
  
  permissive <- tibble_to_matrix(tbl, A:Y) < 0.05
  permissive <- rowSums(permissive) == 0
  
  min_non_zero <- tibble_to_matrix(tbl, A:Y)
  min_non_zero <- min(min_non_zero[min_non_zero > 0], na.rm = TRUE)
  
  tbl_log10 <- mutate_at(tbl, .vars = vars(A:Y), .funs = ~ log10(. + min_non_zero))
  
  # Normalise by wt score at position
  mat <- tibble_to_matrix(tbl, A:Y)
  wt <- sapply(1:nrow(mat), function(x){mat[x, tbl$wt[x]]}) %>% unname()
  mat_norm <- mat / wt
  
  min_non_zero_wt_norm <- min(mat_norm[mat_norm > 0])
  wt_norm <- select(tbl, -(A:Y)) %>%
    bind_cols(., as_tibble(mat_norm)) %>%
    mutate_at(vars(A:Y), .funs = ~ log10(. + min_non_zero_wt_norm))
  
  # Normalise by average profile of wt AA
  mean_profiles <- filter(tbl, !permissive, !invariant) %>%
    group_by(wt) %>%
    summarise_at(vars(A:Y), mean) %>%
    tibble_to_matrix(A:Y, row_names = 'wt')
  
  mat_mean_norm <- mat / mean_profiles[tbl$wt,]
  min_non_zero_mean_norm <- min(mat_mean_norm[mat_mean_norm > 0])
  
  mean_norm <- select(tbl, -(A:Y)) %>%
    bind_cols(., as_tibble(mat_mean_norm)) %>%
    mutate_at(vars(A:Y), .funs = ~ log2(. + min_non_zero_mean_norm))
  
  return(list(standard=tbl, log10=tbl_log10, wt_norm=wt_norm, mean_norm=mean_norm,
              invariant_positions=invariant, permissive_positions=permissive))
}

analyse_sift_clusters <- function(sift, foldx, deep_mut_hclust, score_str='SIFT', n=4, h=NULL, k=NULL){
  if (is.null(h) & is.null(k)){stop('One of h or k must be specified')}
  
  sum_missing <- tibble_to_matrix(sift, A:Y) %>% is.na() %>% sum()
  total_scores <- dim(sift)[1] * 20
  prop_missing <- sum_missing / total_scores
  
  message(str_c(signif(prop_missing, 4), '% of SIFT scores missing (', sum_missing, ' / ', total_scores,')'))
  
  plots <- list()
  
  ## Summary statistics for sift
  plots$summary <- plot_sift_score_summary(sift, score_labeler = function(x){x})
  
  ## Cluster AA profiles
  # kmeans
  kmeans_clusters <- group_by(sift, wt) %>%
    do(kmeans = make_kmeans_clusters(., A:Y, n = n))
  
  kmeans_tbl <- map_dfr(kmeans_clusters$kmeans, .f = ~ .[[1]]) %>%
    mutate(cluster = str_c(wt, '_', cluster))
  
  kmeans_analysis <- cluster_analysis(kmeans_tbl, cluster_str = str_c('kmeans (n = ', n, ')'), er_str = score_str,
                                      foldx = foldx, id_col = uniprot_id)
  
  plots$kmeans <- kmeans_analysis$plots
  plots$kmeans <- plots$kmeans[!sapply(plots$kmeans, is.null)]
  
  plots$kmeans$cluster_metrics <- plot_cluster_metrics(kmeans_tbl, cluster_order = kmeans_analysis$cluster_cor_order,
                                                       median_ic, n_aa, n_seq) %>%
    labeled_ggplot(units = 'cm', height = 30, width = length(unique(kmeans_tbl$cluster)) * 0.5 + 2)
  
  # Hclust
  hclust_clusters <- group_by(sift, wt) %>%
    do(hclust = make_hclust_clusters(., A:Y, h = h))
  
  hclust_tbl <- map_dfr(hclust_clusters$hclust, .f = ~ .[[1]]) %>%
    mutate(cluster = str_c(wt, '_', cluster))
  
  hclust_analysis <- cluster_analysis(hclust_tbl, cluster_str = str_c('hclust (h = ', h, ')'), er_str = score_str,
                                      foldx = foldx, id_col = uniprot_id)
  
  plots$hclust <- hclust_analysis$plots
  plots$hclust <- plots$hclust[!sapply(plots$hclust, is.null)]
  
  plots$hclust$cluster_metrics <- plot_cluster_metrics(hclust_tbl, cluster_order = hclust_analysis$cluster_cor_order,
                                                       median_ic, n_aa, n_seq) %>%
    labeled_ggplot(units = 'cm', height = 30, width = length(unique(hclust_tbl$cluster)) * 0.5 + 2)
  
  ## Compare to deep mut clusters
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
  
  plots$hclust$deep_mut_cluster_comparison_per_aa <- labeled_ggplot(
    p = ggplot(filter(cluster_cors, str_sub(sift_cluster, end = 1) == str_sub(dm_cluster, end = 1)) %>%
                 mutate(aa = str_sub(sift_cluster, end = 1)), 
               aes(x=sift_cluster, y=dm_cluster, fill=cor)) +
      facet_wrap(~aa, nrow = 4, scales = 'free') +
      geom_tile() +
      scale_fill_gradient2() +
      ggtitle('Correlation between cluster profiles determined from SIFT and deep mutational scanning scores') +
      guides(fill=guide_colourbar(title = 'Spearmans Rank Correlation')) +
      xlab('SIFT') +
      ylab('Deep Mutational Scanning') +
      theme(axis.ticks = element_blank(),
            strip.background = element_blank(),
            panel.background = element_blank(),
            axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)),
    units = 'cm', width = 30, height = 20)
  
  top_cors <- arrange(cluster_cors, desc(cor)) %>%
    group_by(dm_cluster) %>%
    summarise(first = nth(sift_cluster, n = 1), second = nth(sift_cluster, n = 2),
              third = nth(sift_cluster, n = 3), fourth = nth(sift_cluster, n = 4),
              fifth = nth(sift_cluster, n = 5), bot_first = nth(sift_cluster, n = -1),
              bot_second = nth(sift_cluster, n = -2), bot_third = nth(sift_cluster, n = -3),
              bot_fourth = nth(sift_cluster, n = -4), bot_fifth = nth(sift_cluster, n = -5)) %>%
    arrange(as.character(dm_cluster))
  
  return(list(hclust=list(clusters=hclust_clusters, tbl=hclust_tbl, analysis=hclust_analysis,
                          dm_cluster_cors=cluster_cors, dm_top_cors=top_cors),
              kmeans=list(clusters=kmeans_clusters, tbl=kmeans_tbl, analysis=kmeans_analysis),
              plots=plots))
}
########

#### Summary distributions ####
plot_sift_score_summary <- function(tbl, score_labeler=NULL){
  if (is.null(score_labeler)){
    score_labeler <- make_log_labeler(base=10, force = 'exp')
  }
  
  plots <- list()
  plots$score_distribution <- ggplot(gather(tbl, key = 'aa', value = 'score', A:Y),
                                             aes(x=score)) +
    facet_wrap(~aa) +
    geom_histogram() +
    scale_x_continuous(labels = score_labeler) +
    scale_y_log10()
  
  plots$per_sub_type_score_distribution <- labeled_ggplot(
    p = ggplot(gather(tbl, key = 'mut', value = 'score', A:Y), aes(x=score)) +
      geom_histogram() +
      facet_grid(rows = vars(mut), cols = vars(wt)) + 
      ggtitle(expression('Distribution of log'[10]*'(SIFT +'~epsilon*') for each substitution (wt: cols, mut: rows)')) +
      scale_x_continuous(labels = score_labeler) +
      scale_y_log10() +
      theme_pubclean() +
      theme(strip.background = element_blank(),
            legend.position = 'right'),
    units = 'cm', height = 44, width = 44)
  
  plots$protein_length <- group_by(tbl, uniprot_id) %>%
    summarise(length = max(position)) %>%
    ggplot(aes(x=length)) +
    geom_histogram()
  
  plots$aa_sub_distribution <- gather(tbl, key = 'mut', value = 'score', A:Y) %>%
    select(uniprot_id, position, wt, mut, score) %>%
    gather(key = 'type', value = 'aa', wt, mut) %>%
    drop_na(score) %>%
    mutate(aa = factor(aa, levels = names(AA_COLOURS))) %>%
    ggplot(aes(x=aa, fill=aa)) +
    geom_bar() +
    scale_fill_manual(values = AA_COLOURS) +
    facet_wrap(~type)
  
  plots$metric_distribution <- ggpairs(tbl, columns = c('median_ic', 'n_aa', 'n_seq'))
  
  return(plots)
}
########

#### Cluster Plots ####
# ... = metric columns to plot
plot_cluster_metrics <- function(tbl, cluster_order=NULL, ...){
  metric_cols <- enquos(...)
  
  if (!is.null(cluster_order)){
    tbl <- mutate(tbl, cluster = factor(cluster, levels = cluster_order))
  }
  
  cluster_metrics <- gather(tbl, 'metric', 'value', !!! metric_cols)
  
  p <- ggplot(cluster_metrics, aes(x=cluster, y=value, fill=wt)) +
      facet_wrap(~metric, ncol = 1, scales = 'free_y', strip.position = 'left') +
      ylab(NULL) +
      geom_boxplot() +
      scale_fill_manual(values = AA_COLOURS) +
      guides(fill=FALSE) +
      theme_pubclean() +
      theme(axis.text.x = element_text(colour = AA_COLOURS[str_sub(levels(cluster_metrics$cluster), end = 1)],
                                       angle = 90, hjust = 1, vjust = 0.5),
            strip.background = element_blank(),
            strip.placement = 'outside')
}
########