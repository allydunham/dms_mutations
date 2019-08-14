#!/usr/bin/env Rscript 
# PFunctions for analysing FoldX scores clustering

# analyse a dataframe of clusters
# ... columns on which clustering was run
analyse_clusters <- function(tbl, str, origin_str, ..., transform_ddg=function(x){x}){
  cols <- enquos(...)
  plots <- list()
  
  profiles <- group_by(tbl, cluster) %>%
    summarise_at(.vars = vars(!!!cols), .funs = mean)
  
  profiles_long <- gather(profiles, key = 'term', value = 'ddG', -cluster) %>%
    add_factor_order(cluster, term, ddG, sym = FALSE) %>%
    group_by(term) %>%
    mutate(trans_ddg = transform_ddg(ddG)) %>%
    ungroup()
  
  sizes <- group_by(tbl, cluster) %>%
    summarise(n = n()) %>%
    mutate(aa = str_sub(cluster, end = 1), cluster = factor(cluster, levels = levels(profiles_long$cluster)))
  
  cors <- transpose_tibble(profiles, cluster, name_col = 'term') %>%
    tibble_correlation(-term) %>%
    rename(cluster1 = cat1, cluster2 = cat2) %>%
    mutate(wt1 = str_sub(cluster1, end = 1),
           wt2 = str_sub(cluster2, end = 1)) %>%
    mutate(pair = mapply(function(x, y){str_c(str_sort(c(x, y)), collapse = '')}, wt1, wt2))

  plots$sizes <- ggplot(sizes, aes(x=cluster, y=n, fill=aa)) +
    geom_col() +
    scale_fill_manual(values = AA_COLOURS) +
    labs(x='Cluster', y='Size', title = str_c('Size of ', str,  ' clusters from ', origin_str)) +
    theme_pubclean() +
    scale_y_log10() +
    guides(fill=FALSE) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
  plots$mean_profile <- labeled_ggplot(
    p = ggplot(profiles_long, aes(x=term, y=cluster, fill=trans_ddg)) +
      geom_raster() +
      scale_fill_gradient2() +
      coord_fixed() +
      theme(axis.ticks = element_blank(),
            panel.background = element_blank(),
            axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5),
            axis.text.y = element_text(colour = AA_COLOURS[str_sub(levels(profiles_long$cluster), end = 1)]),
            plot.title = element_text(hjust = 0.5)) +
      labs(y='Cluster', x='FoldX Term', title = str_c('Mean profile of ', str, ' clusters ', origin_str)),
    units='cm', height = length(levels(profiles_long$cluster)) * 0.5, width = 15)
  
  
  plots$correlation <- labeled_ggplot(
    p = ggplot(cors, aes(x=cluster1, y=cluster2, fill=cor)) +
      geom_raster() +
      scale_fill_gradient2() +
      labs(x='', y='', title = str_c('Correlation of ', str, ' cluster centroids based on ', origin_str)) +
      coord_fixed() +
      theme(axis.ticks = element_blank(),
            panel.background = element_blank(),
            axis.title = element_blank(),
            axis.text.x = element_text(colour = AA_COLOURS[str_sub(levels(cors$cluster1), end = 1)],
                                       angle = 90, vjust = 0.5),
            axis.text.y = element_text(colour = AA_COLOURS[str_sub(levels(cors$cluster2), end = 1)]),
            plot.title = element_text(hjust = 0.5)),
    units='cm',
    height = length(levels(cors$cluster1)) * 0.5,
    width = length(levels(cors$cluster2)) * 0.5)
  
  return(list(sizes=sizes,
              mean_profiles=profiles,
              mean_profiles_long=profiles_long,
              profile_correlations=cors,
              profile_order = levels(profiles_long$cluster),
              correlation_order = levels(cors$cluster1),
              plots = plots))
}
