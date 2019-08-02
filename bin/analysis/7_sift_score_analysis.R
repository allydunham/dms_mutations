#!/usr/bin/env Rscript 
# Perform complimentary analysis on SIFT scores as DMS data

source('src/config.R')
source('src/analysis/sift_score_analysis.R')
source('src/analysis/position_profile_clustering.R')

# Misc data import
deep_mut_hclust <- readRDS('data/rdata/deep_mut_hclust_clusters.RDS')

foldx <- readRDS('data/rdata/human_foldx_reduced.RDS')

# 100 Sampled proteins sift scores
sift <- readRDS('data/rdata/human_sift_reduced.RDS') %>%
  import_sift_scores()

# 50000 Normal/PTM sites sift scores
sift_ptms <- readRDS('data/rdata/human_sift_ptms.RDS') %>%
  import_sift_scores()

# Perform general analyses
# sapply(hclust_clusters$hclust, function(x){plot(x$hclust)})
if (file.exists('data/rdata/sift_cluster_analysis.RDS')){
  analyses <- readRDS('data/rdata/sift_cluster_analysis.RDS')
} else {
  analyses <- list()
}

if (FALSE){
  # Only redo analyses if specifically asked for
  analyses$sift <- analyse_sift_clusters(sift$standard, select(foldx, -chain), deep_mut_hclust, score_str='SIFT', n=4, h=8)
  
  analyses$log10_sift <- analyse_sift_clusters(sift$log10, select(foldx, -chain), deep_mut_hclust,
                                               score_str='log10(SIFT + e)', n=4, h=25)
  
  analyses$log10_sift_filtered <- filter(sift$log10, !sift$invariant_positions, !sift$permissive_positions) %>%
    analyse_sift_clusters(select(foldx, -chain), deep_mut_hclust, score_str='log10(SIFT + e)', n=4, h=25)
  
  analyses$wt_norm <- analyse_sift_clusters(sift$wt_norm, select(foldx, -chain), deep_mut_hclust,
                                                     score_str='log10(SIFT/wt + e)', n=4, h=30)
  
  analyses$wt_norm_filtered <- filter(sift$wt_norm, !sift$invariant_positions, !sift$permissive_positions) %>%
    analyse_sift_clusters(select(foldx, -chain), deep_mut_hclust, score_str='log10(SIFT/wt + e)', n=4, h=25)
  
  analyses$mean_norm_filtered <- filter(sift$mean_norm, !sift$invariant_positions, !sift$permissive_positions) %>%
    analyse_sift_clusters(select(foldx, -chain), deep_mut_hclust, score_str='log2(SIFT/mean sub + e)', n=4, h=100)
  
  # PTM set
  analyses$ptms_sift <- analyse_sift_clusters(sift_ptms$standard, select(foldx, -chain), deep_mut_hclust, score_str='SIFT', n=4, h=8)
  
  analyses$ptms_log10_sift <- analyse_sift_clusters(sift_ptms$log10, select(foldx, -chain), deep_mut_hclust,
                                               score_str='log10(SIFT + e)', n=4, h=25)
  
  analyses$ptms_log10_sift_filtered <- filter(sift_ptms$log10, !sift_ptms$invariant_positions, !sift_ptms$permissive_positions) %>%
    analyse_sift_clusters(select(foldx, -chain), deep_mut_hclust, score_str='log10(SIFT + e)', n=4, h=25)
  
  analyses$ptms_wt_norm <- analyse_sift_clusters(sift_ptms$wt_norm, select(foldx, -chain), deep_mut_hclust,
                                                     score_str='log10(SIFT/wt + e)', n=4, h=30)
  
  analyses$ptms_wt_norm_filtered <- filter(sift_ptms$wt_norm, !sift_ptms$invariant_positions, !sift_ptms$permissive_positions) %>%
    analyse_sift_clusters(select(foldx, -chain), deep_mut_hclust, score_str='log10(SIFT + e)', n=4, h=25)
  
  analyses$ptms_mean_norm_filtered <- filter(sift_ptms$mean_norm, !sift_ptms$invariant_positions, !sift_ptms$permissive_positions) %>%
    analyse_sift_clusters(select(foldx, -chain), deep_mut_hclust, score_str='log2(SIFT/mean sub + e)', n=4, h=100)
}

analyses$ptms_log10_sift$plots$hclust$cluster_modifications <- ggplot(analyses$ptms_log10_sift$hclust$tbl, aes(x=cluster, fill=wt)) +
  facet_grid(rows = vars(modification), cols = vars(wt), scales = 'free') +
  geom_bar(aes(y = ..prop.., group = wt)) +
  scale_fill_manual(values = AA_COLOURS) +
  guides(fill = FALSE) +
  theme(axis.ticks = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5))

# Save plots and analyses
write_rds(analyses, 'data/rdata/sift_cluster_analysis.RDS')
for (n in names(analyses)){
  save_plot_list(analyses[[n]]$plots, root=str_c('figures/7_sift_score_analysis/', n))
}
