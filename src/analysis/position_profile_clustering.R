#!/usr/bin/env Rscript 
# Functions to perform clustering analysis on per position mutational profiles from deep mutagenesis studies

data("BLOSUM62")

#### PCA ####
# Generate PCA of mutational profiles
# TODO move to new tibble_pca func in misc_utils.R
positional_profile_PCA <- function(variant_matrix){
  pca <- prcomp(as.matrix(select(variant_matrix, A:Y)), center = TRUE, scale. = TRUE)
  pca_variants <- bind_cols(select(variant_matrix, -(A:Y)), as_tibble(pca$x))
  
  return(list(profiles=pca_variants, pca=pca))
}

basic_pca_plots <- function(pca){
  plots <- list()
  plots$all_pcs <- plot_all_pcs(pca$profiles, colour_var = 'wt')
  
  
  plots$by_authour <- ggplot(pca$profiles, aes(x=PC1, y=PC2, colour=gene_name)) + 
    facet_wrap(~study) +
    geom_point()
  
  plots$secondary_structure <- plot_all_pcs(pca$profiles, colour_var = 'ss')
  plots$secondary_structure_reduced <- plot_all_pcs(pca$profiles, colour_var = 'ss_reduced')
  
  plots$fields_group_studies <- ggplot(filter(pca$profiles,
                                              authour %in% c('Araya et al.', 'Melamed et al.', 'Starita et al.',
                                                             'Kitzman et al.', 'Weile et al.')),
                                       aes(x=PC1, y=PC2, colour=study)) +
    geom_point()
  
  plots$position_sig <- ggplot(pca$profiles, aes(x=PC1, y=PC2, colour=sig_count)) +
    geom_point()
  
  plots$by_aa <- ggplot(pca$profiles, aes(x=PC1, y=PC2, colour=gene_name)) + 
    facet_wrap(~wt) +
    geom_point()
  
  plots$surface_accesibility <- ggplot(pca$profiles, aes(x=PC1, y=PC2, colour=all_atom_rel)) +
    geom_point() +
    scale_colour_gradientn(colours = c('blue', 'green', 'yellow', 'orange', 'red'))
  
  return(plots)
}

get_avg_aa_pca_profile <- function(pca, aa_col='wt'){
  aa_col_sym <- sym(aa_col)
  avg_profile <- pca$profiles %>%
    group_by(!!aa_col_sym) %>%
    summarise_at(.vars = vars(starts_with('PC')), .funs = list(~ mean(.)))
  
  cor_mat <- select(avg_profile, -!!aa_col_sym) %>% 
    t() %>% 
    set_colnames(avg_profile[[aa_col]]) %>%
    cor()
  
  aa_order <- rownames(cor_mat)[hclust(dist(cor_mat))$order]
  
  cor_tbl <- cor_mat %>%
    as_tibble(rownames = 'AA1') %>%
    gather(key = 'AA2', value = 'cor', -AA1) %>%
    mutate(AA1 = factor(AA1, levels = aa_order),
           AA2 = factor(AA2, levels = aa_order))
  
  return(list(avg_profile=avg_profile, cor_mat=cor_mat, cor_tbl=cor_tbl, aa_order=aa_order))
}

plot_aa_pca_profile_average_cor <- function(pca){
  cors <- tibble_to_matrix(pca$profiles, PC1:PC20,
                           row_names = str_c(pca$profiles$study, pca$profiles$pos, pca$profiles$wt, sep = '~')) %>%
    t() %>%
    cor() %>%
    as_tibble(rownames = 'pos1') %>%
    gather(key = 'pos2', value = 'cor', -pos1) %>%
    mutate(AA1 = str_sub(pos1, start = -1),
           AA2 = str_sub(pos2, start = -1)) %>%
    group_by(AA1, AA2) %>%
    summarise(cor = mean(cor)) %>%
    ungroup()
  
  aa_order <- spread(cors, key = AA2, value = cor) %>%
    tibble_to_matrix(., A:Y, row_names = .$AA1)
  aa_order <- rownames(aa_order)[hclust(dist(aa_order))$order]
  
  cors <- mutate(cors, AA1 = factor(AA1, levels = aa_order), AA2 = factor(AA2, levels = aa_order))
  
  return(
    ggplot(cors, aes(x=AA1, y=AA2, fill=cor)) +
      geom_tile(colour='white') + 
      scale_fill_gradient2() +
      theme(axis.ticks = element_blank(), panel.background = element_blank())
  )
}

pca_factor_cor <- function(pca, .vars){
  pcas_mat <- select(pca$profiles, starts_with('PC')) %>% 
    as.matrix()
  
  factor_mat <- select(pca$profiles, !!!.vars) %>% 
    as.matrix()
  
  cor_mat <- cor(pcas_mat, factor_mat, use = 'pairwise.complete.obs')
  
  cor_tbl <- cor_mat %>%
    as_tibble(rownames = 'PC') %>%
    gather(key = 'factor', value = 'cor', -PC) %>%
    mutate(PC = factor(PC, levels = str_c('PC', 1:dim(cor_mat)[1])))
  
  return(list(tbl=cor_tbl, matrix=cor_mat))
}

pca_factor_heatmap <- function(pca){
  ggplot(pca$tbl, aes(x=PC, y=factor, fill=cor)) + 
    geom_tile(colour='white') + 
    scale_fill_gradient2() +
    theme(axis.ticks = element_blank(), panel.background = element_blank())
}

aa_avg_profile_plot <- function(x){list(avg_aa_profile=ggplot(x$avg_profile, aes(x=PC1, y=PC2, label=wt)) + geom_text())}

aa_profile_heatmap <- function(pca){list(
  aa_profile_heatmap=ggplot(pca$cor_tbl, aes(x=AA1, y=AA2, fill=cor)) + 
    geom_tile(colour='white') + 
    scale_fill_gradient2() +
    theme(axis.ticks = element_blank(), panel.background = element_blank())
)}

per_aa_pcas <- function(aa, variant_matrix){
  variant_matrix <- filter(variant_matrix, wt == aa)
  
  profile_pca <- positional_profile_PCA(variant_matrix)
  surface_cor <- pca_surf_acc_cor(profile_pca)
  
  basic_plots <- basic_pca_plots(profile_pca)
  surface_heatmap <- pca_surface_heatmap(surface_cor)
  
  return(c(basic_plots, pc_surface_acc_heatmap=surface_heatmap))
}
########

#### kmeans ####
make_kmeans_clusters <- function(tbl, cols, n=5, ...){
  cols <- enquo(cols)
  
  mat <- tibble_to_matrix(tbl, !!cols)
  
  km <- kmeans(mat, centers = n, ...)
  
  return(list(tbl=mutate(tbl, cluster = km$cluster),
              kmeans=km))
}
########

#### hclust ####
# Perfrom hclust on columns of a tibble, using parameters in conf or by specific h, k, ... settings if given. conf takes preference
# k overrides h as in the base hclust
make_hclust_clusters <- function(tbl, cols, dist_method = 'manhattan', conf=NULL, h = NULL, k = NULL, max_k=Inf, min_k=0, ...){
  cols <- enquo(cols)
  
  defaults <- list(h=h, k=k, max_k=max_k, min_k=min_k)
  if (is.null(conf)){
    conf <- defaults
  } else {
    conf <- list_modify(defaults, !!!conf)
  }
  
  mat <- tibble_to_matrix(tbl, !!cols)
  hc <- hclust(dist(mat, method = dist_method), ...)
  clus <- cutree(hc, k = conf$k, h = conf$h)
  
  # Use max/min cluster nums if using h (defaults mean any number is allowed)
  if (is.null(conf$k)){
    # too many clusters
    if (max(clus) > conf$max_k){
      clus <- cutree(hc, k = conf$max_k)
    }
    # too few clusters
    if (max(clus) < conf$min_k){
      clus <- cutree(hc, k = conf$min_k)
    }
  }
  
  return(list(tbl = mutate(tbl, cluster = clus),
              hclust = hc))
}

# Generate a sensible name for an hclust run passing a config list and/or individual values for the params (overrides settings)
make_hclust_cluster_str <- function(conf=NULL, ...){
  manual <- list(...)
  
  if (is.null(conf)){
    conf <- list(h=NULL, k=NULL, max_k=NULL, min_k=NULL)
  }
  
  if (length(manual) > 0){
    conf <- list_modify(conf, manual)
  }
  
  conf <- conf[!sapply(conf, is.null)]
  
  return(str_c('hclust ', str_sub(capture.output(dput(conf)), start = 5)))
}
########

#### Cluster analysis ####
# Expects a tbl with a columns:
# study - deep mut study
# pos - position in protein
# wt - wt AA at that position
# cluster - cluster assignment of the position

# backbone_angles = tbl giving psi/phi for each study/pdb_id/chain/aa/position combo
# foldx = tbl giving FoldX derived energy terms for deep mut positions

cluster_analysis <- function(tbl, backbone_angles=NULL, foldx=NULL, cluster_str='<UNKNOWN>', er_str='<UNKNOWN>',
                             id_col=NULL, pos_col=NULL){
  id_col <- enquo(id_col)
  if (rlang::quo_is_null(id_col)){
    id_col <- quo(study)
    id_col_str <- 'study'
  } else {
    id_col_str <- rlang::as_name(id_col)
  }
  
  pos_col <- enquo(pos_col)
  if (rlang::quo_is_null(pos_col)){
    pos_col <- quo(position)
    pos_col_str <- 'position'
  } else {
    pos_col_str <- rlang::as_name(pos_col)
  }
  
  # Ramachandran Plot
  if (!is.null(backbone_angles)){
    angles <- left_join(rename(backbone_angles, !!pos_col:=position, wt=aa),
                        select(tbl, study, !!pos_col, wt, cluster),
                        by = c('study', 'pos', 'wt')) %>%
      drop_na(cluster) %>%
    mutate(cluster_num = str_sub(cluster, start=-1))

    p_ramachandran <- ggplot(angles, aes(x=phi, y=psi, colour=cluster_num)) +
      geom_point() +
      facet_wrap(~wt)
  } else {
    angles <- NULL
    p_ramachandran <- NULL
  }
  
  # Cluster mean profiles
  mean_profiles <- group_by(tbl, cluster) %>%
    summarise_at(.vars = vars(A:Y), .funs = mean)

  mean_prof_long <- gather(mean_profiles, key='mut', value = 'norm_er', -cluster) %>%
    add_factor_order(cluster, mut, norm_er, sym=FALSE)
  
  # Cluster mean profile correlation
  cluster_cors <- transpose_tibble(mean_profiles, cluster, name_col = 'aa') %>%
    tibble_correlation(-aa) %>%
    rename(cluster1 = cat1, cluster2 = cat2) %>%
    mutate(wt1 = str_sub(cluster1, end = 1),
           wt2 = str_sub(cluster2, end = 1)) %>%
    left_join(as_tibble(BLOSUM62, rownames='wt1') %>%
                gather(key = 'wt2', value = 'BLOSUM62', -wt1) %>%
                filter(wt1 %in% Biostrings::AA_STANDARD, wt2 %in% Biostrings::AA_STANDARD),
              by=c('wt1', 'wt2')) %>%
    mutate(pair = mapply(function(x, y){str_c(str_sort(c(x, y)), collapse = '')}, wt1, wt2))
  
  cluster_mean_order <- levels(mean_prof_long$cluster)
  cluster_cor_order <- levels(cluster_cors$cluster1)
  
  mean_prof_long <- mutate(mean_prof_long, cluster = factor(cluster, levels = cluster_cor_order))
  
  p_mean_prof <- labeled_ggplot(
    p=ggplot(mean_prof_long, aes(x=mut, y=cluster, fill=norm_er)) +
    geom_tile() +
    scale_fill_gradient2() +
    coord_fixed() +
    ggtitle(str_c('Cluster centroid', er_str, 'for', cluster_str, 'clusters', sep=' ')) +
    guides(fill=guide_colourbar(title = er_str)) +
    theme(axis.ticks = element_blank(),
          panel.background = element_blank(),
          axis.title = element_blank(),
          axis.text.x = element_text(colour = AA_COLOURS[levels(mean_prof_long$mut)]),
          axis.text.y = element_text(colour = AA_COLOURS[str_sub(levels(mean_prof_long$cluster), end = 1)])),
    units = 'cm', width = 0.5*length(unique(mean_prof_long$mut)) + 4,
    height = 0.5*length(unique(mean_prof_long$cluster)) + 2, limitsize=FALSE)

  # Cluster Sizes
  cluster_sizes <- group_by(tbl, cluster) %>%
    summarise(n = n()) %>%
    mutate(aa = str_sub(cluster, end = 1),
           cluster = factor(cluster, levels = levels(mean_prof_long$cluster)))
  
  p_cluster_size <- labeled_ggplot(
    p = ggplot(cluster_sizes, aes(x=cluster, y=n, fill=aa)) +
    geom_col() +
    xlab('Cluster') +
    ylab('Size') +
    scale_fill_manual(values = AA_COLOURS) +
    scale_y_log10() +
    theme_pubclean() +
    guides(fill=FALSE) +
    theme(axis.text.x = element_text(colour = AA_COLOURS[str_sub(levels(cluster_sizes$cluster), end = 1)],
                                     angle = 90, hjust = 1, vjust = 0.5)),
    units = 'cm', height = 15, width = nrow(cluster_sizes) * 0.5 + 2)
  
  p_centre_cor <- labeled_ggplot(
    p = ggplot(cluster_cors, aes(x=cluster1, y=cluster2, fill=cor)) +
      geom_tile() +
      scale_fill_gradient2() +
      ggtitle(str_c('Correlation of', cluster_str, 'centroids for clusters based on', er_str, sep = ' ')) +
      coord_fixed() +
      theme(axis.ticks = element_blank(),
            panel.background = element_blank(),
            axis.title = element_blank(),
            axis.text.x = element_text(colour = AA_COLOURS[str_sub(levels(cluster_cors$cluster1), end = 1)], angle = 90, vjust = 0.5),
            axis.text.y = element_text(colour = AA_COLOURS[str_sub(levels(cluster_cors$cluster2), end = 1)])),
    units = 'cm', width = 0.5*length(levels(cluster_cors$cluster1)) + 2,
    height = 0.5*length(levels(cluster_cors$cluster2)) + 2, limitsize=FALSE)

  # Vs Blosum
  p_vs_blosum <- plot_cluster_profile_cor_blosum(cluster_cors, 'DE')

  # FoldX params
  if (!is.null(foldx)){
    tbl_fx <- group_by(foldx, !!id_col, !!pos_col, wt) %>%
      summarise_at(.vars = vars(-mut, -pdb_id, -sd), .funs = mean, na.rm=TRUE) %>%
      inner_join(tbl, ., by=c(id_col_str, pos_col_str, 'wt')) %>%
      select(cluster, !!id_col, !!pos_col, wt, total_energy:entropy_complex, everything())
    
    p_foldx_boxes <- labeled_ggplot(
    p=ggplot(gather(tbl_fx, key = 'term', value = 'ddG', total_energy:entropy_complex),
             aes(x=cluster, y=ddG, colour=wt)) +
      scale_colour_manual(values = AA_COLOURS) +
      geom_boxplot() +
      facet_wrap(~term, scales = 'free', ncol = 2) +
      guides(colour=FALSE) +
      ggtitle(str_c('FoldX energy term distribution for ', cluster_str, 'clusters (', er_str, ')')) +
      theme(panel.background = element_blank(),
            axis.title = element_blank(),
            axis.text.x = element_text(colour = AA_COLOURS[str_sub(unique(tbl_fx$cluster), end = 1)],
                                       angle = 90, vjust = 0.5)),
    units = 'cm', width = length(unique(tbl_fx$cluster)) + 5, height = 80)
    
    foldx_cluster_mean_energy <- gather(tbl_fx, key = 'foldx_term', value = 'ddG', total_energy:entropy_complex) %>%
      select(cluster, !!id_col, !!pos_col, wt, foldx_term, ddG, everything()) %>%
      group_by(cluster, foldx_term) %>%
      summarise(ddG = mean(ddG)) %>%
      group_by(foldx_term) %>%
      mutate(max_ddG = max(abs(ddG))) %>%
      filter(max_ddG != 0) %>% # Filter any terms that are all 0
      ungroup() %>%
      mutate(rel_ddG = ddG/max_ddG) %>%
      add_factor_order(cluster, foldx_term, rel_ddG, sym = FALSE)

  p_cluster_avg_foldx_profile <- labeled_ggplot(
    p=ggplot(foldx_cluster_mean_energy,
             aes(x=foldx_term, y=cluster, fill=rel_ddG)) +
      geom_tile() +
      scale_fill_gradient2() +
      ggtitle(str_c('Mean FoldX energy terms for each ', cluster_str, ' cluster (', er_str, ')')) +
      coord_fixed() +
      theme(plot.title = element_text(hjust = 0.5, size=8),
            axis.ticks = element_blank(),
            panel.background = element_blank(),
            axis.title = element_blank(),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            axis.text.y = element_text(colour = AA_COLOURS[str_sub(levels(foldx_cluster_mean_energy$cluster), end = 1)])),
    units='cm', height=0.4 * length(unique(tbl_fx$cluster)) + 5, width=13, limitsize=FALSE)
  } else {
    tbl_fx <- NULL
    foldx_cluster_mean_energy <- NULL
    p_foldx_boxes <- NULL
    p_cluster_avg_foldx_profile <- NULL
  }
  
  return(list(angles=angles,
              cluster_sizes=cluster_sizes,
              mean_profiles=mean_profiles,
              cluster_cor_order=cluster_cor_order,
              cluster_mean_order=cluster_mean_order, 
              foldx=tbl_fx,
              foldx_cluster_mean_energy=foldx_cluster_mean_energy,
              plots=list(cluster_sizes=p_cluster_size,
                         ramachandran=p_ramachandran,
                         mean_profiles=p_mean_prof,
                         mean_profile_vs_blosum=p_vs_blosum,
                         mean_profile_cor=p_centre_cor,
                         foldx_term_distribution=p_foldx_boxes,
                         foldx_cluster_mean_profile=p_cluster_avg_foldx_profile)))


  

}

plot_cluster_profile_cor_blosum <- function(cluster_cors, aa_pair='DE'){
  return(
    ggplot(filter(cluster_cors, cor < 1), aes(x=cor, y=BLOSUM62)) +
      geom_point(aes(colour='All')) +
      geom_point(aes(colour=aa_pair), filter(cluster_cors, cor < 1, pair == aa_pair)) +
      geom_smooth(method = 'lm') +
      scale_colour_manual(values = structure(c('red', 'black'), names=c(aa_pair, 'All'))) +
      guides(colour = guide_legend(title = 'AA Pair'))
  )
}



########