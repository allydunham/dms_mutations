#!/usr/bin/env Rscript 
# Functions to perform clustering analysis on per position mutational profiles from deep mutagenesis studies

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
  profiles <- drop_na(pca$profiles, !!!.vars)
  
  pcas_mat <- select(profiles, starts_with('PC')) %>% 
    as.matrix()
  
  factor_mat <- select(profiles, !!!.vars) %>% 
    as.matrix()
  
  cor_mat <- cor(pcas_mat, factor_mat)
  
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
make_hclust_clusters <- function(tbl, cols, dist_method = 'manhattan', h = NULL, k = NULL, ...){
  cols <- enquo(cols)
  
  mat <- tibble_to_matrix(tbl, !!cols)
  
  hc <- hclust(dist(mat, method = dist_method), ...)
  
  return(list(tbl = mutate(tbl, cluster = cutree(hc, k = k, h = h)),
              hclust = hc))
}
########