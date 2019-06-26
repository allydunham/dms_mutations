#!/usr/bin/env Rscript 
# Functions to analyse chemical environments with respect to deep mutational scanning data

#### PCA of profiles ####
# Plot basic analyses of PCA components
plot_chem_env_basic_pca_plots <- function(pca, discrete_factors=c('aa', 'ss'),
                                          cont_factors=c('all_atom_rel', 'relative_position')){
  scatter_plots <- sapply(c(discrete_factors, cont_factors), function(x){plot_all_pcs(pca$profiles, colour_var = x)}, simplify = FALSE) %>%
    set_names(str_c(names(.), '_pca'))
  
  avg_profile_heatmaps <- sapply(discrete_factors, function(x){plot_avg_factor_pca_profile(pca, variable = x)}, simplify = FALSE) %>%
    set_names(str_c(names(.), '_pca_avg_heatmap'))
  
  return(c(scatter_plots, avg_profile_heatmaps))
}

# Calc PCA
chem_env_pca <- function(tbl, var='nearest_10'){
  pca <- pull(tbl, !!var) %>%
    do.call(rbind, .) %>%
    set_colnames(sort(Biostrings::AA_STANDARD)) %>%
    prcomp()
  
  out_tbl <- bind_cols(tbl, as_tibble(pca$x))
  return(list(profiles=out_tbl, pca=pca))
}

# Avg PCA profile against a factor
plot_avg_factor_pca_profile <- function(pca, variable='aa'){
  variable_sym <- sym(variable)
  avg_prof <- pca$profiles %>%
    group_by(!!variable_sym) %>%
    summarise_at(.vars = vars(starts_with('PC')), .funs = list(~ mean(.)))
  
  clust <- hclust(dist(tibble_to_matrix(avg_prof, PC1:PC20, row_names = avg_prof[[variable]])))
  variable_order <- clust$labels[clust$order]
  pc_order <- str_c('PC', 1:(ncol(avg_prof)-1))
    
  avg_prof_long <- gather(avg_prof, key='PC', value='value', -!!variable_sym) %>%
    mutate(!!variable_sym := factor(!!variable_sym, levels=variable_order),
           PC = factor(PC, levels=pc_order))
    
  return(ggplot(avg_prof_long, aes_string(x=variable, y='PC', fill='value')) +
           geom_tile() +
           scale_fill_gradient2() +
           theme(axis.ticks = element_blank(), panel.background = element_blank()))
}

########

#### tSNE ####
chem_env_tsne <- function(tbl, var='nearest_10', ...){
  mat <- pull(tbl, !!var) %>%
    do.call(rbind, .) %>%
    set_colnames(sort(Biostrings::AA_STANDARD))
  
  mat_dupe_rows <- enumerate_unique_rows(mat)
  mat_deduped <- mat[!mat_dupe_rows$duplicate,]
  
  tsne <- Rtsne(mat_deduped, ...)
  
  out_tbl <- bind_cols(tbl, as_tibble(set_colnames(tsne$Y[mat_dupe_rows$indeces,], c('tSNE1', 'tSNE2'))))
  return(list(profiles=out_tbl, tsne=tsne, dupe_rows=mat_dupe_rows$duplicate, unique_row_indeces=mat_dupe_rows$indeces))
}
########