#!/usr/bin/env Rscript 
# Functions to perform clustering analysis on per position mutational profiles from deep mutagenesis studies


# Generate PCA of mutational profiles
positional_profile_PCA <- function(variant_matrix){
  pca <- prcomp(as.matrix(select(variant_matrix, A:Y)), center = TRUE, scale. = TRUE)
  pca_variants <- bind_cols(select(variant_matrix, -(A:Y)), as_tibble(pca$x))
  
  return(list(variants=pca_variants, pca=pca))
}

basic_pca_plots <- function(pca){
  plots <- list()
  plots$all_pcs <- plot_all_pcs(pca$variants, colour_var = 'wt')
  
  
  plots$by_authour <- ggplot(pca$variants, aes(x=PC1, y=PC2, colour=gene_name)) + 
    facet_wrap(~study) +
    geom_point()
  
  plots$secondary_structure <- plot_all_pcs(pca$variants, colour_var = 'ss')
  plots$secondary_structure_reduced <- plot_all_pcs(pca$variants, colour_var = 'ss_reduced')
  
  plots$fields_group_studies <- ggplot(filter(pca$variants,
                                              authour %in% c('Araya et al.', 'Melamed et al.', 'Starita et al.',
                                                             'Kitzman et al.', 'Weile et al.')),
                                       aes(x=PC1, y=PC2, colour=study)) +
    geom_point()
  
  plots$position_sig <- ggplot(pca$variants, aes(x=PC1, y=PC2, colour=sig_count)) +
    geom_point()
  
  plots$by_aa <- ggplot(pca$variants, aes(x=PC1, y=PC2, colour=gene_name)) + 
    facet_wrap(~wt) +
    geom_point()
  
  plots$surface_accesibility <- ggplot(pca$variants, aes(x=PC1, y=PC2, colour=all_atom_rel)) +
    geom_point() +
    scale_colour_gradientn(colours = c('blue', 'green', 'yellow', 'orange', 'red'))
  
  return(plots)
}

get_avg_aa_pca_profile <- function(pca){
  avg_profile <- pca$variants %>%
    group_by(wt) %>%
    summarise_at(.vars = vars(starts_with('PC')), .funs = list(~ mean(.)))
  
  cor_mat <- select(avg_profile, -wt) %>% 
    t() %>% 
    set_colnames(avg_profile$wt) %>%
    cor()
  
  aa_order <- rownames(cor_mat)[hclust(dist(cor_mat))$order]
  
  cor_tbl <- cor_mat %>%
    as_tibble(rownames = 'AA1') %>%
    gather(key = 'AA2', value = 'cor', -AA1) %>%
    mutate(AA1 = factor(AA1, levels = aa_order),
           AA2 = factor(AA2, levels = aa_order))
  
  return(list(avg_profile=avg_profile, cor_mat=cor_mat, cor_tbl=cor_tbl, aa_order=aa_order))
}

pca_surf_acc_cor <- function(pca){
  variants <- drop_na(pca$variants, all_atom_abs:polar_rel)
  
  variant_mat <- select(variants, starts_with('PC')) %>% 
    as.matrix()
  
  surface_acc_mat <- select(variants, all_atom_abs:polar_rel) %>% 
    as.matrix()
  
  cor_mat <- cor(variant_mat, surface_acc_mat)
  
  cor_tbl <- cor_mat %>%
    as_tibble(rownames = 'PC') %>%
    gather(key = 'Surf', value = 'cor', -PC) %>%
    mutate(PC = factor(PC, levels = str_c('PC', 1:dim(cor_mat)[1])))
  
  return(list(tbl=cor_tbl, matrix=cor_mat))
}

aa_avg_profile_plot <- function(x){list(avg_aa_profile=ggplot(x$avg_profile, aes(x=PC1, y=PC2, label=wt)) + geom_text())}

aa_profile_heatmap <- function(pca){list(
  aa_profile_heatmap=ggplot(pca$cor_tbl, aes(x=AA1, y=AA2, fill=cor)) + 
    geom_tile(colour='white') + 
    scale_fill_gradient2() +
    theme(axis.ticks = element_blank(), panel.background = element_blank())
)}

pca_surface_heatmap <- function(pca){list(
  pc_surf_acc_heatmap=ggplot(pca$tbl, aes(x=PC, y=Surf, fill=cor)) + 
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