#!/usr/bin/env Rscript 
# Functions to analyse chemical environments with respect to deep mutational scanning data

#### Utility ####
# generate a reduced AA profile (with similar AAs grouped) from a 20 long AA count vector
sorted_aa_1_code <- sort(Biostrings::AA_STANDARD)
reduce_aa_profile <- function(x){
  names(x) <- sorted_aa_1_code
  return(sapply(AA_REDUCED_CLASSES, function(aa_group){sum(x[aa_group])}))
}

# Retrieve column names for a chem env profile column
get_profile_names <- function(tbl, col){
  col <- enquo(col)
  
  first_profile <- pull(tbl, !!col)[[1]]
  names <- names(first_profile)
  if (is.null(names)){
    names <- str_c(rlang::quo_text(col), '_', 1:length(first_profile))
  }
  return(names)
}

# Expand a profile column into columns with given names
expand_profile_column <- function(tbl, col, names=NULL){
  col <- enquo(col)
  
  names <- if(is.null(names)) get_profile_names(tbl, !!col) else names
  
  prof <- pull(tbl, !!col) %>%
    do.call(rbind, .) %>%
    set_colnames(names) %>%
    as_tibble()
  
  return(bind_cols(select(tbl, -!!col), prof))
}
########

#### Direct Profile analysis ####
plot_basic_profile_analysis <- function(tbl, prof_col, names=NULL){
  prof_col <- enquo(prof_col)
    
  names <- if(is.null(names)) get_profile_names(tbl, !!prof_col) else names
  
  tbl <- expand_profile_column(tbl, !!prof_col, names=names)
  num_cats <- length(names)
  p_paired = labeled_ggplot(
    ggpairs(tbl, columns = names,
            lower = list(continuous=function(d, m, ...){ggplot(d, m, ...) + geom_count()})),
    height = num_cats * 2, width = num_cats * 2, limitsize=FALSE)
  
  cor_mat <- tibble_to_matrix(tbl, columns = !!names) %>%
    cor()
  
  group_order <- rownames(cor_mat)[hclust(dist(cor_mat))$order]
  
  cors <- as_tibble(cor_mat, rownames = 'cat1') %>%
    gather(key = 'cat2', value = 'cor', -cat1) %>%
    mutate(cat1 = factor(cat1, levels = group_order),
           cat2 = factor(cat2, levels = group_order))
  
  p_cor <- ggplot(cors, aes(x=cat1, y=cat2, fill=cor)) +
      geom_tile() +
      xlab('') +
      ylab('') +
      scale_fill_gradient2() +
      theme(axis.ticks = element_blank(), panel.background = element_blank())
  
  return(list(scatter_pairs = p_paired,
              cor_heatmap = p_cor))
}

# Linear model(s) of ER ~ profile
calc_profile_lm <- function(tbl, target_col, ..., include_sig_count=FALSE){
  target_col <- enquo(target_col)
  prof_cols <- enquos(...)
  
  if (include_sig_count){
    return(
      select(tbl, !!target_col, sig_count, !!! prof_cols) %>%
        drop_na(!!target_col) %>%
        lm(as.formula(str_c(rlang::quo_text(target_col), '~ . + 0')), data = .)
    )
  } else {
    return(
      select(tbl, !!target_col, !!! prof_cols) %>%
        drop_na(!!target_col) %>%
        lm(as.formula(str_c(rlang::quo_text(target_col), '~ .')), data = .)
    )
  }
}

calc_all_profile_lms <- function(tbl, prof_col, target_cols, prof_col_names=NULL, include_sig_count=FALSE){
  prof_col <- enquo(prof_col)
  target_cols <- enquo(target_cols)
  
  prof_col_names <- if(is.null(prof_col_names)) get_profile_names(tbl, !!prof_col) else prof_col_names
  tbl <- expand_profile_column(tbl, !!prof_col, names=prof_col_names)
  all_study_lms <- gather(tbl, key = 'mut_aa', value = 'er', !!target_cols) %>%
    group_by(mut_aa) %>%
    do(model=calc_profile_lm(., er, !!!syms(prof_col_names), include_sig_count = include_sig_count),
       n=sum(!is.na(.$er))) %>%
    mutate(study = 'ALL')
  
  per_study_lms <- gather(tbl, key = 'mut_aa', value = 'er', !!target_cols) %>%
    group_by(mut_aa, study) %>%
    do(model=calc_profile_lm(., er, !!!syms(prof_col_names), include_sig_count = include_sig_count),
       n=sum(!is.na(.$er)))
  
  return(
    bind_rows(all_study_lms, per_study_lms) %>%
      unnest(n) %>%
      mutate(summary = lapply(model, glance)) %>%
      unnest(summary)
  )
}
########

#### PCA of profiles ####
# Plot basic analyses of PCA components
plot_chem_env_basic_pca_plots <- function(pca, discrete_factors=c('aa', 'ss'),
                                          cont_factors=c('all_atom_rel', 'relative_position'),
                                          ...){
  scatter_plots <- sapply(c(discrete_factors, cont_factors),
                          function(x){plot_all_pcs(pca$profiles, colour_var = x, ...)},
                          simplify = FALSE) %>%
    set_names(str_c(names(.), '_pca'))
  
  avg_profile_heatmaps <- sapply(discrete_factors,
                                 function(x){plot_avg_factor_pca_profile(pca, variable = x)},
                                 simplify = FALSE) %>%
    set_names(str_c(names(.), '_pca_avg_heatmap'))
  
  return(c(scatter_plots, avg_profile_heatmaps))
}

# Calc PCA
chem_env_pca <- function(tbl, var='nearest_10', names=NULL){
  if (is.null(names)){
    names <- names(pull(tbl, !!var)[[1]])
    if (is.null(names)){
      stop('Profile vectors have no names and none are supplied')
    }
  }
  
  pca <- pull(tbl, !!var) %>%
    do.call(rbind, .) %>%
    set_colnames(names) %>%
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
  
  clust <- hclust(dist(tibble_to_matrix(avg_prof, starts_with('PC'), row_names = avg_prof[[variable]])))
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
chem_env_tsne <- function(tbl, var='nearest_10', names=NULL, ...){
  if (is.null(names)){
    names <- names(pull(tbl, !!var)[[1]])
    if (is.null(names)){
      stop('Profile vectors have no names and none are supplied')
    }
  }
  
  mat <- pull(tbl, !!var) %>%
    do.call(rbind, .) %>%
    set_colnames(names)
  
  mat_dupe_rows <- enumerate_unique_rows(mat)
  mat_deduped <- mat[!mat_dupe_rows$duplicate,]
  
  tsne <- Rtsne(mat_deduped, ...)
  
  out_tbl <- bind_cols(tbl, as_tibble(set_colnames(tsne$Y[mat_dupe_rows$indeces,], c('tSNE1', 'tSNE2'))))
  return(list(profiles=out_tbl, tsne=tsne, dupe_rows=mat_dupe_rows$duplicate, unique_row_indeces=mat_dupe_rows$indeces))
}
########