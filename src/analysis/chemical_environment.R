#!/usr/bin/env Rscript 
# Functions to analyse chemical environments with respect to deep mutational scanning data

SORTED_AA_1_CODE <- sort(Biostrings::AA_STANDARD)

#### Utility ####
# generate a reduced AA profile (with similar AAs grouped) from a 20 long AA count vector
reduce_aa_profile <- function(x){
  names(x) <- SORTED_AA_1_CODE
  return(sapply(AA_REDUCED_CLASSES, function(aa_group){sum(x[aa_group])}))
}

# Expand a profile column into columns with given names
expand_profile_column <- function(tbl, col, names=NULL, prof_name=FALSE){
  col <- enquo(col)
  
  names <- if(is.null(names)) get_profile_names(tbl, !!col, prof_name) else names
  
  prof <- pull(tbl, !!col) %>%
    do.call(rbind, .) %>%
    set_colnames(names) %>%
    as_tibble()
  
  return(bind_cols(select(tbl, -!!col), prof))
}

# Retrieve column names for a chem env profile column
get_profile_names <- function(tbl, col, prof_name=FALSE){
  col <- enquo(col)
  
  first_profile <- pull(tbl, !!col)[[1]]
  names <- names(first_profile)
  if (is.null(names)){
    names <- str_c(if(prof_name) rlang::quo_text(col) else 'profile', '_', 1:length(first_profile))
  }
  return(names)
}


########

#### Direct Profile analysis ####
########

#### Linear Model ####
# Linear model(s) of ER ~ profile
calc_profile_lm <- function(tbl, target_col, ..., include_intercept=FALSE){
  target_col <- enquo(target_col)
  prof_cols <- enquos(...)
  
  if (!include_intercept){
    return(
      select(tbl, !!target_col, !!! prof_cols) %>%
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

# Currently must use vars to define multiple targets/profile columns
calc_all_profile_lms <- function(tbl, prof_vars, target_vars,
                                 include_intercept=FALSE, per_study=FALSE){
  
  lms <- gather(tbl, key = 'target', value = 'er', !!! target_vars) %>%
    group_by(target) %>%
    do(model=calc_profile_lm(., er, !!! prof_vars, include_intercept = include_intercept),
       n=sum(!is.na(.$er)))
  
  if (per_study){
    lms <- mutate(lms, study == 'ALL') %>%
      bind_rows(.,
                gather(tbl, key = 'target', value = 'er', !!! target_vars) %>%
                  group_by(target, study) %>%
                  do(model=calc_profile_lm(., er, !!! prof_vars, include_intercept = include_intercept),
                     n=sum(!is.na(.$er)))
    )
  }
  
  return(
      unnest(lms, n) %>%
      mutate(summary = lapply(model, glance)) %>%
      unnest(summary)
  )
}

# Plot basic LM analysis
basic_prof_lm_plots <- function(lm_tbl){
  p_summary <- labeled_ggplot(
    ggplot(lm_tbl, aes(x=target, y=r.squared, fill=-log10(p.value), label=str_c('n = ', n))) + 
      geom_col() + 
      geom_text(colour='white', y = 0 + max(lm_tbl$r.squared)/100, check_overlap = FALSE, angle=90, hjust=0, vjust=0.5),
    width=15, height=10)
  
  preds <- pull(lm_tbl, model) %>%
    lapply(augment) %>%
    set_names(lm_tbl$target) %>%
    bind_rows(.id = 'aa')
  
  p_predictions <- ggplot(preds, aes(x=er, y=.fitted, colour=.resid)) + 
    geom_point() + 
    facet_wrap(~aa) +
    guides(colour=guide_colourbar(title='Residual')) +
    ylab('Predicted ER') +
    xlab('ER')
  
  loadings <- select(lm_tbl, target, model, n, r.squared) %>%
    mutate(coef_df = lapply(model, tidy)) %>%
    unnest(coef_df)
  
  p_loadings <- ggplot(loadings, aes(x=target, y=term, fill=estimate)) + 
    geom_tile() +
    scale_fill_gradient2() +
    theme(axis.ticks = element_blank(), panel.background = element_blank()) +
    xlab('Substituted AA') +
    ylab('LM Term') +
    geom_point(data = filter(loadings, p.value < 0.0001), aes(shape='p < 0.0001')) + 
    geom_point(data = filter(loadings, p.value < 0.001, p.value > 0.0001), aes(shape='p < 0.001')) + 
    geom_point(data = filter(loadings, p.value < 0.01, p.value > 0.001), aes(shape='p < 0.01')) +
    scale_shape_manual(values = c('p < 0.0001'=8, 'p < 0.001'=3, 'p < 0.01'=20))

  return(list(fit_summary=p_summary,
              predictions=p_predictions,
              loadings=p_loadings))
}
########

#### PCA of profiles ####
# Plot basic analyses of PCA components
plot_chem_env_basic_pca_plots <- function(tbl, pca,
                                          discrete_factors=c('aa', 'ss'),
                                          cont_factors=c('all_atom_rel', 'relative_position'),
                                          ...){
  scatter_plots <- sapply(c(discrete_factors, cont_factors),
                          function(x){plot_all_pcs(tbl, colour_var = x, ...)},
                          simplify = FALSE)
  
  avg_profile_heatmaps <- sapply(discrete_factors,
                                 function(x){plot_avg_factor_pca_profile(tbl, variable = x)},
                                 simplify = FALSE) %>%
    set_names(str_c(names(.), '_avg_heatmap'))
  
  return(c(scatter_plots, avg_profile_heatmaps))
}

# Avg PCA profile against a factor
plot_avg_factor_pca_profile <- function(tbl, variable='aa'){
  variable_sym <- sym(variable)
  avg_prof <- tbl %>%
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
chem_env_tsne <- function(tbl, ..., tsne_kwargs=list()){
  prof_cols <- enquos(...)
  
  mat <- tibble_to_matrix(tbl, !!! prof_cols)
  
  mat_dupe_rows <- enumerate_unique_rows(mat)
  mat_deduped <- mat[!mat_dupe_rows$duplicate,]
  
  tsne <- do.call(Rtsne, c(list(X=mat_deduped), tsne_kwargs))

  return(list(tsne=tsne, dupe_rows=mat_dupe_rows$duplicate, unique_row_indeces=mat_dupe_rows$indeces))
}
########
