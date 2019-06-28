#!/usr/bin/env Rscript 
# Functions to analyse chemical environments with respect to deep mutational scanning data

SORTED_AA_1_CODE <- sort(Biostrings::AA_STANDARD)

## Overall analysis of a chemical environment profile
analyse_chem_env_profile <- function(chem_env, prof_col, prof_col_names=NULL){
  prof_col <- enquo(prof_col)

  # Fetch profile column names if none given
  prof_col_names <- if(is.null(prof_col_names)) get_profile_names(tbl, !!prof_col) else prof_col_names
  prof_col_syms <- syms(prof_col_names)
  
  chem_env <- expand_profile_column(chem_env, !!prof_col, names=prof_col_names) %>%
    mutate(duplicate_position = duplicated(select(chem_env, pdb_id, chain, position, aa)))
  
  ## Basic analysis plots
  basic_plots <- plot_basic_profile_analysis(filter(chem_env, !duplicate_position), !!! prof_col_syms) %>%
      set_names(str_c('profile_', names(.)))
  
  ## PCA analysis of profile
  pca <- tibble_pca(filter(chem_env, !duplicate_position), !!! prof_col_syms)
  chem_env <- bind_cols(chem_env,
                        as_tibble(scale(select(chem_env, !!! prof_col_syms), pca$center, pca$scale) %*% pca$rotation))
  
  num_pcs <- length(grep('^PC[0-9]*$', names(chem_env)))
  max_plot_pc <- if(num_pcs %% 2 == 0) num_pcs else num_pcs - 1 # Work around to current rigid PC plotting
  plot_row_cols <- get_good_rows_cols(max_plot_pc/2)
  pca_plots <- plot_chem_env_basic_pca_plots(chem_env,
                                             max_pc = max_plot_pc,
                                             nrow = plot_row_cols[1], ncol = plot_row_cols[2],
                                             cont_factors=c('all_atom_rel', 'relative_position', 'sig_count'),
                                             discrete_factors=c('ss', 'ss_reduced', 'aa', 'aa_reduced',
                                                                'pdb_id', 'gene_name', 'species'))
  
  pca_factor_cors <- pca_factor_cor(list(pca=pca, profiles=chem_env),
                                    .vars = vars(all_atom_abs:polar_rel, relative_position, sig_count))
                                                                        
  pca_plots <- c(pca_plots,
                 list(factor_heatmap=pca_factor_heatmap(pca_factor_cors)))
  
  
  ## tSNE analysis of profile
  tsne <- chem_env_tsne(chem_env, !!! prof_col_syms)
  
  chem_env <- bind_cols(chem_env,
                        as_tibble(set_colnames(tsne$tsne$Y[tsne$unique_row_indeces,], c('tSNE1', 'tSNE2'))))      
                  
  tsne_plots <- plot_factors(filter(chem_env, !duplicate_position),
                             tSNE1, tSNE2, quos(ss_reduced=ss_reduced, aa_reduced=aa_reduced, gene_name=gene_name,
                                                sig_count=sig_count, sqrt_suf_acc=sqrt(all_atom_rel)))
  
  ## LM analysis
  # Just profile
  prof_lm <- calc_all_profile_lms(chem_env, prof_vars = vars(!!! prof_col_syms), target_vars = vars(A:Y),
                                  include_intercept = FALSE)
  
  lm_plots <- list()
  lm_plots$lm_summary <- labeled_ggplot(
    ggplot(prof_lm, aes(x=target, y=r.squared, fill=-log10(p.value), label=n)) + 
      facet_wrap(~study) + 
      geom_col() + 
      geom_text(colour='red', check_overlap = FALSE, angle=90, hjust=0, vjust=0.5, nudge_y = 0.05),
    width=15, height=10)
  
  prof_lm_preds <- filter(prof_lm, study=='ALL') %>%
    pull(model) %>%
    lapply(augment) %>%
    set_names(filter(prof_lm, study=='ALL') %>% pull(target)) %>%
    bind_rows(.id = 'aa')
    
  lm_plots$lm_predictions <- ggplot(prof_lm_preds, aes(x=er, y=.fitted, colour=.se.fit)) + 
    geom_point() + 
    facet_wrap(~aa) +
    ylab('Predicted ER') +
    xlab('ER')
  
  # significance of position added
  prof_lm_sig_count <- calc_all_profile_lms(chem_env, prof_vars = vars(!!! prof_col_syms, sig_count), target_vars = vars(A:Y),
                                            include_intercept = TRUE)
  
  lm_plots$lm_sig_count_summary <- labeled_ggplot(
    ggplot(prof_lm_sig_count, aes(x=target, y=r.squared, fill=-log10(p.value), label=n)) + 
      facet_wrap(~study) + 
      geom_col() + 
      geom_text(colour='red', check_overlap = FALSE, angle=90, hjust=0, vjust=0.5, nudge_y = 0.05),
    width=15, height=10)
  
  prof_lm_sig_count_preds <- filter(prof_lm_sig_count, study=='ALL') %>%
    pull(model) %>%
    lapply(augment) %>%
    set_names(filter(prof_lm, study=='ALL') %>% pull(target)) %>%
    bind_rows(.id = 'aa')
  
  lm_plots$lm_sig_count_predictions <- ggplot(prof_lm_sig_count_preds, aes(x=er, y=.fitted, colour=.se.fit)) + 
    geom_point() + 
    facet_wrap(~aa) +
    ylab('Predicted ER') +
    xlab('ER')

  prof_lm_sig_count_loadings <- filter(prof_lm_sig_count, study == 'ALL') %>%
    select(target, model, n, r.squared) %>%
    mutate(coef_df = lapply(model, tidy)) %>%
    unnest(coef_df)
  
  lm_plots$lm_sig_count_loadings <- ggplot(prof_lm_sig_count_loadings, aes(x=target, y=term, fill=estimate)) + 
    geom_tile() +
    scale_fill_gradient2() +
    theme(axis.ticks = element_blank(), panel.background = element_blank()) +
    xlab('Substituted AA') +
    ylab('LM Term') +
    geom_point(data = filter(prof_lm_sig_count_loadings, p.value < 0.0001), aes(shape='p < 0.0001')) + 
    geom_point(data = filter(prof_lm_sig_count_loadings, p.value < 0.001, p.value > 0.0001), aes(shape='p < 0.001')) + 
    geom_point(data = filter(prof_lm_sig_count_loadings, p.value < 0.01, p.value > 0.001), aes(shape='p < 0.01')) +
    scale_shape_manual(values = c('p < 0.0001'=8, 'p < 0.001'=3, 'p < 0.01'=20))
  
  return(
    list(tbl=chem_env,
         pca=pca,
         pca_factor_corelation=pca_factor_cors,
         tsne=tsne,
         lm=prof_lm,
         lm_sig_count=prof_lm_sig_count,
         plots=c(basic_plots,
                 list(pca=pca_plots,
                      tSNE=tsne_plots,
                      lm=lm_plots)
         )
    )
  )
}

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
# ... comprises the tibble columns to make up the profile
plot_basic_profile_analysis <- function(tbl, ...){
  prof_cols <- enquos(...)
    
  num_cats <- length(prof_cols)
  p_paired = labeled_ggplot(
    ggpairs(tbl, columns = sapply(prof_cols, rlang::as_name),
            lower = list(continuous=function(d, m, ...){ggplot(d, m, ...) + geom_count()})),
    height = num_cats * 2, width = num_cats * 2, limitsize=FALSE)
  
  cor_mat <- tibble_to_matrix(tbl, !!! prof_cols) %>%
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
########

#### Linear Model ####
# Linear model(s) of ER ~ profile
calc_profile_lm <- function(tbl, target_col, ..., include_intercept=FALSE){
  target_col <- enquo(target_col)
  prof_cols <- enquos(...)
  
  if (include_intercept){
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
                                 include_intercept=FALSE){
  
  all_study_lms <- gather(tbl, key = 'target', value = 'er', !!! target_vars) %>%
    group_by(target) %>%
    do(model=calc_profile_lm(., er, !!! prof_vars, include_intercept = include_intercept),
       n=sum(!is.na(.$er))) %>%
    mutate(study = 'ALL')
  
  per_study_lms <- gather(tbl, key = 'target', value = 'er', !!! target_vars) %>%
    group_by(target, study) %>%
    do(model=calc_profile_lm(., er, !!! prof_vars, include_intercept = include_intercept),
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
