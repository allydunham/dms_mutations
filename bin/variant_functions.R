#!/usr/bin/env Rscript 
# Script containing functions for loading and analysis of data from Envision, SIFT, FoldX, Polyphen2 & EVCouplings

#### Plotting function ####
## Plot deep mutagenesis data and variant effect predictions
# Generic
plot_predictions <- function(x, ...){
  UseMethod('plot_predictions', x)
}

# Single Variants
plot_predictions.single_variant <- function(x){
  study <- str_c(x$dm$authour, ' ', x$dm$year, ': ', x$dm$gene_name)
  print(str_c('Generating plots for ', study))
  
  plots <- sapply(names(x$foldx), function(id){plot_foldx(x$single_variants, id, study)}, simplify = FALSE) %>%
               unname(.) %>%
               unlist(., recursive = FALSE)
  if ('sift_score' %in% names(x$single_variants)){
    plots <- c(plots, plot_sift(x$single_variants, study))
  }
  if ('pph2_prob' %in% names(x$single_variants)){
    plots <- c(plots, plot_pph(x$single_variants, study))
  }
  if ('envision_prediction' %in% names(x$single_variants)){
    plots <- c(plots, plot_envision(x$single_variants, study))
  }
  if ('evcoup_epistatic' %in% names(x$single_variants)){
    plots <- c(plots, plot_evcoup(x$single_variants, study))
  }
  
  plots <- c(plots, plot_misc(x$single_variants, study))
  
  return(plots)
}

# Multi Variants
plot_predictions.multi_variant <- function(x){
  study <- str_c(x$dm$authour, ' ', x$dm$year, ': ', x$dm$gene_name)
  print(str_c('Generating plots for ', study))

  plots <- sapply(names(x$foldx), function(id){plot_foldx_multi(x$multi_variants, id, study)}, simplify = FALSE) %>%
               unname(.) %>%
               unlist(., recursive = FALSE)
  
  if ('sift_score' %in% names(x$single_variants)){
    plots <- c(plots, plot_sift(x$single_variants, study))
  }
  if ('pph2_prob' %in% names(x$single_variants)){
    plots <- c(plots, plot_pph(x$single_variants, study))
  }
  if ('envision_prediction' %in% names(x$single_variants)){
    plots <- c(plots, plot_envision(x$single_variants, study))
  }
  if ('evcoup_epistatic' %in% names(x$multi_variants)){
    plots <- c(plots, plot_evcoup_multi(x$multi_variants, study))
  }
  
  plots <- c(plots, plot_misc(x$single_variants, study))
  
  return(plots)
}

# Plot SIFT Scores
# Expects a df with columns: score, raw_score, sift_score, sift_prediction
plot_sift <- function(tbl, study=''){
  # Score vs Sift score
  p_score <- ggplot(tbl, aes(x=score, y=sift_score)) +
    xlab(MUT_SCORE_NAME) + 
    ylab('SIFT Score') +
    ggtitle(study)
  
  if (dim(tbl)[1] < 1000){
    p_score <- p_score + geom_point() 
  } else {
    p_score <- p_score + geom_bin2d()
  }
  
  p_raw_score <- p_score + aes(x=raw_score) + xlab(RAW_MUT_SCORE_NAME)
  
  p_log_score <- p_score + scale_y_log10()
  
  # Sift predictions vs score
  p_prediction <- ggplot(drop_na(tbl, sift_prediction), aes(x=sift_prediction, y=score)) +
    geom_boxplot() +
    stat_summary(geom = 'text', fun.data = function(x){return(c(y = -2.5, label = length(x)))}) +
    stat_compare_means(comparisons = list(c('DELETERIOUS','TOLERATED'))) +
    #stat_summary(geom = 'text', colour='red', fun.data = function(x){return(c(y = mean(x), label = signif(mean(x), 3)))}) +
    ggtitle(study) +
    xlab('SIFT Prediction') +
    ylab(MUT_SCORE_NAME)
 
  return(list(score_vs_sift=p_score,
              raw_score_vs_sift=p_raw_score,
              sift_prediction_vs_score=p_prediction,
              log2_sift_vs_score=p_log_score))
}

# Plot Envision Scores
# Expects df with columns: score, raw_score, envision_prediction, log2_envision_prediction 
plot_envision <- function(tbl, study=''){
  p_score <- ggplot(tbl, aes(x=score, y=envision_prediction)) +
    xlab(MUT_SCORE_NAME) + 
    ylab('Envision Prediction') +
    ggtitle(study)
  
  if (dim(tbl)[1] < 1000){
    p_score <- p_score + geom_point() 
  } else {
    p_score <- p_score + geom_bin2d()
  }
  
  p_raw_score <- p_score + aes(x=raw_score) + xlab(RAW_MUT_SCORE_NAME)
  
  p_log2 <- p_score + aes(y=log2_envision_prediction) + ylab('Log2 Envision Prediction')
  
  return(list(score_vs_envision=p_score,
              raw_score_vs_envision=p_raw_score,
              score_vs_log2_envision=p_log2))
}

# Plot Polyphen2 Scores
# Expects df with columns score, raw_score, pph2_prediction, pph2_class, pph2_prob, pph2_FPR, pph2_TPR, pph2_FDR
plot_pph <- function(tbl, study=''){
  p_score <- ggplot(tbl, aes(x=score, y=pph2_prob)) +
    xlab(MUT_SCORE_NAME) + 
    ylab('PolyPhen2 Probability') +
    ggtitle(study)
  
  if (dim(tbl)[1] < 1000){
    p_score <- p_score + geom_point() 
  } else {
    p_score <- p_score + geom_bin2d()
  }
  
  p_raw_score <- p_score + aes(x=raw_score) + xlab(RAW_MUT_SCORE_NAME)
  
  p_class <- ggplot(drop_na(tbl, pph2_class), aes(x=pph2_class, y=score)) +
    geom_boxplot() +
    stat_summary(geom = 'text', fun.data = function(x){return(c(y = -2.5, label = length(x)))}) +
    stat_compare_means(comparisons = list(c('deleterious','neutral'))) +
    #stat_summary(geom = 'text', colour='red', fun.data = function(x){return(c(y = mean(x), label = signif(mean(x), 3)))}) +
    ggtitle(study) +
    xlab('PolyPhen2 Class') +
    ylab(MUT_SCORE_NAME)
  
  p_prediction <- p_class +
    aes(x=pph2_prediction) + 
    stat_compare_means(comparisons = list(c("probably damaging", "possibly damaging", "benign"))) +
    xlab('PolyPhen2 Prediction')
  
    
  p_pph_vs_sift <- ggplot(drop_na(tbl, pph2_prediction, sift_prediction),
                          aes(x=pph2_prediction, y=..count.., fill=sift_prediction)) +
    geom_bar(position = 'dodge') +
    xlab('Polyphen2 Prediction') + 
    ylab('Count') +
    guides(fill=guide_legend(title = 'SIFT Prediction')) +
    ggtitle(study)
  
  return(list(score_vs_pph=p_score,
              raw_score_vs_pph=p_raw_score,
              pph_class_vs_score=p_class,
              pph_prediction_vs_score=p_prediction,
              pph_vs_sift_class=p_pph_vs_sift))
}

# Plot FoldX Scores
# Expects df with columns score, raw_score, foldx_PDB_ddG for PDB id
plot_foldx <- function(tbl, pdb_id, study=''){
  p_score <- ggplot(tbl, aes_string(x='score', y=str_c('foldx_', pdb_id, '_ddG'))) +
    xlab(MUT_SCORE_NAME) + 
    ylab('FoldX ddG') +
    ggtitle(study)
  
  if (dim(tbl)[1] < 1000){
    p_score <- p_score + geom_point() 
  } else {
    p_score <- p_score + geom_bin2d()
  }
  
  p_raw_score <- p_score + aes(x=raw_score) + xlab(RAW_MUT_SCORE_NAME)
  
  l <- list(p_score, p_raw_score)
  names(l) <- str_c('foldx_', pdb_id, c('_ddG_vs_score', '_ddG_vs_raw_score'))
  return(l)
}

# Plot additional FoldX plots for multiple variants
plot_foldx_multi <- function(tbl, pdb_id, study){
  l <- plot_foldx(tbl, pdb_id, study)
  
  p_vars_ddg <- ggplot(tbl, aes_string(x='count', y=str_c('foldx_', pdb_id, '_ddG'))) +
    geom_boxplot() +
    geom_smooth(method = 'lm', colour='red', aes(group=1)) +
    stat_summary(geom = 'text', fun.data = function(x){return(c(y = -3, label = length(x)))}) +
    xlab('Number of Mutations') + 
    ylab('ddG') +
    ggtitle(study)
  
  p_singles <- ggplot(filter(tbl, count==1),
                      aes_string(x='score', y=str_c('foldx_', pdb_id, '_ddG'))) +
    geom_point() +
    xlab(MUT_SCORE_NAME) +
    ylab('ddG') +
    ggtitle(study)
  
  l_multi <- list(p_vars_ddg, p_singles)
  names(l_multi) <- str_c('foldx_', pdb_id, c('_mutation_count_vs_ddG', '_ddG_vs_score_singles'))
  return(c(l_multi, l))
}

# Plot EVCouplings Scores
# Expects a df with variants, score, raw_score, evcoup_epistatic, evcoup_independent
plot_evcoup <- function(tbl, study=''){
  p_epistatic <- ggplot(tbl, aes(x=score, y=evcoup_epistatic)) +
    xlab(MUT_SCORE_NAME) +
    ylab('EVCouplings Epistatic Score') +
    ggtitle(study)
  
  p_ind <- ggplot(tbl, aes(x=score, y=evcoup_independent)) +
    xlab(MUT_SCORE_NAME) +
    ylab('EVCouplings Independent Score') +
    ggtitle(study)
  
  if (dim(tbl)[1] < 1000){
    p_epistatic <- p_epistatic + geom_point()
    p_ind <- p_ind + geom_point()
  } else {
    p_epistatic <- p_epistatic + geom_bin2d()
    p_ind <- p_ind + geom_bin2d()
  }
  
  p_raw_epistatic <- p_epistatic + aes(x=raw_score) + xlab(RAW_MUT_SCORE_NAME)
  p_raw_ind <- p_ind + aes(x=raw_score) + xlab(RAW_MUT_SCORE_NAME)
  
  return(list(evcoup_independent_vs_score=p_ind,
              evcoup_epistatic_vs_score=p_epistatic,
              evcoup_independent_vs_raw_score=p_raw_ind,
              evcoup_epistatic_vs_raw_score=p_raw_epistatic))
}

# Add plots only looking at single variants for EVCoup multi variant data
plot_evcoup_multi <- function(tbl, study){
  l <- plot_evcoup(tbl, study)
  l_sing <- plot_evcoup(filter(tbl, count == 1), study)
  names(l_sing) <- str_c(names(l_sing), '_single_vars')
  return(c(l_sing, l))
}

# Plot misc statistics
plot_misc <- function(tbl, study=''){
  tbl <- tbl %>%
    mutate(position=as.integer(str_sub(variants, 2, -2)),
           aa1 = str_sub(variants, 1, 1),
           aa2= str_sub(variants, -1))
  
  p_position <- ggplot(tbl, aes(group=position, x=position, y=score)) +
    geom_boxplot() +
    xlab('AA Position') + 
    ylab(MUT_SCORE_NAME) +
    ggtitle(study)
  
  return(list(aa_position_vs_score = p_position))
}

# Plot contingency tables
# Expects a tbl with two categorical cols and one score/proportion col + maybe a grouping col 
plot_contingency_table <- function(tbl, cat1, cat2, var, group=NULL,
                                   cat1_name=NULL, cat2_name=NULL, var_name=NULL){
  tbl$rounded_var <- round(tbl[[var]], 2)
  
  p <-ggplot(tbl, aes_string(x=cat1, y=var, fill=cat2)) +
    geom_col() + 
    geom_text(aes(label = rounded_var), 
              position = position_stack(vjust = 0.5)) + 
    xlab(ifelse(is.null(cat1_name), cat1, cat1_name)) +
    ylab(ifelse(is.null(var_name), var, var_name)) +
    guides(fill = guide_legend(title = ifelse(is.null(cat2_name), cat2, cat2_name)))
  
  if (!is.null(group)){
    p <- p + facet_wrap(as.formula(str_c('~', group)), scales = 'free')
  }
  
  return(p)
}

# Boxplots of predicted score vs experimental classification
plot_exp_pred_boxes <- function(tbl, y, x='exp_prediction', y_name=NULL, x_name='Experimental Prediction',
                                group='study', group_name = NULL){
  r <- range(tbl[y], na.rm = TRUE)
  y_n_lab <- min(tbl[y], na.rm = TRUE) - 0.1 * (r[2] - r[1])
  ymin <- min(tbl[y], na.rm = TRUE) - 0.1 * (r[2] - r[1])
  ymax <- max(tbl[y], na.rm = TRUE) + 0.2 * (r[2] - r[1])
  
  p_box_all <- ggplot(tbl, aes_string(x=x, y=y)) + 
    geom_boxplot() +
    xlab(x_name) +
    ylab(ifelse(is.null(y_name), y, y_name)) +
    stat_summary(geom = 'text', fun.data = function(x){return(c(y = y_n_lab,
                                                                label = length(x)))}) +
    stat_compare_means(comparisons = list(c(MUT_CATEGORIES$deleterious, MUT_CATEGORIES$neutral))) +
    stat_summary(geom = 'text', colour='red', fun.data = function(x){return(c(y = mean(x), label = signif(mean(x), 3)))}) +
    ggtitle('Overall') +
    ylim(ymin, ymax)
  
  if (is.null(group)){
    return(p_box_all)
  } else {
    p_box_facet <- p_box_all + 
      facet_wrap(as.formula(str_c('~', group))) +
      ggtitle(str_c('Per ', ifelse(is.null(group_name), group, group_name)))
    return(ggarrange(p_box_all, p_box_facet, ncol = 2, widths = c(1,2)))
  }
  
}

#### Misc Functions ####
format_pdb_variants <- function(x, pdb_offset=0){
  x <- str_split(x, ',')[[1]]
  str_sub(x, 2, 2) <- ''
  str_sub(x, 2, -2) <- as.character(as.integer(str_sub(x, 2, -2)) + pdb_offset)
  return(str_c(x, collapse = ','))
}

#### Convert experimental scores into predictions ####
# Interface wrapper to exchange current version
exp_mut_class <- function(score, study){
  return(exp_manual_thresh(score, study))
}

# Manual thresholds per study
MANUAL_THRESHOLDS <- c(araya_2012_yap1=-2, ashenberg_2016_nucleoprotein=-1, brenan_2016_erk2=-0.5, findlay_2018_brca1=-1,
                       firnberg_2014_tem1=-1.5, giacomelli_2018_tp53.p53_null_etoposide=-1, giacomelli_2018_tp53.p53_null_nutlin3=-0.75,
                       giacomelli_2018_tp53.p53_wt_nutlin3=-1, hietpas_2011_hsp90=-0.5, jiang_2013_hsp90=-0.5, kitzman_2015_gal4=-2.5,
                       lee_2018_ha=-3.75, melamed_2013_pab1=-1.25, melnikov_2014_aph3_ii.ami=-1, melnikov_2014_aph3_ii.g418=-0.75,
                       melnikov_2014_aph3_ii.kan=-2, melnikov_2014_aph3_ii.neo=-2, melnikov_2014_aph3_ii.paro=-2,
                       melnikov_2014_aph3_ii.ribo=-1.25, mishra_2016_hsp90=-0.5, olson_2014_proteing=-5, roscoe_2013_ubi4=-0.25,
                       roscoe_2014_ubi4.excess_e1=-1, roscoe_2014_ubi4.limiting_e1=-1.5, starita_2013_ube4b=-2.5,
                       starita_2015_brca1.bard1_ring_binding=-0.25, starita_2015_brca1.e3_activity=-1, wagenaar_2014_braf=0,
                       weile_2017_sumo1=-3, weile_2017_ube2i=-2)

exp_manual_thresh <- function(score, study){
  res <- rep(MUT_CATEGORIES$neutral, length(score))
  for (n in names(MANUAL_THRESHOLDS)){
    res[study == n & score < MANUAL_THRESHOLDS[n]] <- MUT_CATEGORIES$deleterious
  }
  return(res)
}

# Assume any deviation from neutral is bad
exp_not_neut <- function(x){
  res <- rep(MUT_CATEGORIES$neutral, length(x))
  res[abs(x) > 0.5] <- MUT_CATEGORIES$deleterious
  return(res)
}
