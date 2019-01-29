#!/usr/bin/env Rscript 
# Script containing functions for loading and analysis of data from Envision, SIFT, FoldX, Polyphen2 & EVCouplings

# Constants
MUT_SCORE_NAME <- 'Standardised Mutagenesis Score'
RAW_MUT_SCORE_NAME <- 'Raw Mutagenesis Score'

#### Plotting function ####
## Plot deep mutagenesis data and variant effect predictions
# Generic
plot_predictions <- function(x, ...){
  UseMethod('plot_predictions', x)
}

# Single Variants
plot_predictions.single_variant <- function(x){
  study <- str_c(x$dm$authour, ' ', x$dm$year, ': ', x$dm$gene_name)
  
  plots <- c(plot_sift(x$single_variants, study),
             sapply(names(x$foldx), function(id){plot_foldx(x$single_variants, id, study)}, simplify = FALSE) %>%
               unname(.) %>%
               unlist(., recursive = FALSE),
             plot_pph(x$single_variants, study),
             plot_envision(x$single_variants, study),
             plot_evcoup(x$single_variants, study))
  return(plots)
}

# Multi Variants
plot_predictions.multi_variant <- function(x){
  study <- str_c(x$dm$authour, ' ', x$dm$year, ': ', x$dm$gene_name)

  plots <- c(plot_sift(x$single_variants, study),
             sapply(names(x$foldx), function(id){plot_foldx(x$multi_variants, id, study)}, simplify = FALSE) %>%
               unname(.) %>%
               unlist(., recursive = FALSE),
             plot_pph(x$single_variants, study),
             plot_envision(x$single_variants, study),
             plot_evcoup(x$multi_variants, study))
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
    p_score <- p_score + geom_density2d()
  }
  
  p_raw_score <- p_score + aes(x=raw_score) + xlab(RAW_MUT_SCORE_NAME)
    
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
              sift_prediction_vs_score=p_prediction))
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
# Expects df with columns score, raw_score, pph2_prediction, pph2_class, pph2_prob pph2_FPR pph2_TPR pph2_FDR
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
    aes(x=pph2_prediction)
    stat_compare_means(comparisons = list(c("probably damaging", "possibly damaging", "benign"))) +
    xlab('PolyPhen2 Prediction')
  
    
  p_pph_vs_sift <- ggplot(drop_na(tbl, pph2_prediction, sift_prediction),
                          aes(x=pph2_prediction, y=..count.., fill=sift_prediction)) +
    geom_bar(position = 'dodge') +
    xlab('Polyphen2 Prediction') + 
    ylab('Count') +
    guides(fill=guide_legend(title = 'SIFT Prediction'))
  
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
  names(l) <- c(str_c('foldx_', pdb_id, '_ddG_vs_score'), str_c('foldx_', pdb_id, '_ddG_vs_raw_score'))
  return(l)
}

# Plot EVCouplings Scores
# Expects a df with variants, score, raw_score, evcoup_epistatic, evcoup_independent
plot_evcoup <- function(tbl, study=''){
  p_epistatic <- ggplot(tbl, aes(x=score, y=evcoup_epistatic)) +
    xlab(MUT_SCORE_NAME) +
    ylab('EVCouplings Epistatic Score')
  
  p_ind <- ggplot(tbl, aes(x=score, y=evcoup_independent)) +
    xlab(MUT_SCORE_NAME) +
    ylab('EVCouplings Independent Score')
  
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

# Plot misc statistics
plot_misc <- function(tbl, study=''){
  tbl <- tbl %>%
    mutate(pos=as.integer(str_sub(variants, 2, -2)),
           aa1 = str_sub(variants, 1, 1),
           aa2= str_sub(variants, -1))
  
  p_position <- ggplot(t, aes(group=pos, x=pos, y=score)) +
    geom_boxplot() +
    xlab('AA Position') + 
    ylab(MUT_SCORE_NAME) +
    ggtitle(study)
}

#### Misc Functions ####
remove_pdb_chains <- function(x){
  x <- str_split(x, ',')[[1]]
  str_sub(x, 2, 2) <- ''
  return(str_c(x, collapse = ','))
}
