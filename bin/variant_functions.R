#!/usr/bin/env Rscript 
# Script containing functions for loading and analysis of data from Envision, SIFT, FoldX, Polyphen2 & EVCouplings
library(tidyverse)
library(ggpubr)

MUT_SCORE_NAME <- 'Standardised Mutagenesis Score'
RAW_MUT_SCORE_NAME <- 'Raw Mutagenesis Score'

## Plot deep mutagenesis data and variant effect predictions
# Generic
plot_predictions <- function(x, ...){
  UseMethod('plot_predictions', x)
}

# Single Variants
plot_predictions.single_variant <- function(x){
  study <- str_c(x$dm$authour, ' ', x$dm$year, ': ', x$dm$gene_name)
  
  plots <- c(plot_sift(x$single_variants, study),
             plot_foldx(x$single_variants, study),
             plot_pph(x$single_variants, study),
             plot_envision(x$single_variants, study),
             plot_evcoup(x$single_variants, study))
  return(plots)
}

# Multi Variants
plot_predictions.multi_variant <- function(x){
  study <- str_c(x$dm$authour, ' ', x$dm$year, ': ', x$dm$gene_name)
  
  plots <- c(plot_sift(x$single_variants, study),
             plot_foldx(x$multi_variants, study),
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
# Expects df with columns variant, score, raw_score, 
plot_envision <- function(tbl, study=''){
  p_score <- ggplot(tbl, aes(x=score, y=envision_prediction)) +
    xlab(MUT_SCORE_NAME) + 
    ylab('Envision Prediction') +
    ggtitle(study)
  
  if (dim(tbl)[1] < 1000){
    p_score <- p_score + geom_point() 
  } else {
    p_score <- p_score + geom_density2d()
  }
  
  p_raw_score <- p_score + aes(x=raw_score) + xlab(RAW_MUT_SCORE_NAME)
  
  return(list(score_vs_envision=p_score,
              raw_score_vs_envision=p_raw_score))
}

# Plot Polyphen2 Scores
# Expects df with columns variant, score, raw_score, pph2_class, pph2_prob pph2_FPR pph2_TPR pph2_FDR
plot_pph <- function(tbl, study=''){
  p_score <- ggplot(tbl, aes(x=score, y=pph_prob)) +
    xlab(MUT_SCORE_NAME) + 
    ylab('PolyPhen2 Probability') +
    ggtitle(study)
  
  if (dim(tbl)[1] < 1000){
    p_score <- p_score + geom_point() 
  } else {
    p_score <- p_score + geom_density2d()
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
  
  return(list(score_vs_pph=p_score,
              raw_score_vs_pph=p_raw_score,
              pph_class_vs_score=p_class))
}

# Plot FoldX Scores
plot_foldx <- function(tbl, study=''){
  
}

# Plot EVCouplings Scores
# Expects a df with variant, score, raw_score, evcoup_epistatic, evcoup_independent
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
}

# Plot misc statistics
plot_misc <- function(tbl, study=''){
  tbl <- tbl %>%
    mutate(pos=as.integer(str_sub(variant, 2, -2)),
           aa1 = str_sub(variant, 1, 1),
           aa2= str_sub(variant, -1))
  
  p_position <- ggplot(t, aes(group=pos, x=pos, y=score)) +
    geom_boxplot() +
    xlab('AA Position') + 
    ylab(MUT_SCORE_NAME) +
    ggtitle(study)
}