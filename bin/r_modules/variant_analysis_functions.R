#!/usr/bin/env Rscript 
# Script containing functions for loading and analysis of data from Envision, SIFT, FoldX, Polyphen2 & EVCouplings

#### General Functions ####
# Function generating a per position mutational profile for a given study (i.e. enrichment ratio of all substitutions at a position)
make_var_matrix <- function(x, score='score'){
  variants <- select(x$single_variants, variants, score=!!score) %>%
    mutate(wt = str_sub(variants, end=1),
           mut = str_sub(variants, start=-1),
           pos = as.integer(str_sub(variants, start=2, end=-2))) %>%
    select(-variants)
  
  if ('=' %in% variants$mut){
    variants$mut[variants$mut=='='] <- variants$wt[variants$mut=='=']
  }
  
  variants <- spread(variants, key = 'mut', value = 'score') %>%
    arrange(pos) %>%
    select(pos, wt, A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y)
  
  return(variants)
}

# Impute missing values in variant profiles matrix
# impute function is applied to each subsititution type (e.g. A->C) accross all 
# positions in each study to build a per study profile
impute_variant_profiles <- function(variant_matrix, background_matrix=NULL, impute_function=median){
  # Per study/per AA median mutational profiles
  if (is.null(background_matrix)){
    background_matrix <- variant_matrix
  }
  per_study_mean_profiles <- select(background_matrix, study, pos, wt, A:Y) %>%
    group_by(study, wt) %>%
    summarise_at(.vars = vars(-pos), .funs = impute_function, na.rm=TRUE) %>%
    replace(is.na(.), 0)
  
  # Matrix of mutational profiles with missing values imputed from the per study/AA medians
  imputed_matrix <- gather(variant_matrix, key='mut', value='score', A:Y) %>%
    left_join(., gather(per_study_mean_profiles, key='mut', value='imp_score', A:Y), by = c("study", 'mut', 'wt'))
  
  imputed_matrix$score[is.na(imputed_matrix$score)] <- imputed_matrix$imp_score[is.na(imputed_matrix$score)]
  
  imputed_matrix <- imputed_matrix %>%
    select(-imp_score) %>%
    spread(key=mut, value=score)
  
  return(imputed_matrix)
}

#### Plots ####
# Plot study histograms
plot_study_histogram <- function(tbl, thresh_tbl=NULL, x='score', fill='authour', facet='~study', thresh='thresh'){
  p <- ggplot(tbl, aes_string(x=x, fill='authour')) + 
    geom_histogram() +
    facet_wrap(facet, scales = 'free') +
    xlab(MUT_SCORE_NAME) +
    ylab('Count')
  
  if(!is.null(thresh_tbl)){
    p <- p + geom_vline(aes_string(xintercept=thresh), data = thresh_tbl, colour='red')
  }
  
  return(p)
}

# Plot factor frequencies
plot_factor_density <- function(tbl, facet, x='score', col='authour'){
  p <- ggplot(tbl, aes_string(x=x, colour=col)) + 
    geom_density(trim=TRUE) +
    facet_wrap(facet, scales = 'free') +
    xlab(MUT_SCORE_NAME) +
    ylab('Density') +
    scale_color_discrete()
  
  return(p)
}

# Plot distribution of enrichment score along length of proteins
# Expects a tibble with columns pos and score, 
# with optional column study from which title is derived if unique and not specified
plot_score_distribution_per_position <- function(tbl, title=NULL){
  if (is.null(title)){
    if ('study' %in% names(tbl) & col_unique_counts(tbl)['study'] == 1){
      title <- unique(tbl$study)
    } else {
      title = ''
    }
  }
  
  p <- ggplot(tbl, aes(group=pos, x=pos, y=score)) +
    geom_boxplot() +
    xlab('AA Position') + 
    ylab(MUT_SCORE_NAME) +
    ggtitle(title)
  
  return(p)
}

#### Predicting from experimental scores ####
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


#### Misc helper functions ####
col_unique_counts <- function(tbl){
  return(apply(tbl, 2, function(x){length(unique(x))}))
}
