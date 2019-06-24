#!/usr/bin/env Rscript 
# Functions for plotting summary plots of enrichment scores gathered in the studies

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
plot_factor_density <- function(tbl, facet, x='norm_score', col='authour'){
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