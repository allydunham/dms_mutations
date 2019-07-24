#!/usr/bin/env Rscript 
# Functions for analysis of SIFT Scores in the manner of dms data

#### Summary distributions ####
plot_sift_score_summary <- function(tbl){
  plots <- list()
  plots$score_distribution <- ggplot(gather(tbl, key = 'aa', value = 'score', A:Y),
                                             aes(x=score + min(score[score > 0], na.rm = TRUE))) +
    facet_wrap(~aa) +
    geom_histogram() +
    scale_x_continuous(labels = make_log_labeler(base=10, force = 'exp')) +
    scale_y_log10()
  
  plots$per_sub_type_score_distribution <- labeled_ggplot(
    p = ggplot(gather(tbl, key = 'mut', value = 'score', A:Y), aes(x=score)) +
      geom_histogram() +
      facet_grid(rows = vars(mut), cols = vars(wt)) + 
      ggtitle(expression('Distribution of log'[10]*'(SIFT +'~epsilon*') for each substitution (wt: cols, mut: rows)')) +
      scale_x_continuous(labels = make_log_labeler(base = 10, force = 'exp')) +
      scale_y_log10() +
      theme_pubclean() +
      theme(strip.background = element_blank(),
            legend.position = 'right'),
    units = 'cm', height = 44, width = 44)
  
  plots$protein_length <- group_by(tbl, acc) %>%
    summarise(length = max(pos)) %>%
    ggplot(aes(x=length)) +
    geom_histogram()
  
  plots$aa_sub_distribution <- gather(tbl, key = 'mut', value = 'score', A:Y) %>%
    select(acc, pos, wt, mut, score) %>%
    gather(key = 'type', value = 'aa', wt, mut) %>%
    drop_na(score) %>%
    mutate(aa = factor(aa, levels = names(AA_COLOURS))) %>%
    ggplot(aes(x=aa, fill=aa)) +
    geom_bar() +
    scale_fill_manual(values = AA_COLOURS) +
    facet_wrap(~type)
  
  plots$metric_distribution <- ggpairs(tbl, columns = c('median_ic', 'n_aa', 'n_seq'))
  
  return(plots)
}
########