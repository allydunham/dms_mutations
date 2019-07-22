#!/usr/bin/env Rscript 
# Functions for analysis of SIFT Scores in the manner of dms data

#### Summary distributions ####
plot_sift_score_summary <- function(tbl){
  plots <- list()
  plots$score_distribution <- ggplot(gather(tbl, key = 'aa', value = 'score', A:Y),
                                             aes(x=score + min(score[score > 0], na.rm = TRUE))) +
    facet_wrap(~aa) +
    geom_histogram() +
    scale_x_continuous(trans = scales::pseudo_log_trans(sigma=10^-6, base=10),
                       breaks = c(0, 10^-4, 10^-2, 1),
                       labels = c('0', expr(10^-4), expr(10^-2), '1')) +
    scale_x_log10()
  
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
  
  plots$aa_wt_distribution <- mutate(tbl, wt = factor(wt, levels = names(AA_COLOURS))) %>%
    ggplot(aes(x=wt, fill=wt)) +
    geom_bar() +  
    scale_fill_manual(values = AA_COLOURS)
  
  plots$metric_distribution <- ggpairs(tbl, columns = c('median_ic', 'n_aa', 'n_seq'))
  
  return(plots)
}
########