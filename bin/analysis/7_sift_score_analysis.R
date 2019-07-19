#!/usr/bin/env Rscript 
# Perform complimentary analysis on SIFT scores as DMS data

source('src/config.R')

sift <- readRDS('data/rdata/human_sift_reduced.RDS')

plots <- list()

#### Summary statistics for sift ####
plots$summary <- list()
plots$summary$score_distribution <- ggplot(gather(sift, key = 'aa', value = 'score', A:Y),
                                           aes(x=score + min(score[score > 0], na.rm = TRUE))) +
  facet_wrap(~aa) +
  geom_histogram() +
  scale_y_log10() +
  scale_x_log10()

plots$summary$protein_length <- group_by(sift, acc) %>%
  summarise(length = max(pos)) %>%
  ggplot(aes(x=length)) +
  geom_histogram()

plots$summary$aa_sub_distribution <- gather(sift, key = 'mut', value = 'score', A:Y) %>%
  select(acc, pos, wt, mut, score) %>%
  gather(key = 'type', value = 'aa', wt, mut) %>%
  drop_na(score) %>%
  mutate(aa = factor(aa, levels = names(AA_COLOURS))) %>%
  ggplot(aes(x=aa, fill=aa)) +
    geom_bar() +
    scale_fill_manual(values = AA_COLOURS) +
    facet_wrap(~type)

plots$summary$aa_wt_distribution <- mutate(sift, wt = factor(wt, levels = names(AA_COLOURS))) %>%
  ggplot(aes(x=wt, fill=wt)) +
  geom_bar() +
  scale_fill_manual(values = AA_COLOURS)

plots$summary$metric_distribution <- ggpairs(sift, columns = c('median_ic', 'n_aa', 'n_seq'))
########

# Save plots
save_plot_list(plots, root='figures/7_sift_score_analysis')
