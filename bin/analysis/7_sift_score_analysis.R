#!/usr/bin/env Rscript 
# Perform complimentary analysis on SIFT scores as DMS data

source('src/config.R')
source('src/analysis/sift_score_analysis.R')
source('src/analysis/position_profile_clustering.R')

sift <- readRDS('data/rdata/human_sift_reduced.RDS')
sift_no_missing <- readRDS('data/rdata/human_sift_no_missing_reduced.RDS')

sum_missing <- tibble_to_matrix(sift, A:Y) %>% is.na() %>% sum()
total_scores <- dim(sift)[1] * 20
prop_missing <- sum_missing / total_scores

message(str_c(signif(prop_missing, 4), '% of SIFT scores missing (', sum_missing, ' / ', total_scores,')'))

plots <- list()

#### Summary statistics for sift ####
plots$summary <- plot_sift_score_summary(sift)
p <- plot_sift_score_summary(sift_no_missing)
names(p) <- str_c(names(p), '_no_missing')
plots$summary <- c(plots$summary, p)

# Compare positions where information is missing or not
plots$summary$missing_comparison <- list()
comb <- bind_rows(yes=gather(sift, key = 'mut', value = 'score', A:Y),
                  no=gather(sift_no_missing, key = 'mut', value = 'score', A:Y),
                  .id = 'missing')

plots$summary$missing_comparison$per_sub_type_score_distribution <- labeled_ggplot(
  p = ggplot(comb, aes(x=score, fill=missing)) +
    geom_histogram(alpha = 0.5) +
    facet_grid(rows = vars(mut), cols = vars(wt)) + 
    scale_fill_manual(values = c(yes='red', no='black')) +
    scale_x_continuous(trans = 'pseudo_log') +
    scale_y_log10() +
    theme_pubclean() +
    theme(strip.background = element_blank(),
          legend.position = 'right'),
  units = 'cm', height = 44, width = 44)

plots$summary$missing_comparison$score_distribution <- labeled_ggplot(
  p = ggplot(gather(comb, key = 'type', value = 'aa', wt, mut) %>% mutate(type = factor(type, levels = c('wt', 'mut'))),
             aes(x=score, fill=missing)) +
    geom_histogram(alpha = 0.5) +
    facet_grid(rows = vars(type), cols = vars(aa)) + 
    scale_fill_manual(values = c(yes='red', no='black')) +
    scale_x_continuous(trans = scales::pseudo_log_trans(sigma=10^-6, base=10),
                       breaks = c(0, 10^-4, 10^-2, 1),
                       labels = c('0', expr(10^-4), expr(10^-2), '1')) +
    scale_y_log10() +
    theme_pubclean() +
    theme(strip.background = element_blank(),
          legend.position = 'right',
          axis.text.x = element_text(angle = 45)),
  units = 'cm', height = 10, width = 44)

plots$summary$missing_comparison$metric_distribution <- labeled_ggplot(
  p = ggplot(gather(comb, key = 'metric', value = 'value', median_ic, n_aa, n_seq),
             aes(x=value, fill=missing)) +
    geom_histogram(alpha = 0.5) +
    facet_wrap(~metric, nrow = 3, scales = 'free') + 
    scale_fill_manual(values = c(yes='red', no='black')) +
    scale_y_log10() +
    theme_pubclean() +
    theme(strip.background = element_blank(),
          legend.position = 'right'),
  units = 'cm', height = 30, width = 15)

plots$summary$missing_comparison$aa_freq_distribution <- labeled_ggplot(
  p = ggplot(gather(comb, key = 'type', value = 'aa', wt, mut) %>%
               drop_na(score) %>%
               mutate(type = factor(type, levels = c('wt', 'mut'))),
             aes(x=aa, fill=missing, y=..prop.., group=missing)) +
    geom_bar(position = 'dodge') + 
    facet_wrap(~type, nrow = 2) +
    scale_fill_manual(values = c(yes='red', no='black')) +
    theme_pubclean() +
    theme(strip.background = element_blank(),
          legend.position = 'right'),
  units = 'cm', height = 20, width = 20)
########

#### Cluster AA profiles ####
h <- 8
hclust_clusters <- group_by(sift_no_missing, wt) %>%
  do(hclust = make_hclust_clusters(., A:Y, h = h))

hclust_tbl <- map_dfr(hclust_clusters$hclust, .f = ~ .[[1]]) %>%
  mutate(cluster = str_c(wt, '_', cluster))

########

# Save plots
save_plot_list(plots, root='figures/7_sift_score_analysis')
