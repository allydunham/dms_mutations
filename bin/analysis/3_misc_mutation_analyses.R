#!/usr/bin/env Rscript 
# Script containing miscelaneous analyses of deep mutagenesis studies too small for their own script

source('src/config.R')

variant_matrices <- readRDS('data/rdata/all_study_position_matrices.RDS')

plots <- list()

#### Mutation severity vs Blosum ####
data(BLOSUM62)
sub_profile_scores <- variant_matrices$norm_all_variants %>%
  select(wt, A:Y) %>%
  group_by(wt) %>%
  summarise_all(.funs = mean, na.rm=TRUE) %>%
  gather(key = 'mut', value = 'er', -wt) %>%
  left_join(., BLOSUM62 %>%
              as_tibble(rownames = 'wt') %>%
              gather(key='mut', value = 'blosum62', -wt),
            by=c('wt', 'mut'))

sub_scores <- variant_matrices$norm_all_variants %>%
  select(study, pos, wt, A:Y) %>%
  gather(key = 'mut', value = 'er', -wt, -study, -pos) %>%
  left_join(., BLOSUM62 %>%
              as_tibble(rownames = 'wt') %>%
              gather(key='mut', value = 'blosum62', -wt),
            by=c('wt', 'mut'))

plots$aa_profile_heatmap <- ggplot(sub_profile_scores, aes(x=wt, y=mut, fill=er)) + 
  geom_tile(colour='white') + 
  scale_fill_gradient2() +
  theme(axis.ticks = element_blank(), panel.background = element_blank())

plots$er_vs_blosum <- ggplot(sub_scores, aes(x=blosum62, y=er, group=blosum62)) + geom_boxplot()
plots$blosum$profile_er_vs_blosum <- ggplot(sub_profile_scores, aes(x=blosum62, y=er, group=blosum62)) + 
  geom_boxplot()
########

save_plot_list(plots, root='figures/3_misc_mutation_analyses/')
