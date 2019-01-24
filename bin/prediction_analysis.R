#!/usr/bin/env Rscript 
# Script to load processed deep mutagenesis data for analysis

library(tidyverse)

source('bin/dm_functions.R')
source('bin/variant_functions.R')

# Import data
formatted_deep_data <-readRDS('data/formatted_deep_mut_data.RDS')
deep_variant_data <- readRDS('data/variant_data.RDS')

# Basic plots
deep_variant_plots <- sapply(deep_variant_data, plot_predictions, simplify = FALSE)

#### Araya et al. 2012 ####
df <- deep_variant_data$araya_2012_hYAP65$multi_variants %>%
  mutate(count = sapply(variants, function(x){dim(str_split(x, ',', simplify = TRUE))[2]}))

deep_variant_plots$araya_2012_hYAP65$foldx_4REX_ddG_vs_mut_count <- ggplot(df, aes(x=factor(count), y=foldx_4REX_ddG)) +
  geom_boxplot() +
  geom_smooth(method = 'lm', colour='red', aes(group=1)) +
  stat_summary(geom = 'text', fun.data = function(x){return(c(y = -3, label = length(x)))}) +
  xlab('Number of Mutations') + 
  ylab('ddG') +
  ggtitle('Araya et al. 2012: hYAP65 Variants')

deep_variant_plots$araya_2012_hYAP65$foldx_4REX_ddG_vs_score_single_vars <- ggplot(filter(df, count==1),
                                                                           aes(x=score, y=foldx_4REX_ddG)) +
  geom_point() +
  xlab(MUT_SCORE_NAME) +
  ylab('ddG') +
  ggtitle('Araya et al. 2012: hYAP65')

#### Save plots ####
for (study_name in names(deep_variant_plots)){
  fig_root <- str_c('figures/variant_analysis/', study_name)
  dir.create(fig_root, showWarnings = FALSE, recursive = TRUE)
  for (plot_name in names(deep_variant_plots[[study_name]])){
    ggsave(str_c(fig_root, '/', plot_name, '.pdf'), deep_variant_plots[[study_name]][[plot_name]], width = 7, height = 5)
  }
}