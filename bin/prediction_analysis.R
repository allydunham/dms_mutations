#!/usr/bin/env Rscript 
# Script to load processed deep mutagenesis data for analysis

source('bin/libraries.R')
source('bin/dm_functions.R')
source('bin/variant_functions.R')

# Import data
formatted_deep_data <-readRDS('data/formatted_deep_mut_data.RDS')
deep_variant_data <- readRDS('data/variant_data.RDS')

# Basic plots
deep_variant_plots <- sapply(deep_variant_data, plot_predictions, simplify = FALSE)

# General analysis

### Save plots ####
for (study_name in names(deep_variant_plots)){
  print(str_c('Writing figures for ', study_name))
  fig_root <- str_c('figures/variant_analysis/', study_name)
  dir.create(fig_root, showWarnings = FALSE, recursive = TRUE)
  for (plot_name in names(deep_variant_plots[[study_name]])){
    plot_path <- str_c(fig_root, '/', plot_name, '.pdf')
    # Only write figures that don't exist, delete figs to regenerate
    if (!file.exists(plot_path)){
      ggsave(plot_path, deep_variant_plots[[study_name]][[plot_name]], width = 7, height = 5)
    }
  }
}
