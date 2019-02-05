#!/usr/bin/env Rscript 
# Script to load processed deep mutagenesis data for analysis

source('bin/libraries.R')
source('bin/dm_functions.R')
source('bin/variant_functions.R')

# Import data
formatted_deep_data <-readRDS('data/formatted_deep_mut_data.RDS')
deep_variant_data <- readRDS('data/variant_data.RDS')

# Basic plots per study
deep_variant_plots <- sapply(deep_variant_data, plot_predictions, simplify = FALSE)

### General plots across studies ###
deep_variant_plots$combined_plots <- list()

#### DM Data ####
all_dm <- bind_rows(lapply(deep_variant_data, function(x){x$dm$variant_data}), .id='study') %>%
  select(study, variants, score, raw_score) %>%
  mutate(variants = str_replace_all(variants, 'p.', ''))

deep_variant_plots$combined_plots$dm_hists <- ggplot(all_dm, aes(x=score)) + 
  geom_histogram() +
  facet_wrap(~study, scales = 'free')

#### Envision ####
select_envision <- function(x){
  if ('envision_prediction' %in% names(x$single_variants)){
    return(select(x$single_variants, variants, score, envision_prediction, log2_envision_prediction))
  } else {
    return(NULL)
  }
}

envision <- bind_rows(lapply(deep_variant_data, select_envision), .id='study')

envision_cor_test <- group_by(envision, study) %>%
  do(tidy(cor.test(.$score, .$log2_envision_prediction)))
max_abs_t <- max(abs(envision_cor_test$statistic))
  
deep_variant_plots$combined_plots$envision_correlation <- ggplot(envision_cor_test, aes(y=estimate, x=study, fill=statistic)) + 
  geom_col() +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high), width=0.5) +
  ggtitle('Correlation between log2(Envision predictions) and mutagenesis scores') +
  xlab('') +
  ylab('Pearson Correlation Coefficient') +
  scale_fill_gradientn(colours = c('red', 'white', 'blue'), limits = c(-max_abs_t, max_abs_t)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

#### SIFT ####
select_sift <- function(x){
  if ('sift_prediction' %in% names(x$single_variants)){
    return(select(x$single_variants, variants, score, sift_prediction, sift_score, sift_median))
  } else {
    return(NULL)
  }
}

sift <- bind_rows(lapply(deep_variant_data, select_sift), .id='study') %>%
  mutate(exp_prediction = score < -1,
         sift_prediction = str_to_lower(sift_prediction))

deep_variant_plots$combined_plots$sift_prediction_accuracy <- ggplot(drop_na(sift, exp_prediction),
                                                                     aes(x=sift_prediction, fill=exp_prediction)) +
  facet_wrap(~study, scales = 'free') +
  geom_bar(position = 'dodge') + 
  guides(fill=guide_legend(title = 'Score < -0.5'))

#### FoldX ####

#### Polyphen2 ####

#### EVCouplings ####

#### Save plots #####
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
