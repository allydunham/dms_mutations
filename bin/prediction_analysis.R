#!/usr/bin/env Rscript 
# Script to load processed deep mutagenesis data for analysis

source('bin/config.R')

# Import data
formatted_deep_data <-readRDS('data/formatted_deep_mut_data.RDS')
deep_variant_data <- readRDS('data/variant_data.RDS')

# Basic plots per study
deep_variant_plots <- sapply(deep_variant_data, plot_predictions, simplify = FALSE)

### General plots across studies ###
deep_variant_plots$all_studies <- list()

#### DM Data ####
all_dm <- bind_rows(lapply(deep_variant_data, function(x){x$dm$variant_data}), .id='study') %>%
  select(study, variants, score, raw_score) %>%
  mutate(variants = str_replace_all(variants, 'p.', ''))

deep_variant_plots$all_studies$dm_hists <- labeled_ggplot(p = ggplot(all_dm, aes(x=score)) + 
                                                            geom_histogram() +
                                                            facet_wrap(~study, scales = 'free') +
                                                            xlab(MUT_SCORE_NAME) +
                                                            ylab('Count'),
                                                          width = 12, height = 9)

#### Envision ####
select_envision <- function(x){
  if ('envision_prediction' %in% names(x$single_variants)){
    return(select(x$single_variants, variants, score, envision_prediction, log2_envision_prediction))
  } else {
    return(NULL)
  }
}

envision <- bind_rows(lapply(deep_variant_data, select_envision), .id='study') %>%
  mutate(exp_prediction = predict_exp_function(score))

deep_variant_plots$all_studies$envision_experimental_boxplot <- labeled_ggplot(plot_exp_pred_boxes(envision,
                                                                                                   y = 'envision_prediction',
                                                                                                   y_name = 'Envision Prediction'),
                                                                               width = 12,
                                                                               height = 8)

envision_cor_test <- group_by(envision, study) %>%
  do(tidy(cor.test(.$score, .$log2_envision_prediction)))
max_abs_t <- max(abs(envision_cor_test$statistic))
  
deep_variant_plots$all_studies$envision_correlation <- ggplot(envision_cor_test, aes(y=estimate, x=study, fill=statistic)) + 
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
  mutate(exp_prediction = predict_exp_function(score),
         sift_prediction = str_to_lower(gsub('TOLERATED', MUT_CATEGORIES$neutral, sift_prediction)))

deep_variant_plots$all_studies$sift_experimental_boxplot <- labeled_ggplot(plot_exp_pred_boxes(sift,
                                                                                               y = 'sift_score',
                                                                                               y_name = 'SIFT Score'),
                                                                           width = 12,
                                                                           height = 8)

sift_summary <- group_by(sift, study, sift_prediction) %>%
  drop_na(sift_prediction, exp_prediction) %>%
  summarise(num_vars = n(),
            deleterious = sum(exp_prediction == MUT_CATEGORIES$deleterious),
            neutral = sum(exp_prediction == MUT_CATEGORIES$neutral)) %>%
  gather(key = 'exp_prediction', value = 'count', deleterious, neutral) %>%
  group_by(study, sift_prediction) %>%
  mutate(prop = count/sum(count)) %>%
  arrange(study, sift_prediction)

deep_variant_plots$all_studies$sift_prediction_counts <- labeled_ggplot(p = plot_contingency_table(sift_summary,
                                                                                                     cat1 = 'sift_prediction',
                                                                                                     cat2 = 'exp_prediction',
                                                                                                     var = 'count', group = 'study',
                                                                                                     cat1_name = 'SIFT Prediction',
                                                                                                     cat2_name = 'Exp. Prediction',
                                                                                                     var_name = 'Count'),
                                                                          width = 12, height = 9)

deep_variant_plots$all_studies$sift_prediction_accuracy <- labeled_ggplot(p = plot_contingency_table(sift_summary,
                                                                                                     cat1 = 'sift_prediction',
                                                                                                     cat2 = 'exp_prediction',
                                                                                                     var = 'prop', group = 'study',
                                                                                                     cat1_name = 'SIFT Prediction',
                                                                                                     cat2_name = 'Exp. Prediction',
                                                                                                     var_name = 'Proportion'),
                                                                          width = 12, height = 9)

#### FoldX ####
select_foldx <- function(x){
  if (any(grepl('foldx_', names(x$multi_variants)))){
    tbl <- x$multi_variants
  } else if (any(grepl('foldx_', names(x$single_variants)))){
    tbl <- x$single_variants
  } else {
    return(NULL)
  }
  tbl <- select(tbl, variants, score, starts_with('foldx_')) %>%
    rename_at(vars(contains('foldx_')), .funs = funs(gsub('foldx_', '', .))) %>%
    gather(key = 'k', value = 'v', -variants, -score) %>%
    separate(k, into = c('pdb_id', 'k'), sep = '_') %>%
    spread(key = k, value = v)
  return(tbl)
}

foldx <- bind_rows(lapply(deep_variant_data, select_foldx), .id='study') %>%
  mutate(count = factor(sapply(variants, function(x){dim(str_split(x, ',', simplify = TRUE))[2]})),
         single = count == 1,
         exp_prediction = predict_exp_function(score))

deep_variant_plots$all_studies$foldx_experimental_boxplot <- labeled_ggplot(plot_exp_pred_boxes(foldx,
                                                                                                y = 'ddG'),
                                                                            width = 12,
                                                                            height = 8)

deep_variant_plots$all_studies$foldx_experimental_boxplot_single <- labeled_ggplot(plot_exp_pred_boxes(filter(foldx, single),
                                                                                                       y = 'ddG'),
                                                                                   width = 12,
                                                                                   height = 8)

# Correlations
foldx_cor_test <- group_by(foldx, study, pdb_id, single) %>%
  drop_na(ddG) %>%
  do(tidy(cor.test(abs(.$score), .$ddG))) %>%
  mutate(multi = ifelse(single, '', 'multi')) %>%
  unite(name, study, pdb_id, multi)

max_abs_t <- max(abs(foldx_cor_test$statistic))
deep_variant_plots$all_studies$foldx_correlation <- ggplot(foldx_cor_test, aes(y=estimate, x=name, fill=statistic)) + 
  geom_col() +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high), width=0.5) +
  ggtitle('Correlation between FoldX ddG and abs(experimental scores)') +
  xlab('') +
  ylab('Pearson Correlation Coefficient') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# Predictions
foldx_summary <- mutate(foldx, foldx_prediction = ifelse(abs(ddG) > 2, MUT_CATEGORIES$deleterious, MUT_CATEGORIES$neutral)) %>%
  group_by(study, foldx_prediction)  %>%
  summarise(num_vars = n(),
            deleterious = sum(exp_prediction == MUT_CATEGORIES$deleterious),
            neutral = sum(exp_prediction == MUT_CATEGORIES$neutral)) %>%
  gather(key = 'exp_prediction', value = 'count', deleterious, neutral) %>%
  group_by(study, foldx_prediction) %>%
  mutate(prop = count/sum(count)) %>%
  arrange(study, foldx_prediction)

deep_variant_plots$all_studies$foldx_prediction_counts <- labeled_ggplot(p = plot_contingency_table(foldx_summary,
                                                                                                    cat1 = 'foldx_prediction',
                                                                                                    cat2 = 'exp_prediction',
                                                                                                    var = 'count', group = 'study',
                                                                                                    cat1_name = 'abs(ddG) > 2',
                                                                                                    cat2_name = 'Exp. Prediction',
                                                                                                    var_name = 'Count'),
                                                                         width = 12, height = 9)

deep_variant_plots$all_studies$foldx_prediction_accuracy <- labeled_ggplot(p = plot_contingency_table(foldx_summary,
                                                                                                      cat1 = 'foldx_prediction',
                                                                                                      cat2 = 'exp_prediction',
                                                                                                      var = 'prop', group = 'study',
                                                                                                      cat1_name = 'abs(ddG) > 2',
                                                                                                      cat2_name = 'Exp. Prediction',
                                                                                                      var_name = 'Proportion'),
                                                                           width = 12, height = 9)


#### Polyphen2 ####
select_pph <- function(x){
  if ('pph2_class' %in% names(x$single_variants)){
    return(select(x$single_variants, variants, score, pph2_prediction, pph2_class, pph2_prob, pph2_FPR, pph2_TPR, pph2_FDR))
  } else {
    return(NULL)
  }
}

pph <- bind_rows(lapply(deep_variant_data, select_pph), .id='study') %>%
  mutate(exp_prediction = predict_exp_function(score))

deep_variant_plots$all_studies$pph_experimental_boxplot <- labeled_ggplot(plot_exp_pred_boxes(pph,
                                                                                              y = 'pph2_prob',
                                                                                              y_name = 'Polyphen2 Probability'),
                                                                          width = 12,
                                                                          height = 8)


pph_summary <- group_by(pph, study, pph2_class) %>%
  drop_na(pph2_class, exp_prediction) %>%
  summarise(num_vars = n(),
            deleterious = sum(exp_prediction == MUT_CATEGORIES$deleterious),
            neutral = sum(exp_prediction == MUT_CATEGORIES$neutral)) %>%
  gather(key = 'exp_prediction', value = 'count', deleterious, neutral) %>%
  group_by(study, pph2_class) %>%
  mutate(prop = count/sum(count)) %>%
  arrange(study, pph2_class)

deep_variant_plots$all_studies$pph_prediction_counts <- labeled_ggplot(p = plot_contingency_table(pph_summary,
                                                                                                   cat1 = 'pph2_class',
                                                                                                   cat2 = 'exp_prediction',
                                                                                                   var = 'count', group = 'study',
                                                                                                   cat1_name = 'Polyphen2 Prediction',
                                                                                                   cat2_name = 'Exp. Prediction',
                                                                                                   var_name = 'Count'),
                                                                        width = 12, height = 9)

deep_variant_plots$all_studies$pph_prediction_accuracy <- labeled_ggplot(p = plot_contingency_table(pph_summary,
                                                                                                     cat1 = 'pph2_class',
                                                                                                     cat2 = 'exp_prediction',
                                                                                                     var = 'prop', group = 'study',
                                                                                                     cat1_name = 'Polyphen2 Prediction',
                                                                                                     cat2_name = 'Exp. Prediction',
                                                                                                     var_name = 'Proportion'),
                                                                          width = 12, height = 9)

#### EVCouplings ####
select_evcoup <- function(x){
  if ('evcoup_epistatic' %in% names(x$multi_variants)){
    return(select(x$multi_variants, variants, score, evcoup_epistatic, evcoup_independent))
  } else if ('evcoup_epistatic' %in% names(x$single_variants)){
    return(select(x$single_variants, variants, score, evcoup_epistatic, evcoup_independent))
  } else {
    return(NULL)
  }
}

evcoup <- bind_rows(lapply(deep_variant_data, select_evcoup), .id='study') %>%
  mutate(exp_prediction = predict_exp_function(score),
         evcoup_prediction = ifelse(evcoup_epistatic < -6, MUT_CATEGORIES$deleterious, MUT_CATEGORIES$neutral))

deep_variant_plots$all_studies$evcoup_experimental_boxplot <- labeled_ggplot(plot_exp_pred_boxes(evcoup,
                                                                                                 y = 'evcoup_epistatic',
                                                                                                 y_name = 'EVCouplings Epistatic Score'),
                                                                          width = 12,
                                                                          height = 8)

# Correlation
evcoup_cor_test <- group_by(evcoup, study) %>%
  drop_na(score, evcoup_epistatic)

evcoup_cor_test <- full_join(do(evcoup_cor_test, tidy(cor.test(abs(.$score), .$evcoup_epistatic))),
                            summarise(evcoup_cor_test, multi = any(grepl(',', variants))),
                            by = 'study')

max_abs_t <- max(abs(evcoup_cor_test$statistic))
deep_variant_plots$all_studies$evcoup_correlation <- ggplot(evcoup_cor_test, aes(y=estimate, x=study, fill=multi)) + 
  geom_col() +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high), width=0.5) +
  geom_text(aes(label = ifelse(multi, '*', '')), position = position_stack(vjust = 1.1)) + 
  ggtitle('Correlation between EVCouplings Epistatic Score and abs(experimental scores)') +
  xlab('') +
  ylab('Pearson Correlation Coefficient') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# Predictions
evcoup_summary <- group_by(evcoup, study, evcoup_prediction)  %>%
  summarise(num_vars = n(),
            deleterious = sum(exp_prediction == MUT_CATEGORIES$deleterious),
            neutral = sum(exp_prediction == MUT_CATEGORIES$neutral)) %>%
  gather(key = 'exp_prediction', value = 'count', deleterious, neutral) %>%
  group_by(study, evcoup_prediction) %>%
  mutate(prop = count/sum(count)) %>%
  arrange(study, evcoup_prediction)

deep_variant_plots$all_studies$evcoup_prediction_counts <- labeled_ggplot(p = plot_contingency_table(evcoup_summary,
                                                                                                  cat1 = 'evcoup_prediction',
                                                                                                  cat2 = 'exp_prediction',
                                                                                                  var = 'count', group = 'study',
                                                                                                  cat1_name = 'EVCouplings Prediction',
                                                                                                  cat2_name = 'Exp. Prediction',
                                                                                                  var_name = 'Count'),
                                                                       width = 12, height = 9)

deep_variant_plots$all_studies$evcoup_prediction_accuracy <- labeled_ggplot(p = plot_contingency_table(evcoup_summary,
                                                                                                    cat1 = 'evcoup_prediction',
                                                                                                    cat2 = 'exp_prediction',
                                                                                                    var = 'prop', group = 'study',
                                                                                                    cat1_name = 'EVCouplings Prediction',
                                                                                                    cat2_name = 'Exp. Prediction',
                                                                                                    var_name = 'Proportion'),
                                                                         width = 12, height = 9)


#### Save plots #####
for (study_name in names(deep_variant_plots)){
  print(str_c('Writing figures for ', study_name))
  fig_root <- str_c('figures/variant_analysis/', study_name)
  dir.create(fig_root, showWarnings = FALSE, recursive = TRUE)
  for (plot_name in names(deep_variant_plots[[study_name]])){
    plot_path <- str_c(fig_root, '/', plot_name, '.pdf')
    # Only write figures that don't exist, delete figs to regenerate
    if (!file.exists(plot_path)){
      smart_save(deep_variant_plots[[study_name]][[plot_name]], plot_path, override=FALSE, width = 7, height = 5)
    }
  }
}
