#!/usr/bin/env Rscript 
# Script to analyse tool predictions on deep mutagenesis data

source('src/config.R')
source('src/analysis/data_processing.R')
source('src/analysis/tool_prediction_analysis.R')

# Import data
deep_variant_data <- readRDS('data/rdata/processed_variant_data.RDS')
all_variants <- readRDS('data/rdata/all_study_variants.RDS')
meta_df <- readRDS('data/rdata/study_meta_data.RDS')
variant_matrices <- readRDS('data/rdata/all_study_position_matrices.RDS')

# Generate basic plots per study
per_study_plots <- sapply(deep_variant_data, plot_predictions, simplify = FALSE)
save_plot_list(per_study_plots, root='figures/2_tool_prediction_analysis/per_study/')

# General plots across all studies
combined_plots <- list()

#### Envision ####
combined_plots$envision <- list()

select_envision <- function(x){
  if ('envision_prediction' %in% names(x$single_variants)){
    return(select(x$single_variants, variants, score, norm_score, envision_prediction, log2_envision_prediction))
  } else {
    return(NULL)
  }
}

envision <- bind_rows(lapply(deep_variant_data, select_envision), .id='study') %>%
  mutate(exp_prediction = exp_mut_class(score, study))

combined_plots$envision$pred_vs_score <- ggplot(envision, aes(x=norm_score,
                                                                           y=log2_envision_prediction)) +
  geom_point(shape=20) +
  xlab('Normalised Experimental Score') +
  ylab('Log2 Envision Prediction') +
  facet_wrap(~study, scales = 'free')

combined_plots$envision$experimental_boxplot <- labeled_ggplot(plot_exp_pred_boxes(envision,
                                                                                                y = 'envision_prediction',
                                                                                                y_name = 'Envision Prediction'),
                                                                            width = 12,
                                                                            height = 8)

envision_cor_test <- group_by(envision, study) %>%
  do(tidy(cor.test(.$score, .$log2_envision_prediction)))
max_abs_t <- max(abs(envision_cor_test$statistic))

combined_plots$envision$correlation <- ggplot(envision_cor_test, aes(y=estimate, x=study, fill=statistic)) + 
  geom_col() +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high), width=0.5) +
  ggtitle('Correlation between log2(Envision predictions) and mutagenesis scores') +
  xlab('') +
  ylab('Pearson Correlation Coefficient') +
  scale_fill_gradientn(colours = c('red', 'white', 'blue'), limits = c(-max_abs_t, max_abs_t)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
########

#### SIFT ####
combined_plots$sift <- list()

select_sift <- function(x){
  if ('sift_prediction' %in% names(x$single_variants)){
    return(select(x$single_variants, variants, score, norm_score, sift_prediction, sift_score, sift_median))
  } else {
    return(NULL)
  }
}

sift <- bind_rows(lapply(deep_variant_data, select_sift), .id='study') %>%
  mutate(exp_prediction = exp_mut_class(score, study),
         sift_prediction = str_to_lower(gsub('TOLERATED', MUT_CATEGORIES$neutral, sift_prediction)))

combined_plots$sift$pred_vs_score <- ggplot(sift, aes(x=norm_score,
                                                                   y=sift_score)) +
  geom_point(shape=20) +
  xlab('Normalised Experimental Score') +
  ylab('SIFT Score') +
  scale_y_log10() +
  facet_wrap(~study, scales = 'free')

combined_plots$sift$experimental_boxplot <- labeled_ggplot(plot_exp_pred_boxes(sift,
                                                                                            y = 'sift_score',
                                                                                            y_name = 'SIFT Score'),
                                                                        width = 12,
                                                                        height = 8)

# Correlation
sift_cor_test <- group_by(sift, study) %>%
  drop_na(score, sift_score) %>%
  do(tidy(cor.test(.$score, log10(.$sift_score + min(.$sift_score[.$sift_score > 0], na.rm = TRUE)))))
max_abs_t <- max(abs(sift_cor_test$statistic), na.rm = TRUE)

combined_plots$sift$correlation <- ggplot(sift_cor_test, aes(y=estimate, x=study, fill=statistic)) + 
  geom_col() +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high), width=0.5) +
  ggtitle('Correlation between log10(SIFT Score) and mutagenesis scores') +
  xlab('') +
  ylab('Pearson Correlation Coefficient') +
  scale_fill_gradientn(colours = c('red', 'white', 'blue'), limits = c(-max_abs_t, max_abs_t)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

sift_abs_cor_test <- group_by(sift, study) %>%
  drop_na(score, sift_score) %>%
  do(tidy(cor.test(abs(.$score), -log10(.$sift_score + min(.$sift_score[.$sift_score > 0], na.rm = TRUE)))))
max_abs_t <- max(abs(sift_cor_test$statistic), na.rm = TRUE)

combined_plots$sift$correlation_abs <- ggplot(sift_abs_cor_test, aes(y=estimate, x=study, fill=statistic)) + 
  geom_col() +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high), width=0.5) +
  ggtitle('Correlation between -log10(SIFT Score) and abs(mutagenesis scores)') +
  xlab('') +
  ylab('Pearson Correlation Coefficient') +
  scale_fill_gradientn(colours = c('white', 'blue'), limits = c(0, max_abs_t)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# Predictions
sift_summary <- group_by(sift, study, sift_prediction) %>%
  drop_na(sift_prediction, exp_prediction) %>%
  summarise(num_vars = n(),
            deleterious = sum(exp_prediction == MUT_CATEGORIES$deleterious),
            neutral = sum(exp_prediction == MUT_CATEGORIES$neutral)) %>%
  gather(key = 'exp_prediction', value = 'count', deleterious, neutral) %>%
  group_by(study, sift_prediction) %>%
  mutate(prop = count/sum(count)) %>%
  arrange(study, sift_prediction)

combined_plots$sift$prediction_counts <- labeled_ggplot(
  p = plot_contingency_table(sift_summary,
                             cat1 = 'sift_prediction',
                             cat2 = 'exp_prediction',
                             var = 'count', group = 'study',
                             cat1_name = 'SIFT Prediction',
                             cat2_name = 'Exp. Prediction',
                             var_name = 'Count'),
  width = 12, height = 9)

combined_plots$sift$prediction_accuracy <- labeled_ggplot(
  p = plot_contingency_table(sift_summary,
                             cat1 = 'sift_prediction',
                             cat2 = 'exp_prediction',
                             var = 'prop', group = 'study',
                             cat1_name = 'SIFT Prediction',
                             cat2_name = 'Exp. Prediction',
                             var_name = 'Proportion'),
  width = 12, height = 9)
########

#### FoldX ####
combined_plots$foldx <- list()

select_foldx <- function(x){
  if (any(grepl('foldx_', names(x$multi_variants)))){
    tbl <- x$multi_variants
  } else if (any(grepl('foldx_', names(x$single_variants)))){
    tbl <- x$single_variants
  } else {
    return(NULL)
  }
  tbl <- select(tbl, variants, score, norm_score, starts_with('foldx_')) %>%
    rename_at(vars(contains('foldx_')), .funs = ~ gsub('foldx_', '', .)) %>%
    gather(key = 'k', value = 'v', -variants, -score, -norm_score) %>%
    separate(k, into = c('pdb_id', 'k'), sep = '_') %>%
    spread(key = k, value = v)
  return(tbl)
}

foldx <- bind_rows(lapply(deep_variant_data, select_foldx), .id='study') %>%
  mutate(count = factor(sapply(variants, function(x){dim(str_split(x, ',', simplify = TRUE))[2]})),
         single = count == 1,
         exp_prediction = exp_mut_class(score, study))

combined_plots$foldx$ddG_vs_score <- ggplot(filter(foldx, count==1), aes(x=norm_score, y=ddG)) +
  geom_point(shape=20) +
  xlab('Normalised Experimental Score') +
  ylab('FoldX ddG') +
  facet_wrap(~study, scales = 'free')

combined_plots$foldx$experimental_boxplot <- labeled_ggplot(plot_exp_pred_boxes(foldx,
                                                                                             y = 'ddG'),
                                                                         width = 12,
                                                                         height = 8)

combined_plots$foldx$experimental_boxplot_single <- labeled_ggplot(plot_exp_pred_boxes(filter(foldx, single),
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
combined_plots$foldx$correlation <- ggplot(foldx_cor_test, aes(y=estimate, x=name, fill=statistic)) + 
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

combined_plots$foldx$foldx_prediction_counts <- labeled_ggplot(
  p = plot_contingency_table(foldx_summary,
                             cat1 = 'foldx_prediction',
                             cat2 = 'exp_prediction',
                             var = 'count', group = 'study',
                             cat1_name = 'abs(ddG) > 2',
                             cat2_name = 'Exp. Prediction',
                             var_name = 'Count'),
  width = 12, height = 9)

combined_plots$foldx$foldx_prediction_accuracy <- labeled_ggplot(
  p = plot_contingency_table(foldx_summary,
                             cat1 = 'foldx_prediction',
                             cat2 = 'exp_prediction',
                             var = 'prop', group = 'study',
                             cat1_name = 'abs(ddG) > 2',
                             cat2_name = 'Exp. Prediction',
                             var_name = 'Proportion'),
  width = 12, height = 9)
########

#### Polyphen2 ####
combined_plots$polyphen2 <- list()

select_pph <- function(x){
  if ('pph2_class' %in% names(x$single_variants)){
    return(select(x$single_variants, variants, score, norm_score, pph2_prediction, pph2_class, pph2_prob, pph2_FPR, pph2_TPR, pph2_FDR))
  } else {
    return(NULL)
  }
}

pph <- bind_rows(lapply(deep_variant_data, select_pph), .id='study') %>%
  mutate(exp_prediction = exp_mut_class(score, study)) %>%
  drop_na(pph2_prob)

combined_plots$polyphen2$prob_vs_score <- ggplot(pph, aes(x=norm_score, y=pph2_prob + min(pph2_prob[pph2_prob > 0]))) +
  geom_point(shape=20) +
  xlab('Normalised Experimental Score') +
  ylab('log10(Polyphen2 Deleterious Probability)') +
  facet_wrap(~study) +
  scale_y_log10()

combined_plots$polyphen2$experimental_boxplot <- labeled_ggplot(
  plot_exp_pred_boxes(pph,
                      y = 'pph2_prob',
                      y_name = 'Polyphen2 Probability'),
  width = 12,
  height = 8)

# Correlation
pph_cor_test <- group_by(pph, study) %>%
  drop_na(score, pph2_prob) %>%
  do(tidy(cor.test(.$score, -log10(.$pph2_prob + min(.$pph2_prob[.$pph2_prob > 0], na.rm = TRUE)))))
max_abs_t <- max(abs(pph_cor_test$statistic), na.rm = TRUE)

combined_plots$polyphen2$correlation <- ggplot(pph_cor_test, aes(y=estimate, x=study, fill=statistic)) + 
  geom_col() +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high), width=0.5) +
  ggtitle('Correlation between -log10(Polyphen2 Probability) and mutagenesis scores') +
  xlab('') +
  ylab('Pearson Correlation Coefficient') +
  scale_fill_gradientn(colours = c('red', 'white', 'blue'), limits = c(-max_abs_t, max_abs_t)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# Predictions
pph_summary <- group_by(pph, study, pph2_class) %>%
  drop_na(pph2_class, exp_prediction) %>%
  summarise(num_vars = n(),
            deleterious = sum(exp_prediction == MUT_CATEGORIES$deleterious),
            neutral = sum(exp_prediction == MUT_CATEGORIES$neutral)) %>%
  gather(key = 'exp_prediction', value = 'count', deleterious, neutral) %>%
  group_by(study, pph2_class) %>%
  mutate(prop = count/sum(count)) %>%
  arrange(study, pph2_class)

combined_plots$polyphen2$prediction_counts <- labeled_ggplot(
  p = plot_contingency_table(pph_summary,
                             cat1 = 'pph2_class',
                             cat2 = 'exp_prediction',
                             var = 'count', group = 'study',
                             cat1_name = 'Polyphen2 Prediction',
                             cat2_name = 'Exp. Prediction',
                             var_name = 'Count'),
  width = 12, height = 9)

combined_plots$polyphen2$prediction_accuracy <- labeled_ggplot(
  p = plot_contingency_table(pph_summary,
                             cat1 = 'pph2_class',
                             cat2 = 'exp_prediction',
                             var = 'prop', group = 'study',
                             cat1_name = 'Polyphen2 Prediction',
                             cat2_name = 'Exp. Prediction',
                             var_name = 'Proportion'),
  width = 12, height = 9)
########

#### EVCouplings ####
combined_plots$evcouplings <- list()

select_evcoup <- function(x){
  if ('evcoup_epistatic' %in% names(x$multi_variants)){
    return(select(x$multi_variants, variants, score, norm_score, evcoup_epistatic, evcoup_independent))
  } else if ('evcoup_epistatic' %in% names(x$single_variants)){
    return(select(x$single_variants, variants, score, norm_score, evcoup_epistatic, evcoup_independent))
  } else {
    return(NULL)
  }
}

evcoup <- bind_rows(lapply(deep_variant_data, select_evcoup), .id='study') %>%
  mutate(exp_prediction = exp_mut_class(score, study),
         evcoup_prediction = ifelse(evcoup_epistatic < -6, MUT_CATEGORIES$deleterious, MUT_CATEGORIES$neutral))

combined_plots$evcouplings$evcoup_vs_score <- ggplot(
  evcoup, aes(x=norm_score,
              y=evcoup_epistatic)) +
  geom_point(shape=20) +
  xlab('Normalised Experimental Score') +
  ylab('EVCouplings Epistatic Score') +
  facet_wrap(~study, scales = 'free')

combined_plots$evcouplings$experimental_boxplot <- labeled_ggplot(
  plot_exp_pred_boxes(evcoup,
                      y = 'evcoup_epistatic',
                      y_name = 'EVCouplings Epistatic Score'),
  width = 12,
  height = 8)

# Correlation
evcoup_cor_test <- group_by(evcoup, study) %>%
  drop_na(score, evcoup_epistatic)

evcoup_cor_test <- full_join(do(evcoup_cor_test, tidy(cor.test(abs(.$score), -1*.$evcoup_epistatic))),
                             summarise(evcoup_cor_test, multi = any(grepl(',', variants))),
                             by = 'study')

max_abs_t <- max(abs(evcoup_cor_test$statistic))
combined_plots$evcouplings$correlation <- ggplot(evcoup_cor_test, aes(y=estimate, x=study, fill=multi)) + 
  geom_col() +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high), width=0.5) +
  geom_text(aes(label = ifelse(multi, '*', '')), position = position_stack(vjust = 1.1)) + 
  ggtitle('Correlation between -(EVCouplings Epistatic Score) and abs(experimental scores)') +
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

combined_plots$evcouplings$prediction_counts <- labeled_ggplot(
  p = plot_contingency_table(evcoup_summary,
                             cat1 = 'evcoup_prediction',
                             cat2 = 'exp_prediction',
                             var = 'count', group = 'study',
                             cat1_name = 'EVCouplings Prediction',
                             cat2_name = 'Exp. Prediction',
                             var_name = 'Count'),
  width = 12, height = 9)

combined_plots$evcouplings$prediction_accuracy <- labeled_ggplot(
  p = plot_contingency_table(evcoup_summary,
                             cat1 = 'evcoup_prediction',
                             cat2 = 'exp_prediction',
                             var = 'prop', group = 'study',
                             cat1_name = 'EVCouplings Prediction',
                             cat2_name = 'Exp. Prediction',
                             var_name = 'Proportion'),
  width = 12, height = 9)
########

#### Mixture Plots ####
# Sift/FoldX Score Correlation with exp score

fx <- ungroup(foldx_cor_test) %>%
  filter(single) %>%
  mutate(study=sapply(name, function(x){str_sub(x, end=-7)}),
         pdb_id=sapply(name, function(x){str_sub(x, start=-5, end=-2)})) %>%
  select(study, pdb_id, estimate, statistic, p.value, parameter, conf.low, conf.high)

si <- select(sift_abs_cor_test, study, estimate, statistic, p.value, parameter, conf.low, conf.high)

sift_foldx_cor <- bind_rows(SIFT=filter(si, study %in% fx$study),
                            FoldX=fx, .id = 'tool') %>%
  mutate(study_pretty = sapply(study, format_study),
         p_cat = cut(p.value, breaks = c(0, 1e-12, 1e-06, 1e-3, 0.01, 0.05, 1)))

laber <- function(df){
  return(lapply(df[,1], function(x){)[x]}))
}

combined_plots$sift_foldx_correlations <- labeled_ggplot(
  p=ggplot(sift_foldx_cor, aes(y=estimate, x=study_pretty, group=pdb_id, fill=p_cat)) + 
    facet_wrap(~tool, ncol = 1,
               labeller = labeller(tool=as_labeller(x=c(SIFT="-log[10]~SIFT", FoldX="FoldX~Delta*Delta*G"),
                                                    default = label_parsed))) +
    geom_col(position = position_dodge()) +
    # To ylim(0,...) elegantly use sapply(conf.low, function(x){max(x, 0)})
    geom_errorbar(aes(ymin=conf.low, ymax=conf.high), width=0.5, position = position_dodge(0.9)) +
    geom_hline(yintercept = 0) +
    ggtitle("Correlation between |ER| and predictions") +
    xlab('') +
    ylab('Pearson Correlation') +
    scale_fill_viridis_d(guide=guide_legend(title='p-value'), direction = -1, drop=FALSE) +
    theme_pubclean() +
    theme(plot.title = element_text(size=30, hjust = 0.3),
          axis.text = element_text(size=11),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title.y = element_text(size=20),
          legend.title = element_text(size=15),
          legend.text = element_text(size=10),
          legend.position = 'right',
          strip.background = element_blank(),
          strip.text = element_text(size = 20)),
  width = 20, height = 20, units = 'cm')

combined_plots$sift_foldx_correlations <- labeled_ggplot(combined_plots$sift_foldx_correlations,
                                                         width=9, height=7)

# All tool correlations
ev <- select(evcoup_cor_test, study, estimate, statistic, p.value, parameter, conf.low, conf.high, multi)
en <- select(envision_cor_test, study, estimate, statistic, p.value, parameter, conf.low, conf.high)
pp <- select(pph_cor_test, study, estimate, statistic, p.value, parameter, conf.low, conf.high)

tool_titles <- c(si='SIFT: -log10(SIFT Score) vs abs(Exp Score)',
                 fx='FoldX: ddG vs abs(Exp Score)',
                 ev='EVCouplings: Epistatic Score vs abs(Exp Score)',
                 en='Envision: log2(Envision Score) vs Exp Score',
                 pp='Polyphen2: -log10(PPH Score) vs Exp Score')

all_cor <- bind_rows(si=si, fx=fx, ev=ev, en=en, pp=pp, .id = 'tool') %>%
  mutate(study_pretty = sapply(study, format_study),
         p_cat = cut(p.value, breaks = c(0, 1e-12, 1e-06, 1e-3, 0.01, 0.05, 1), include.lowest = TRUE),
         tool_title=tool_titles[tool])

combined_plots$all_correlations <- ggplot(all_cor, aes(y=estimate, x=study_pretty,
                                                                    group=pdb_id, fill=p_cat)) + 
  facet_wrap(~tool_title, ncol = 1) +
  geom_col(position = position_dodge()) +
  # To ylim(0,...) elegantly use sapply(conf.low, function(x){max(x, 0)})
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high), width=0.5, position = position_dodge(0.9)) +
  geom_hline(yintercept = 0) +
  ggtitle('Correlation between |ER| and tool scores') +
  #  ylim(0, 0.75) +
  xlab('') +
  ylab('Pearson Correlation Coefficient') +
  scale_fill_viridis_d(guide=guide_legend(title='p-value'), direction = -1, drop=FALSE) +
  theme_pubclean() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = 'right',
        strip.background = element_blank())
  
combined_plots$all_correlations <- labeled_ggplot(combined_plots$all_correlations, width=18, height=25, units = 'cm')
########

#### SIFT/FoldX Position Profiles vs DMS Position Profiles ####
prediction_martices <- list()
prediction_martices$sift <- bind_rows(lapply(deep_variant_data, make_var_matrix, score='sift_score'), .id = 'study') %>%
  mutate(sift = TRUE)
prediction_martices$foldx <- bind_rows(lapply(deep_variant_data, make_foldx_var_matrix), .id = 'study')  %>%
  mutate(foldx = TRUE)
prediction_martices$exp <- select(variant_matrices$norm_all_variants, study, pos, wt, A:Y)

pairwise_dists <- mapply(compute_per_gene_pairwise_profile_dists,
                         prediction_martices,
                         str_c(names(prediction_martices),'_dist'),
                         SIMPLIFY = FALSE) %>%
  reduce(left_join, by = c('study', 'pos1', 'aa1', 'pos2', 'aa2'))

combined_plots$foldx$per_position_profile_cor <- ggplot(drop_na(pairwise_dists, foldx_dist), aes(x=exp_dist, y=foldx_dist, colour=study)) + 
  geom_density_2d() +
  geom_smooth(method = 'lm', colour='black') +
  facet_wrap(~study) + 
  theme(legend.position = 'none')

combined_plots$sift$per_position_profile_cor <- ggplot(drop_na(pairwise_dists, sift_dist), aes(x=exp_dist, y=sift_dist, colour=study)) + 
  geom_density_2d() +
  geom_smooth(method = 'lm', colour='black') +
  facet_wrap(~study) + 
  theme(legend.position = 'none')
########

#### ROC/PR Curves ####
study_cols <- c('variants', 'score', 'norm_score', 'evcoup_epistatic', 'evcoup_independent', 'pph2_prob', 'envision_prediction', 'sift_score')
variants <- bind_rows(lapply(deep_variant_data, function(x){select(x$single_variants, one_of(study_cols), starts_with('foldx_'))}), .id='study') %>%
  mutate(exp_pred = exp_mut_class(score, study),
         neutral = exp_pred == 'neutral',
         deleterious = exp_pred == 'deleterious') %>%
  mutate(xfold_ddG = tibble_to_matrix(., starts_with('foldx_')) %>% rowMeans(., na.rm = TRUE)) %>%
  select(-starts_with('foldx_'), foldx_ddG=xfold_ddG)
variants$foldx_ddG[is.nan(variants$foldx_ddG)] <- NA
variants <- group_by(variants, study)

study_roc_vals <- bind_rows(
  .id = 'tool',
  EVCouplings=do(variants, calc_true_false_over_range(., deleterious, evcoup_epistatic, comparison_func = function(x, y){x < y})),
  SIFT=do(variants, calc_true_false_over_range(., deleterious, sift_score, comparison_func = function(x, y){x < y})),
  FoldX=do(variants, calc_true_false_over_range(., deleterious, foldx_ddG, comparison_func = function(x, y){x > y})),
  PolyPhen2=do(variants, calc_true_false_over_range(., deleterious, pph2_prob, comparison_func = function(x, y){x > y})),
  Envision=do(variants, calc_true_false_over_range(., deleterious, envision_prediction, comparison_func = function(x, y){x < y}))
) %>%
  mutate(TPR = TP/(TP + FN),
         FPR = FP/(FP + TN),
         precision = TP/(TP + FP))
study_roc_vals$precision[is.nan(study_roc_vals$precision)] <- 0

combined_plots$roc_curve <- ggplot(study_roc_vals, aes(x=FPR, y=TPR, colour=study)) +
  geom_line() +
  geom_abline(slope = 1) +
  facet_wrap(~tool, nrow = 1) +
  guides(colour=FALSE) +
  theme_pubclean(base_size = 10)

combined_plots$pr_curve <- ggplot(study_roc_vals, aes(x=TPR, y=precision, colour=study)) +
  geom_line() +
  facet_wrap(~tool, nrow = 1) +
  guides(colour=FALSE) +
  xlab('Recall') +
  ylab('Precision') +
  theme_pubclean(base_size = 10)

combined_plots$roc_curve_si_fx <- labeled_ggplot(
  p=ggplot(filter(study_roc_vals, tool %in% c('SIFT', 'FoldX')), aes(x=FPR, y=TPR, colour=study)) +
    geom_line() +
    geom_segment(x = 0, y = 0, xend = 1, yend = 1, linetype='dotted', colour='black') +
    facet_wrap(~tool) +
    guides(colour=FALSE) +
    theme_pubclean() +
    theme(strip.background = element_blank(),
          strip.text = element_text(size = 25),
          axis.title = element_text(size = 25),
          axis.text = element_text(size=13)),
  height=13.75, width=19.5, units = 'cm')

combined_plots$pr_curve_si_fx <- labeled_ggplot(
  p=ggplot(filter(study_roc_vals, tool %in% c('SIFT', 'FoldX')), aes(x=TPR, y=precision, colour=study)) +
    geom_line() +
    facet_wrap(~tool) +
    guides(colour=FALSE) +
    xlab('Recall') +
    ylab('Precision') +
    theme_pubclean() +
    theme(strip.background = element_blank(),
          strip.text = element_text(size = 25),
          axis.title = element_text(size = 25),
          axis.text = element_text(size=13)),
  height=13.75, width=19.5, units = 'cm')

combined_plots$pr_roc_curves <- labeled_ggplot(
  p=gridExtra::arrangeGrob(combined_plots$roc_curve + 
                             labs(tag = 'A') + 
                             theme(strip.background = element_blank(),
                                   axis.text.x = element_text(angle = 90, vjust = 0.5)),
                           combined_plots$pr_curve + 
                             labs(tag = 'B') + 
                             theme(strip.background = element_blank(),
                                   axis.text.x = element_text(angle = 90, vjust = 0.5)),
                           ncol = 1),
  units = 'cm', width = 18, height = 10
)
########

# Save all study plots 
save_plot_list(combined_plots, root='figures/2_tool_prediction_analysis/')
