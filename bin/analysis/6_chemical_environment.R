#!/usr/bin/env Rscript 
# Script analysing chemical environment in proteins, in the context of our deep mutagenesis data

source('src/config.R')
source('src/analysis/chemical_environment.R')
source('src/analysis/position_profile_clustering.R')

foldx <- readRDS('data/rdata/all_foldx.RDS')

variant_matrices <- readRDS('data/rdata/all_study_position_matrices.RDS')

chemical_environments <- readRDS('data/rdata/position_chemical_environments.RDS') %>%
  left_join(., select(variant_matrices$norm_all_variants, study, position=pos, aa=wt, sig_count, A:Y), by = c("study", "position", "aa")) %>%
  mutate_at(.vars = vars(one_of(c('nearest_10', 'within_10.0'))),
            .funs = list(reduced=~sapply(., reduce_aa_profile, simplify = FALSE)))

# chem_env_profiles <- c('nearest_10', 'nearest_10_reduced', 'within_10.0', 'within_10.0_reduced')
# chem_env_profile_col_names <- list(str_c('prof_', SORTED_AA_1_CODE),
#                                    str_c('prof_', names(AA_REDUCED_CLASSES)),
#                                    str_c('prof_', SORTED_AA_1_CODE),
#                                    str_c('prof_', names(AA_REDUCED_CLASSES)))
chem_env_profiles <- 'within_10.0'
chem_env_profile_col_names <- list(str_c('prof_', SORTED_AA_1_CODE))

chem_env_analyses <- list()
plots <- list()

for (i in 1:length(chem_env_profiles)){
  prof_plots <- list()
  prof_col <- chem_env_profiles[i]
  prof_col_names <- chem_env_profile_col_names[[i]]
  prof_col_syms <- syms(prof_col_names)
  
  chem_env <- select(chemical_environments, study, pdb_id, chain, position, relative_position, aa, aa_reduced, sig_count, 
                     !!prof_col, gene_type, gene_name, test_class, species, authour, group, single, thresh, norm_thresh, factor,
                     ss, ss_reduced, all_atom_abs:polar_rel, A:Y) %>%
    expand_profile_column(!!prof_col, names=prof_col_names) %>%
    mutate(duplicate_position = duplicated(select(chemical_environments, pdb_id, chain, position, aa)))
  
  ## Correlation between profile components
  prof_component_cors <- tibble_correlation(chem_env, !!! prof_col_syms, filter_diag = TRUE) %>%
    mutate_at(vars(cat1, cat2), .funs = ~ str_sub(., start=-1))
  prof_plots$profile_correlation_heatmap <- ggplot(prof_component_cors, aes(x=cat1, y=cat2, fill=cor)) +
    geom_tile() +
    xlab('') +
    ylab('') +
    scale_fill_gradient2() +
    theme(axis.ticks = element_blank(), panel.background = element_blank())
  
  ## Paired plots of of profile components 
  # commented out as long calculation but not particularly interesting, but kept for info
  # num_cats <- length(prof_col_names)
  # prof_plots$profile_pairs = labeled_ggplot(
  #   ggpairs(tbl, columns = prof_col_names, lower = list(continuous=function(d, m, ...){ggplot(d, m, ...) + geom_count()})),
  #   height = num_cats * 2, width = num_cats * 2, limitsize=FALSE)
  
  ## PCA analysis of profile
  pca <- tibble_pca(filter(chem_env, !duplicate_position), !!! prof_col_syms)
  chem_env <- bind_cols(chem_env, as_tibble(scale(select(chem_env, !!! prof_col_syms), pca$center, pca$scale) %*% pca$rotation))
  
  num_pcs <- length(grep('^PC[0-9]*$', names(chem_env)))
  max_plot_pc <- if(num_pcs %% 2 == 0) num_pcs else num_pcs - 1 # Work around to current rigid PC plotting
  plot_row_cols <- get_good_rows_cols(max_plot_pc/2)
  pca_plots <- plot_chem_env_basic_pca_plots(chem_env,
                                             max_pc = max_plot_pc,
                                             nrow = plot_row_cols[1], ncol = plot_row_cols[2],
                                             cont_factors=c('all_atom_rel', 'relative_position', 'sig_count'),
                                             discrete_factors=c('ss', 'ss_reduced', 'aa', 'aa_reduced',
                                                                'pdb_id', 'gene_name', 'species'))
  
  pca_factor_cors <- pca_factor_cor(list(pca=pca, profiles=chem_env),
                                    .vars = vars(all_atom_abs:polar_rel, relative_position, sig_count))
  
  prof_plots$pca <- c(pca_plots, list(factor_heatmap=pca_factor_heatmap(pca_factor_cors)))
  
  ## tSNE analysis of profile
  tsne <- tibble_tsne(chem_env, !!! prof_col_syms)
  
  chem_env <- bind_cols(chem_env, as_tibble(set_colnames(tsne$tsne$Y[tsne$unique_row_indeces,], c('tSNE1', 'tSNE2'))))      
  
  prof_plots$tsne <- plot_factors(filter(chem_env, !duplicate_position),
                                  tSNE1, tSNE2, quos(ss_reduced=ss_reduced, aa_reduced=aa_reduced, gene_name=gene_name,
                                                sig_count=sig_count, sqrt_suf_acc=sqrt(all_atom_rel)))
  
  ## LMs predicting ER for a given AA from chem env profile
  prof_plots$lm <- list()
  
  # Only use profile
  prof_lm <- calc_all_profile_lms(chem_env, prof_vars = vars(!!! prof_col_syms), target_vars = vars(A:Y),
                                  include_intercept = TRUE, per_study = FALSE)
  
  prof_plots$lm[['er_profile_only']] <- basic_prof_lm_plots(prof_lm)
  
  # Add significance of position
  prof_lm_sig <- calc_all_profile_lms(chem_env, prof_vars = vars(!!! prof_col_syms, sig_count), target_vars = vars(A:Y),
                                      include_intercept = FALSE, per_study = FALSE)
  
  prof_plots$lm[['er_profile_sig_count']] <- basic_prof_lm_plots(prof_lm_sig)
  
  # Predict sig count from profile
  sig_count_lm <- calc_profile_lm(chem_env, target_col=sig_count, !!!prof_col_syms, include_intercept = FALSE)
  summ <- summary(sig_count_lm)
  
  prof_plots$lm[['sig_count']] <- list()
  prof_plots$lm$sig_count$preds <- ggplot(augment(sig_count_lm), aes(x=sig_count, y=.fitted, colour=.resid)) +
    geom_point() + 
    geom_smooth(method = 'lm') +
    guides(colour=guide_colourbar(title='Residual')) +
    ylab('Predicted # Significant Substitutions') +
    xlab('# Significant Substitutions') +
    geom_text(label = str_c('r^2 = ', signif(summ$r.squared, 4)), x=15, y=1, hjust = 0)
  
  loadings <- tidy(sig_count_lm) %>% mutate(x=1)
  prof_plots$lm$sig_count$loadings <- ggplot(loadings, aes(x=x, y=term, fill=estimate)) + 
    geom_tile() +
    scale_fill_gradient2() +
    theme(axis.ticks = element_blank(), panel.background = element_blank()) +
    xlab('') +
    ylab('LM Term') + 
    scale_shape_manual(values = c('p < 0.0001'=8, 'p < 0.001'=3, 'p < 0.01'=20))
    
  if (any(loadings$p.value < 0.0001)){
    prof_plots$lm$sig_count$loadings <- prof_plots$lm$sig_count$loadings +
      geom_point(data = filter(loadings, p.value < 0.0001), aes(shape='p < 0.0001'))
  }
  if (any(loadings$p.value < 0.001 & loadings$p.value > 0.0001)){
    prof_plots$lm$sig_count$loadings <- prof_plots$lm$sig_count$loadings +
      geom_point(data = filter(loadings, p.value < 0.001, p.value > 0.0001), aes(shape='p < 0.001'))
  }
  if (any(loadings$p.value < 0.01 & loadings$p.value > 0.001)){
    prof_plots$lm$sig_count$loadings <- prof_plots$lm$sig_count$loadings +
      geom_point(data = filter(loadings, p.value < 0.01, p.value > 0.001), aes(shape='p < 0.01'))
  }
    
  ## Save results and plots for specific analysis
  plots[[prof_col]] <- prof_plots
  chem_env_analyses[[prof_col]] <- list(tbl=chem_env,
                                        pca=pca,
                                        tsne=tsne,
                                        lm_er=prof_lm,
                                        lm_er_sig=prof_lm_sig,
                                        lm_sig=sig_count_lm)
}

#### Top/Bottom environments for substitutions ####
# example workflow using Lysine
# Can be incorporated in general analysis function and extended to all AAs when ready

## Calc top/bottom environments
top_5_percent <- group_by(chem_env_analyses$within_10.0$tbl, study, position, aa) %>%
  summarise_all(.funs = first) %>%
  ungroup() %>%
  top_frac(0.05, K) %>%
  mutate(frac='top 5%')
bot_5_percent <- group_by(chem_env_analyses$within_10.0$tbl, study, position, aa) %>%
  summarise_all(.funs = first) %>%
  ungroup() %>%
  top_frac(-0.05, K) %>%
  mutate(frac='bottom 5%')

top_bot <- bind_rows(top_5_percent, bot_5_percent) %>%
  select(study, position, aa, pdb_id, gene_name, ss, ss_reduced, all_atom_rel, aa_reduced, relative_position,
         sig_count, A:Y, prof_A:prof_Y, PC1:PC20, tSNE1, tSNE2, frac)

# Plot profile count distributions
plots$within_10.0$top_bot_K_sub_profiles <- select(top_bot, frac, prof_A:prof_Y) %>%
  gather('aa', 'count', -frac) %>%
  mutate(aa = str_sub(aa, start = -1)) %>% 
  ggplot(., aes(x=aa, y=count, colour=frac)) + 
  geom_count(position = position_dodge(0.75)) +
  scale_fill_manual(values = c(`bottom 5%`='red', `top 5%`='blue'))

## Calc correlation between AA counts in profiles
top_bot_cors <- group_by(top_bot, frac) %>% 
  do(tibble_correlation(., prof_A:prof_Y)) %>%
  mutate_at(vars(cat1, cat2), .funs = ~ as.character(str_sub(., start = -1))) %>%
  rename(AA1=cat1, AA2=cat2)

plots$within_10.0$top_bot_K_prof_cors <- ggplot(top_bot_cors, aes(x=AA1, y=AA2, fill=cor)) +
  facet_wrap(~frac, nrow = 1) +
  geom_tile() +
  xlab('') +
  ylab('') +
  scale_fill_gradient2() +
  theme(axis.ticks = element_blank(), panel.background = element_blank(),
        strip.background = element_blank()) +
  ggtitle('Correlation between AAs in profiles for top and bottom 5% of K substitutions')

# Relative correlation change top - bot
top_bot_rel_cors <- spread(top_bot_cors, frac, cor) %>%
  mutate(rel_cor = `top 5%` - `bottom 5%`) %>%
  add_factor_order(AA1, AA2, rel_cor)

plots$within_10.0$lm$top_bot_k_profs_rel_cors <- ggplot(top_bot_rel_cors, aes(x=AA1, y=AA2, fill=rel_cor)) +
  geom_tile() +
  xlab('') +
  ylab('') +
  scale_fill_gradient2() +
  theme(axis.ticks = element_blank(), panel.background = element_blank(),
        strip.background = element_blank()) +
  ggtitle('Correlation in top 5% K substitutions - correlation in bottom 5%')
########

#### Examine best and worst predicted environments from LMs ####
lm_preds_tbl <- pull(chem_env_analyses$within_10.0$lm_er,model) %>%
  sapply(function(x){predict(x, newdata=chem_env_analyses$within_10.0$tbl)}) %>%
  set_colnames(str_c('pred_', chem_env_analyses$within_10.0$lm_er$target)) %>%
  as_tibble() %>%
  bind_cols(chem_env_analyses$within_10.0$tbl, .)

pred_top_5_percent <- group_by(lm_preds_tbl, study, position, aa) %>%
  summarise_all(.funs = first) %>%
  ungroup() %>%
  drop_na(sig_count) %>%
  top_frac(0.05, pred_K) %>%
  mutate(frac='top 5%')
pred_bot_5_percent <- group_by(lm_preds_tbl, study, position, aa) %>%
  summarise_all(.funs = first) %>%
  ungroup() %>%
  drop_na(sig_count) %>%
  top_frac(-0.05, pred_K) %>%
  mutate(frac='bottom 5%')
pred_top_bot <- bind_rows(pred_top_5_percent, pred_bot_5_percent) %>%
  select(frac, study, position, aa, pdb_id, gene_name, ss, ss_reduced, all_atom_rel, aa_reduced, relative_position,
         sig_count, A:Y, prof_A:prof_Y, PC1:PC20, tSNE1, tSNE2, pred_A:pred_Y)

plots$within_10.0$lm$er_profile_only$top_bot_pred_K_sub_profiles <- select(pred_top_bot, frac, prof_A:prof_Y) %>%
  gather('aa', 'count', -frac) %>%
  mutate(aa = str_sub(aa, start = -1)) %>% 
  ggplot(., aes(x=aa, y=count, colour=frac)) + 
  geom_count(position = position_dodge(0.75)) +
  scale_fill_manual(values = c(`bottom 5%`='red', `top 5%`='blue'))
plots$within_10.0$lm$er_profile_only$top_bot_pred_K_sig_count <- ggplot(pred_top_bot, aes(x=frac, y=sig_count, colour=frac)) + 
  geom_boxplot() +
  scale_fill_manual(values = c(`bottom 5%`='red', `top 5%`='blue'))

# FoldX predictions at top and bottom predictions
chem_env_foldx <- select(lm_preds_tbl, study, pdb_id, position, wt=aa, gene_name, A:Y, pred_A:pred_Y, prof_A:prof_Y) %>%
  gather(key = 'mut', value = 'er', A:Y, pred_A:pred_Y) %>%
  mutate(pred = ifelse(grepl('pred\\_', mut), 'pred', 'exp'),
         mut = str_sub(mut, start=-1)) %>%
  spread(pred, er) %>%
  select(study, gene_name, pdb_id, position, wt, mut, exp, pred, prof_A:prof_Y) %>%
  left_join(., foldx, by=c('study', 'pdb_id', 'position', 'wt', 'mut')) %>%
  drop_na(total_energy)

chem_env_foldx_top_bot <- bind_rows(
  `top 5%` = sapply(unique(chem_env_foldx$mut), function(x){filter(chem_env_foldx, mut == x) %>% top_frac(0.05, pred)},
                    simplify = FALSE) %>% bind_rows(),
  `bottom 5%` = sapply(unique(chem_env_foldx$mut), function(x){filter(chem_env_foldx, mut == x) %>% top_frac(-0.05, pred)},
                    simplify = FALSE) %>% bind_rows(),
  .id = 'frac')

plots$within_10.0$lm$er_profile_only$fx_top_bot_dist <- labeled_ggplot(
  p = ggplot(gather(chem_env_foldx_top_bot, key = 'term', value = 'energy', total_energy:entropy_complex),
                            aes(x=mut, y=energy, colour=frac)) +
  geom_boxplot() +
  facet_wrap(~term, scales = 'free_y', ncol = 1),
  width = 10, height = 50, limitsize = FALSE
)

chem_env_foldx_top_bot_summary <- gather(chem_env_foldx_top_bot, key = 'term', value = 'energy', total_energy:entropy_complex) %>%
  group_by(frac, mut, term) %>%
  summarise(mean = mean(energy),
            sd = sd(energy),
            median = median(energy),
            upper = quantile(energy, 0.75),
            lower = quantile(energy, 0.25)) 

chem_env_foldx_top_bot_rel_summary <- gather(chem_env_foldx_top_bot_summary, key = 'sum', value = 'val', -frac, -mut, -term) %>%
  spread(key = frac, value = val) #%>%
#  mutate(top_bot_ratio = `top 5%`/`bottom 5%`) %>%
#  select(mut, term, sum, top_bot_ratio) %>%
#  mutate(top_bot_ratio = sapply(top_bot_ratio, function(x){if (is.finite(x)) x else NA})) %>%
#  spread(key = sum, value = top_bot_ratio)

plots$within_10.0$lm$er_profile_only$fx_top_bot_means <- labeled_ggplot(
  p=ggplot(filter(chem_env_foldx_top_bot_rel_summary, sum=='mean'),
                                                             aes(x=`bottom 5%`, y=`top 5%`, label=mut)) +
  geom_text() +
  facet_wrap(~term, scales = 'free') +
  geom_abline(slope=1, intercept = 0, linetype='dotted') +
  ggtitle('Mean FoldX energy of terms in the top/bottom 5% of predicted substitutions for each AA (LM without sig_count)'),
width=12, height=12)

plots$within_10.0$lm$er_profile_only$fx_top_bot_medians <- labeled_ggplot(
  p=ggplot(filter(chem_env_foldx_top_bot_rel_summary, sum=='median'),
           aes(x=`bottom 5%`, y=`top 5%`, label=mut)) +
    geom_text() +
    facet_wrap(~term, scales = 'free') +
    geom_abline(slope=1, intercept = 0, linetype='dotted') +
    ggtitle('Median FoldX energy of terms in the top/bottom 5% of predicted substitutions for each AA (LM without sig_count)'),
  width=12, height=12)

########

#### Cluster chemical environment profiles ####
clust <- drop_na(chem_env_analyses$within_10.0$tbl, sig_count) %>%
  tibble_to_matrix(prof_A:prof_Y) %>%
  dist(method = 'manhattan') %>%
  hclust()
clusters <- cutree(clust, h = 25.5)

tbl <- drop_na(chem_env_analyses$within_10.0$tbl, sig_count) %>%
  mutate(cluster = as.factor(clusters))

p <- ggplot(tbl, aes(x=tSNE1, y=tSNE2, colour = cluster)) + 
  geom_point()
  
p <- gather(tbl, key = 'mut', value = 'er', A:Y) %>%
  ggplot(aes(x=cluster, y=er)) + 
  facet_wrap(~mut) +
  geom_boxplot() +
  stat_summary(fun.data = function(x){c(label=length(x), y = 0.5)}, geom = "text", colour = 'red')
########

#### Tidy plots ####
# Select subset of models using all studies and the best, worst and two intermediate models
plots$within_10.0$lm$er_profile_sig_count$poster_subset_predictions <- pull(chem_env_analyses$within_10.0$lm_er_sig, model) %>%
  lapply(augment) %>%
  set_names(pull(chem_env_analyses$within_10.0$lm_er_sig,target)) %>%
  bind_rows(.id = 'aa') %>%
  filter(aa %in% c('E', 'Y', 'N', 'C')) %>%
  mutate(aa = c(E='Glutamate', Y='Tyrosine', N='Asparagine', C='Cysteine')[aa]) %>%
  mutate(aa = factor(aa, levels = c('Glutamate', 'Tyrosine', 'Asparagine', 'Cysteine'))) %>%
  ggplot(aes(x=er, y=.fitted)) + 
      geom_point(colour='cornflowerblue', shape=20) + 
      coord_equal() +
      facet_wrap(~aa, nrow = 2, ncol = 2) +
      ylab('Predicted ER') +
      xlab('ER') +
      theme_pubclean() +
      theme(strip.background = element_blank(), strip.text = element_text(size=30),
            axis.title = element_text(size=30), axis.text = element_text(size=20),
            panel.grid.major = element_line(colour = 'lightgrey', linetype = 'dotted'))
plots$within_10.0$lm$er_profile_sig_count$poster_subset_predictions <- labeled_ggplot(p = plots$within_10.0$lm$er_profile_sig_count$poster_subset_predictions,
                                                                                      units='in', width=10, height=10)

# Summary plot for within 10 sig count lm
loadings <- select(chem_env_analyses$within_10.0$lm_er_sig, target, model, n, r.squared) %>%
  mutate(coef_df = lapply(model, tidy)) %>%
  unnest(coef_df) %>%
  mutate(term = str_remove(term, 'prof_') %>% str_replace('sig_count', '#Sig'))

p_loadings <- ggplot(loadings, aes(x=term, y=target, fill=estimate)) + 
  geom_raster() +
  scale_fill_gradient2() +
  geom_point(data = filter(loadings, p.value < 0.0001), aes(shape='p < 0.0001')) + 
  geom_point(data = filter(loadings, p.value < 0.001, p.value > 0.0001), aes(shape='p < 0.001')) + 
  geom_point(data = filter(loadings, p.value < 0.01, p.value > 0.001), aes(shape='p < 0.01')) +
  scale_shape_manual(values = c('p < 0.0001'=8, 'p < 0.001'=3, 'p < 0.01'=20)) +
  guides(shape = guide_legend(title = '', order = 2), fill = guide_colourbar(title = 'Coefficient', order = 1)) +
  labs(y = 'Substituted AA', x = 'LM Term', tag = 'A') +
  theme(axis.ticks = element_blank(), panel.background = element_blank(), legend.position = 'left',
        legend.key = element_rect(fill = 'white'), legend.text = element_text(size = 8))

p_r_squared <- ggplot(chem_env_analyses$within_10.0$lm_er_sig, aes(x=target, y=r.squared,
                                                                   label=str_c('n = ', n, ', t = ', round(statistic, digits = 1)))) + 
  geom_col(fill = 'gray87') + 
  geom_text(colour='black', y = 0 + max(chem_env_analyses$within_10.0$lm_er_sig$r.squared)/100,
            check_overlap = FALSE, hjust=0, vjust=0.5, size = 3) +
  coord_flip() +
  theme(panel.background = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
        axis.title.y = element_blank()) +
  labs(y = 'R squared', tag = 'B') +
  lims(y = c(0, 1))

plots$within_10.0$lm$er_profile_sig_count$clean_loadings <- labeled_ggplot(
  p = ggarrange(p_loadings, p_r_squared, nrow = 1, align = 'h', widths = c(2, 1), legend = 'bottom', common.legend = TRUE),
  units = 'cm', width = 20, height = 12)

#######

# Save plots
save_plot_list(plots, root='figures/6_chemical_environment/')
