#!/usr/bin/env Rscript 
# Script analysing chemical environment in proteins, in the context of our deep mutagenesis data

source('src/config.R')
source('src/analysis/chemical_environment.R')
source('src/analysis/position_profile_clustering.R')

variant_matrices <- readRDS('data/rdata/all_study_position_matrices.RDS')

chemical_environments <- readRDS('data/rdata/position_chemical_environments.RDS') %>%
  left_join(., select(variant_matrices$norm_all_variants, study, position=pos, aa=wt, sig_count, A:Y), by = c("study", "position", "aa")) %>%
  mutate_at(.vars = vars(one_of(c('nearest_10', 'within_10.0'))),
            .funs = list(reduced=~sapply(., reduce_aa_profile, simplify = FALSE)))

chem_env_profiles <- 'within_10.0' #c('nearest_10', 'nearest_10_reduced', 'within_10.0', 'within_10.0_reduced')
chem_env_profile_col_names <- list(str_c('prof_', SORTED_AA_1_CODE))
                              # list(str_c('prof_', SORTED_AA_1_CODE),
                              #      str_c('prof_', names(AA_REDUCED_CLASSES)),
                              #      str_c('prof_', SORTED_AA_1_CODE),
                              #      str_c('prof_', names(AA_REDUCED_CLASSES)))

chem_env_analyses <- mapply(function(prof_col, prof_names){analyse_chem_env_profile(chemical_environments, prof_col = !!prof_col,
                                                                                   prof_col_names = prof_names)},
                            SIMPLIFY =  FALSE, chem_env_profiles, chem_env_profile_col_names)

# saveRDS(chem_env_analyses, 'data/rdata/chemical_environment_analyses.RDS')
# 
# # Uncomment and comment above to switch to cached version <- seems to be even slower to load...
# # chem_env_analyses <- readRDS('data/rdata/chemical_environment_analyses.RDS')

plots <- sapply(chem_env_analyses, function(x){x$plots}, simplify = FALSE)

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
plots$within_10.0$lm$top_bot_K_sub_profiles <- select(top_bot, frac, prof_A:prof_Y) %>%
  gather('aa', 'count', -frac) %>%
  mutate(aa = str_sub(aa, start = -1)) %>% 
  ggplot(., aes(x=aa, y=count, colour=frac)) + 
  geom_count(position = position_dodge(0.75)) +
  scale_fill_manual(values = c(`bottom 5%`='red', `top 5%`='blue'))

## Calc correlation between AA counts in profiles
top_5_per_cor <- tibble_to_matrix(top_5_percent, prof_A:prof_Y) %>%
  set_colnames(str_sub(colnames(.), start = -1)) %>%
  cor() %>%
  as_tibble(rownames = 'AA1') %>%
  gather(key = 'AA2', value = 'cor', -AA1)

bot_5_per_cor <- tibble_to_matrix(bot_5_percent, prof_A:prof_Y) %>%
  set_colnames(str_sub(colnames(.), start = -1)) %>%
  cor() %>%
  as_tibble(rownames = 'AA1') %>%
  gather(key = 'AA2', value = 'cor', -AA1)

# Correlation between AA pairs in profiles for top/bot K subs
top_bot_cors <- bind_rows(`top 5%`=top_5_per_cor, `bottom 5%`=bot_5_per_cor, .id = 'frac')
p_top_bot_prof_cors <- ggplot(top_bot_cors, aes(x=AA1, y=AA2, fill=cor)) +
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
  mutate(rel_cor = `top 5%` - `bottom 5%`)

top_bot_rel_cor_mat <- select(top_bot_rel_cors, AA1, AA2, rel_cor) %>%
  spread(key = AA2, value = rel_cor, ) %>%
  tibble_to_matrix(A:Y, row_names = .$AA1)

aa_order <- rownames(top_bot_rel_cor_mat)[hclust(dist(top_bot_rel_cor_mat))$order]

top_bot_rel_cors <- mutate(top_bot_rel_cors,
                              AA1 = factor(AA1, levels = aa_order),
                              AA2 = factor(AA2, levels = aa_order))

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
lm_preds_tbl <- filter(chem_env_analyses$within_10.0$lm_sig_count, study=='ALL') %>%
  pull(model) %>%
  sapply(function(x){predict(x, newdata=chem_env_analyses$within_10.0$tbl)}) %>%
  set_colnames(str_c('sig_lm_pred_',
                     filter(chem_env_analyses$within_10.0$lm_sig_count, study=='ALL')$target)) %>%
  as_tibble() %>%
  bind_cols(chem_env_analyses$within_10.0$tbl, .)

# get top/bot for K sub pred

# analyse

########

#### Tidy plots ####
# Select subset of models using all studies and the best, worst and two intermediate models
lm_preds_subset <- filter(chem_env_analyses$within_10.0$lm_sig_count, study=='ALL') %>%
  pull(model) %>%
  lapply(augment) %>%
  set_names(filter(chem_env_analyses$within_10.0$lm_sig_count, study=='ALL') %>% pull(target)) %>%
  bind_rows(.id = 'aa') %>%
  filter(aa %in% c('E', 'Y', 'N', 'C')) %>%
  mutate(aa = c(E='Glutamate', Y='Tyrosine', N='Asparagine', C='Cysteine')[aa]) %>%
  mutate(aa = factor(aa, levels = c('Glutamate', 'Tyrosine', 'Asparagine', 'Cysteine')))

plots$within_10.0$lm$subset_lm_sig_count_predictions <- labeled_ggplot(
  p = ggplot(lm_preds_subset, aes(x=er, y=.fitted)) + 
    geom_point(colour='cornflowerblue', shape=20) + 
    coord_equal() +
    facet_wrap(~aa, nrow = 2, ncol = 2) +
    ylab('Predicted ER') +
    xlab('ER') +
    theme_pubclean() +
    theme(strip.background = element_blank(), strip.text = element_text(size=30),
          axis.title = element_text(size=30), axis.text = element_text(size=20),
          panel.grid.major = element_line(colour = 'lightgrey', linetype = 'dotted')),
  width=10, height=10)

#######

# Save plots
save_plot_list(plots, root='figures/6_chemical_environment/')
