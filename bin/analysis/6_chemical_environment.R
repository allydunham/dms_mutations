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
# # Uncomment and comment above to switch to cached version
# # chem_env_analyses <- readRDS('data/rdata/chemical_environment_analyses.RDS')

plots <- sapply(chem_env_analyses, function(x){x$plots}, simplify = FALSE)

#### Specific analyses ####
## Top/bottom environments
# example workflow using Lysine substitutions
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

p_avg_profile <- select(top_bot, frac, prof_A:prof_Y) %>%
  gather('aa', 'count', -frac) %>%
  mutate(aa = str_sub(aa, start = -1)) %>% 
  group_by(frac, aa) %>%
  summarise(mean=mean(count), sd=sd(count)) %>%
  ggplot(., aes(x=aa, fill=frac, y=mean)) + 
  geom_col(position = 'dodge') + 
  geom_errorbar(aes(ymin = pmax(mean - sd, 0), ymax = mean + sd), position = 'dodge')

########

# Save plots
save_plot_list(plots, root='figures/6_chemical_environment/')
