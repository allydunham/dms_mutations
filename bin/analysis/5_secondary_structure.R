#!/usr/bin/env Rscript 
# Script analysing deep mutagenesis data in the context of secondary structure

source('src/config.R')
source('src/analysis/secondary_structure.R')
source('src/analysis/position_profile_clustering.R')

a_helix_propensity <- readRDS('data/rdata/a_helix_propensity.RDS')
secondary_structure <- readRDS('data/rdata/position_secondary_structure.RDS')
variant_matrices <- readRDS('data/rdata/all_study_position_matrices.RDS')
imputed_matrices <- readRDS('data/rdata/all_study_imputed_position_matrices.RDS')

plots <- list()

#### AA frequencies in secondary structure ####
labeled_secondary_structures <- label_secondary_structure(secondary_structure, ss_col = 'ss')

## Overall
ss_aa_frequencies <- mutate(labeled_secondary_structures, pos_class = ifelse(is.na(beta_sheet),
                                                                             ifelse(is.na(alpha_helix),
                                                                                    'background',
                                                                                    'alpha_helix'),
                                                                             'beta_sheet')) %>%
  group_by(pos_class) %>%
  group_modify(.f = ~ tibble(!!!table(.$aa)/sum(table(.$aa)))) %>%
  ungroup() %>%
  tibble_to_matrix(., -pos_class, row_names = 'pos_class') %>% 
  t() %>% 
  as_tibble(rownames = 'aa') %>%
  mutate(alpha_helix_rel = log2(alpha_helix/background),
         beta_sheet_rel = log2(beta_sheet/background))

plots$ss_aa_freqs <- ggplot(select(ss_aa_frequencies, -alpha_helix_rel, -beta_sheet_rel) %>%
                              gather(key = 'ss', value = 'freq', -aa),
                            aes(x=aa, y=ss, fill=freq)) +
  geom_tile() +
  theme(axis.ticks = element_blank(), panel.background = element_blank())

plots$rel_ss_aa_freqs <- ggplot(select(ss_aa_frequencies, aa, alpha_helix=alpha_helix_rel, beta_sheet=beta_sheet_rel) %>%
                                  gather(key = 'ss', value = 'log2_freq', -aa) %>%
                                  filter(!aa == 'Z'),
                                aes(x=aa, y=ss, fill=log2_freq)) +
  geom_tile() +
  scale_fill_gradient2() +
  theme(axis.ticks = element_blank(), panel.background = element_blank())


## Per position alpha helix frequencies
max_common_alpha_helix_length <- 20
per_position_ah_aa_freqs <- filter(labeled_secondary_structures, !is.na(alpha_helix)) %>%
  select(study, pos, aa, alpha_helix, alpha_helix_position) %>%
  group_by(alpha_helix_position) %>%
  group_modify(.f = ~ tibble(!!!table(.$aa)/sum(table(.$aa)))) %>%
  ungroup() %>%
  select(alpha_helix_position, !!!sort(Biostrings::AA_STANDARD)) %>%
  filter(alpha_helix_position <= max_common_alpha_helix_length) %>% # Don't have many examples after this
  tibble_to_matrix(., -alpha_helix_position) %>%
  apply(1, function(x){log2(x/ss_aa_frequencies$background[1:20])}) %>% # 1:20 removes Z (which is tiny proportion anyway)
  t() %>%
  as_tibble() %>%
  mutate(alpha_helix_position = 1:max_common_alpha_helix_length) %>%
  gather(key = 'aa', value = 'log2_freq', -alpha_helix_position)

plots$rel_alpha_helix_pos_aa_freqs <- ggplot(per_position_ah_aa_freqs,
                                             aes(x=alpha_helix_position, y=aa, fill=log2_freq)) +
  geom_tile() + 
  scale_fill_gradient2() +
  theme(axis.ticks = element_blank(), panel.background = element_blank())

## Per position beta_sheet frequencies
max_common_beta_sheet_length <- 8
per_position_bs_aa_freqs <- filter(labeled_secondary_structures, !is.na(beta_sheet)) %>%
  select(study, pos, aa, beta_sheet, beta_sheet_position) %>%
  group_by(beta_sheet_position) %>%
  group_modify(.f = ~ tibble(!!!table(.$aa)/sum(table(.$aa)))) %>%
  ungroup() %>%
  select(beta_sheet_position, !!!sort(Biostrings::AA_STANDARD)) %>%
  filter(beta_sheet_position <= max_common_beta_sheet_length) %>% # Don't have many examples after this
  tibble_to_matrix(., -beta_sheet_position) %>%
  apply(1, function(x){log2(x/ss_aa_frequencies$background[1:20])}) %>% # 1:20 removes Z (which is tiny proportion anyway)
  t() %>%
  as_tibble() %>%
  mutate(beta_sheet_position = 1:max_common_beta_sheet_length) %>%
  gather(key = 'aa', value = 'log2_freq', -beta_sheet_position)

plots$rel_beta_sheet_pos_aa_freqs <- ggplot(per_position_bs_aa_freqs,
                                            aes(x=beta_sheet_position, y=aa, fill=log2_freq)) +
  geom_tile() + 
  scale_fill_gradient2() +
  theme(axis.ticks = element_blank(), panel.background = element_blank())
########

#### SS against deep mut results ####
plots <- list_modify(plots, !!!sapply(variant_matrices,
                                      plot_sec_strct_freq_enrichment_correlation,
                                      overall_freqs=ss_aa_frequencies,
                                      ah_per_pos_freqs=per_position_ah_aa_freqs,
                                      bs_per_pos_freqs=per_position_bs_aa_freqs,
                                      simplify = FALSE))

plots <- list_modify(plots, sapply(variant_matrices, plot_secondary_structure_profile,
                                   a_helix_propensity=a_helix_propensity, simplify = FALSE))

plots$alpha_helix_lengths <- ggplot(variant_matrices$all_variants %>%
                                      group_by(alpha_helix) %>%
                                      summarise(length = max(alpha_helix_position)),
                                    aes(x = length)) + geom_histogram()

plots$beta_sheet_lengths <- ggplot(variant_matrices$all_variants %>%
                                     group_by(beta_sheet) %>%
                                     summarise(length = max(beta_sheet_position)),
                                   aes(x = length)) + geom_histogram()

plots <- list_modify(plots, !!!sapply(variant_matrices, plot_alpha_helix_dist_plots, simplify = FALSE))
plots <- list_modify(plots, !!!sapply(variant_matrices, simplify = FALSE,
                                      function(x){list(beta_sheet_side_boxplot=plot_beta_sheet_orientation(x))}))
plots <- list_modify(plots, !!!sapply(imputed_matrices, plot_sec_struct_pca, simplify = FALSE))

# Save plots
save_plot_list(plots, root='figures/5_secondary_structure')
