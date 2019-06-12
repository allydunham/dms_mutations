#!/usr/bin/env Rscript 
# Script to analyse deep mutagenesis data itself

source('bin/config.R')

#### Import and process data ####
dms_data <- readRDS('data/variant_data.RDS')

deep_variant_plots <- list(score_distributions=list())

# Meta data about each study, including chosen threshold for deleterousness
meta_df <- tibble(study = names(dms_data),
                  gene_type = sapply(dms_data, function(x){get_meta(x$dm, 'gene_type')}),
                  gene_name = sapply(dms_data, function(x){get_meta(x$dm, 'gene_name')}),
                  test_class = sapply(dms_data, function(x){get_meta(x$dm, 'test_class')}),
                  species = sapply(dms_data, function(x){get_meta(x$dm, 'species')}),
                  authour = sapply(dms_data, function(x){get_meta(x$dm, 'authour')}),
                  group = sapply(dms_data, function(x){get_meta(x$dm, 'group')}),
                  single = sapply(dms_data, function(x){"single_variant" %in% class(x)}),
                  thresh=sapply(dms_data, function(x){x$manual_threshold}),
                  factor=sapply(dms_data, function(x){x$norm_factor}),
                  norm_thresh=thresh/factor)

# Helix propensity
a_helix_propensity <- read_tsv('meta/alpha_helix_energy', comment = '#') %>%
  rename_all(str_to_lower) %>%
  mutate(aa1 = aa_3_to_1(aa)) %>%
  select(aa1, aa3 = aa, exptl:luque)

# Surface accessibility scores for all variants
surface_accesibility <- sapply(dms_data, function(x){if(!identical(NA, x$surface_accesibility)){x$surface_accesibility$combined}},
                                   simplify = FALSE) %>%
  bind_rows(.id = 'study')

# Secondary structure predictions for all positions
secondary_structure <- sapply(dms_data, function(x){if(!identical(NA, x$secondary_structure)){x$secondary_structure}},
                              simplify = FALSE) %>%
  bind_rows(.id = 'study')

# Dataframe of all individually scored variant/sets of variants in all studies
all_variants <- bind_rows(lapply(dms_data, function(x){x$dm$variant_data}), .id = 'study') %>%
  select(study, variants, score, raw_score, norm_score) %>%
  mutate(variants = str_replace_all(variants, 'p.', '')) %>%
  left_join(., meta_df, by='study')

# Dataframe with matrix of mutational profiles for all positions and all studies
variant_matrices <- list()
variant_matrices$all_variants <- bind_rows(lapply(dms_data, make_var_matrix, score='score'), .id = 'study') %>%
  filter(!wt %in% c('Z', 'B')) %>%
  drop_na(wt, pos) %>%
  left_join(., meta_df, by='study') %>%
  left_join(., rename(surface_accesibility, wt=res1), by=c('study', 'wt', 'pos')) %>%
  left_join(., rename(secondary_structure, wt=aa, sec_struct=ss) %>%
              rename_at(vars(-study, -pos, -wt, -sec_struct), .funs= ~ str_c('ss_prob_', .)),
            by=c('study', 'wt', 'pos')) %>%
  select(-factor, -norm_thresh) %>%
  mutate(sig_count = mutate_at(., .vars = vars(A:Y), .funs = list(~ . < thresh)) %>%
           select(A:Y) %>%
           rowSums(na.rm = TRUE))

variant_matrices$norm_all_variants <- bind_rows(lapply(dms_data, make_var_matrix, score='norm_score'), .id = 'study') %>%
  filter(!wt %in% c('Z', 'B')) %>%
  drop_na(wt, pos) %>%
  left_join(., meta_df, by='study') %>%
  left_join(., rename(surface_accesibility, wt=res1), by=c('study', 'wt', 'pos')) %>%
  left_join(., rename(secondary_structure, wt=aa, sec_struct=ss) %>%
              rename_at(vars(-study, -pos, -wt, -sec_struct), .funs= ~ str_c('ss_prob_', .)),
            by=c('study', 'wt', 'pos')) %>%
  select(-factor, -thresh, thresh = norm_thresh) %>%
  mutate(sig_count = mutate_at(., .vars = vars(A:Y), .funs = list(~ . < thresh)) %>%
           select(A:Y) %>%
           rowSums(na.rm = TRUE))

# Label secondary structure before filtering data
variant_matrices <- sapply(variant_matrices, label_secondary_structure, simplify = FALSE)

# Filter interesting subsets
variant_matrices$sig_positions <- filter(variant_matrices$all_variants, sig_count > 0)
variant_matrices$norm_sig_positions <- filter(variant_matrices$norm_all_variants, sig_count > 0)
variant_matrices$single_variants <- filter(variant_matrices$all_variants, single)
variant_matrices$sig_single_variants <- filter(variant_matrices$sig_positions, single)
variant_matrices$norm_single_variants <- filter(variant_matrices$norm_all_variants, single)
variant_matrices$norm_sig_single_variants <- filter(variant_matrices$norm_sig_positions, single)

imputed_matrices <- sapply(variant_matrices, impute_variant_profiles, background_matrix=variant_matrices$all_variants, simplify=FALSE)

#### Enrichment Score Distributions ####
deep_variant_plots$dm_hists <- labeled_ggplot(p = plot_study_histogram(all_variants, meta_df),
                                              width = 14, height = 9)

deep_variant_plots$dm_norm_hists <- labeled_ggplot(p = plot_study_histogram(all_variants, meta_df,
                                                                            x='norm_score', thresh = 'norm_thresh'),
                                                   width = 14, height = 9)

factor_density_config <- list(gene_type=list(facet='~gene_type'), test_class=list(facet='~test_class'),
                              species=list(facet='~species'), gene=list(facet='~gene_name'))
deep_variant_plots$score_factor_densities <- sapply(factor_density_config, function(x){
  labeled_ggplot(
    do.call(plot_factor_density, c(list(tbl=all_variants), x)),
    width=14, height=9
    )
}, simplify = FALSE)

per_position_boxplots <- gather(variant_matrices$all_variants, key='mut', value='score', A:Y) %>%
  select(study, pos, wt, mut, score) %>%
  group_by(study) %>%
  do(plots = plot_score_distribution_per_position(.)) 

deep_variant_plots$score_distributions$per_position_boxplots <- per_position_boxplots$plots
names(deep_variant_plots$score_distributions$per_position_boxplots) <- per_position_boxplots$study
########

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

deep_variant_plots$blosum <- list()

deep_variant_plots$blosum$aa_profile_heatmap <- ggplot(sub_profile_scores, aes(x=wt, y=mut, fill=er)) + 
  geom_tile(colour='white') + 
  scale_fill_gradient2() +
  theme(axis.ticks = element_blank(), panel.background = element_blank())
        
deep_variant_plots$blosum$er_vs_blosum <- ggplot(sub_scores, aes(x=blosum62, y=er, group=blosum62)) + geom_boxplot()
deep_variant_plots$blosum$profile_er_vs_blosum <- ggplot(sub_profile_scores, aes(x=blosum62, y=er, group=blosum62)) + 
  geom_boxplot()

########

#### Variant Score PCAs ####
variant_pcas <- sapply(imputed_matrices, positional_profile_PCA, simplify = FALSE)
  
deep_variant_plots$pcas <- list()

# plot PCAs against basic covariates
deep_variant_plots$pcas <- sapply(variant_pcas, basic_pca_plots, simplify = FALSE)

# Average PCA loading for each AA
avg_AA_pca_profiles <- sapply(variant_pcas, get_avg_aa_pca_profile, simplify = FALSE)

deep_variant_plots <- list_modify(deep_variant_plots, pcas=sapply(avg_AA_pca_profiles, aa_avg_profile_plot, simplify = FALSE))

deep_variant_plots <- list_modify(deep_variant_plots, pcas=sapply(avg_AA_pca_profiles, aa_profile_heatmap, simplify = FALSE))

# Corelation between PCs and surface accesibility
pca_surf_cor <- sapply(variant_pcas, pca_surf_acc_cor, simplify = FALSE)

deep_variant_plots <- list_modify(deep_variant_plots, pcas=sapply(pca_surf_cor, pca_surface_heatmap, simplify = FALSE))

# Per AA PCAs #
AAs <- sort(unique(variant_matrices$all_variants$wt))
deep_variant_plots <- list_modify(deep_variant_plots, pcas=sapply(imputed_matrices, function(x){
  list(per_aa_pcas=sapply(AAs, per_aa_pcas, variant_matrix=x, simplify = FALSE))
}, simplify = FALSE))

# Predict properties from PCs



########

#### Profile vs SIFT/FoldX ####
prediction_martices <- list()
prediction_martices$sift <- bind_rows(lapply(dms_data, make_var_matrix, score='sift_score'), .id = 'study') %>%
  mutate(sift = TRUE)
prediction_martices$foldx <- bind_rows(lapply(dms_data, make_foldx_var_matrix), .id = 'study')  %>%
  mutate(foldx = TRUE)
prediction_martices$exp <- select(variant_matrices$norm_all_variants, study, pos, wt, A:Y)

pairwise_dists <- mapply(compute_per_gene_pairwise_profile_dists,
                         prediction_martices,
                         str_c(names(prediction_martices),'_dist'),
                         SIMPLIFY = FALSE) %>%
  reduce(left_join, by = c('study', 'pos1', 'aa1', 'pos2', 'aa2'))

deep_variant_plots$pred_dists <- list()
deep_variant_plots$pred_dists$foldx <- ggplot(drop_na(pairwise_dists, foldx_dist), aes(x=exp_dist, y=foldx_dist, colour=study)) + 
  geom_density_2d() +
  geom_smooth(method = 'lm', colour='black') +
  facet_wrap(~study) + 
  theme(legend.position = 'none')

deep_variant_plots$pred_dists$sift <- ggplot(drop_na(pairwise_dists, sift_dist), aes(x=exp_dist, y=sift_dist, colour=study)) + 
  geom_density_2d() +
  geom_smooth(method = 'lm', colour='black') +
  facet_wrap(~study) + 
  theme(legend.position = 'none')

########

#### Secondary structure ####
deep_variant_plots$secondary_structure <- lapply(variant_matrices, plot_secondary_structure_profile, a_helix_propensity=a_helix_propensity)
deep_variant_plots$secondary_structure$alpha_helix_lengths <- ggplot(variant_matrices$all_variants %>%
                                                              group_by(alpha_helix) %>%
                                                              summarise(length = max(alpha_helix_position)),
                                                            aes(x = length)) + geom_histogram()

deep_variant_plots$secondary_structure$beta_sheet_lengths <- ggplot(variant_matrices$all_variants %>%
                                                              group_by(beta_sheet) %>%
                                                              summarise(length = max(beta_sheet_position)),
                                                            aes(x = length)) + geom_histogram()

deep_variant_plots <- list_modify(deep_variant_plots, secondary_structure=sapply(variant_matrices, plot_alpha_helix_dist_plots, simplify = FALSE))
deep_variant_plots <- list_modify(deep_variant_plots, secondary_structure=sapply(
  variant_matrices, simplify = FALSE, function(x){list(beta_sheet_side_boxplot=plot_beta_sheet_orientation(x))}))
deep_variant_plots <- list_modify(deep_variant_plots, secondary_structure=sapply(imputed_matrices, plot_sec_struct_pca, simplify = FALSE))
  
########

#### Save plots #####
save_plot_list(deep_variant_plots, root='figures/variant_analysis/')
