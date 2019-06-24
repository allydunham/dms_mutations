#!/usr/bin/env Rscript 
# Functions to analyse secondary structure with respect to deep mutational scanning data

plot_secondary_structure_profile <- function(tbl, a_helix_propensity=NULL){
  # AH Positional profile
  ah_pos_avg <- group_by(tbl, alpha_helix_position) %>%
    summarise_at(.vars = vars(A:Y), .funs = mean, na.rm=TRUE) %>%
    gather(key = 'mut', value = 'score', -alpha_helix_position)
  
  p_ah_pos_profile <- ggplot(filter(ah_pos_avg, alpha_helix_position < 20),
                             aes(x=alpha_helix_position, y=mut, fill=score)) +
    geom_tile() +
    scale_fill_gradient2()
  
  # AH Propensity Plot
  if (!is.null(a_helix_propensity)){
    ah_avg_profile <- filter(tbl, !is.na(alpha_helix)) %>%
      select(A:Y) %>%
      gather(key = 'mut', value = 'score') %>%
      group_by(mut) %>%
      summarise(mean = mean(score, na.rm=TRUE),
                sd = sd(score, na.rm = TRUE),
                n = sum(!is.na(score))) %>%
      left_join(., select(a_helix_propensity, mut = aa1, exptl), by='mut')
    
    p_ah_propensity <- ggplot(ah_avg_profile, aes(x = mean, y = exptl, label = mut)) +
      geom_text() +
      geom_smooth(method = 'lm') +
      xlab('Avg Mutant Enrichment Score') + 
      ylab('Experimental Helix Propensity')
  } else {
    p_ah_propensity <- NULL
  }
  
  # AH Substitution Matrix
  ah_mat <- filter(tbl, !is.na(alpha_helix)) %>%
    gather(key = 'mut', value = 'score', A:Y) %>%
    select(study, pos, wt, mut, score) %>%
    group_by(wt, mut) %>%
    summarise(score = mean(score, na.rm=TRUE))
  
  p_ah_subs_mat <- ggplot(ah_mat, aes(x=wt, y=mut, fill=score)) +
    geom_tile() +
    scale_fill_gradient2()
  
  # BS Positional Profile
  bs_mat <- group_by(tbl, beta_sheet_position) %>%
    summarise_at(.vars = vars(A:Y), .funs = mean, na.rm=TRUE) %>%
    gather(key = 'mut', value = 'score', -beta_sheet_position)
  
  
  p_bs_pos_profile <- ggplot(filter(bs_mat, beta_sheet_position < 9), aes(x=beta_sheet_position, y=mut, fill=score)) +
    geom_tile() +
    scale_fill_gradient2()
  
  # BS Substitution Matrix
  bs_mat <- filter(tbl, !is.na(beta_sheet)) %>%
    gather(key = 'mut', value = 'score', A:Y) %>%
    select(study, pos, wt, mut, score) %>%
    group_by(wt, mut) %>%
    summarise(score = mean(score, na.rm=TRUE))
  
  p_bs_subs_mat <- ggplot(bs_mat, aes(x=wt, y=mut, fill=score)) +
    geom_tile() +
    scale_fill_gradient2()
  
  return(list(alpha_helix_profile = p_ah_pos_profile,
              alpha_helix_propensity = p_ah_propensity,
              alpha_helix_substitution_matrix = p_ah_subs_mat,
              beta_sheet_profile = p_bs_pos_profile,
              beta_sheet_substitution_matrix = p_bs_subs_mat))
}

plot_alpha_helix_dist_plots <- function(tbl, seq_intervals=5, cols=quo(A:Y)){
  ah_tbl <- select(tbl, alpha_helix, alpha_helix_position, pos, wt, !!cols) %>%
    drop_na(alpha_helix) %>%
    mutate(alpha_helix_angle = mod((alpha_helix_position - 1) * 100, 360)) %>%
    group_by(alpha_helix) %>%
    do(get_ah_distances(., cols = cols))
  
  cut_size <- max(ah_tbl$seq_dist) / seq_intervals
  ah_tbl <- mutate(ah_tbl,
                   seq_dist_cut = cut(seq_dist, breaks = seq_intervals, labels = FALSE)*cut_size-cut_size/2)
  
  p_angle <- ggplot(ah_tbl, aes(x=angle_dist, y=profile_dist)) +
    geom_violin(aes(group=angle_dist, fill=..n..)) +
    geom_smooth(method = 'lm') +
    xlab('Angular Distance') +
    ylab('Enrichment Profile Distance')
  
  p_seq_point <- ggplot(ah_tbl, aes(x=seq_dist, y=profile_dist, colour=angle_dist)) +
    geom_point() +
    geom_smooth(method = 'lm') +
    xlab('Sequence Distance') +
    ylab('Enrichment Profile Distance')
  
  p_seq_viol <- ggplot(ah_tbl, aes(x=seq_dist, y=profile_dist)) +
    geom_violin(aes(x=seq_dist_cut, group=seq_dist_cut, fill=..n..), scale = 'width') +
    geom_smooth(method = 'lm') +
    xlab('Sequence Distance') +
    ylab('Enrichment Profile Distance')
  
  return(list(angular_dist_vs_profile_dist = p_angle,
              seq_dist_vs_profile_dist = p_seq_point,
              seq_dist_vs_profile_dist_viol = p_seq_viol))
}

# Calculate distances between enrichment profiles and angular distances for all positions within an a-helix
get_ah_distances <- function(tbl, cols=quo(A:Y)){
  prof_dists <- select(tbl, !!cols) %>%
    as.matrix() %>%
    dist(method = 'manhattan')
  
  angle_dists <- select(tbl, alpha_helix_angle) %>%
    as.matrix() %>%
    dist(method = 'manhattan')
  
  helix_len <- dim(tbl)[1]
  dists <- tibble(pos1 = rep(1:helix_len, helix_len),
                  pos2 = rep(1:helix_len, each=helix_len)) %>%
    filter(pos1 != pos2, pos1 > pos2) %>%
    mutate(profile_dist = prof_dists[helix_len*(pos2-1) - pos2*(pos2-1)/2 + pos1 - pos2],
           angle_dist = angle_dists[helix_len*(pos2-1) - pos2*(pos2-1)/2 + pos1 - pos2],
           angle_dist = pmin(angle_dist, abs(360 - angle_dist)),
           seq_dist = pos1 - pos2)
  
  return(dists)
}

plot_beta_sheet_orientation <- function(tbl, cols=quo(A:Y)){
  bs_tbl <- filter(tbl, !is.na(beta_sheet)) %>%
    select(study, pos, wt, beta_sheet, beta_sheet_position, !!cols) %>%
    group_by(beta_sheet) %>%
    do(get_bs_cors(., cols = cols))
  
  p_bs_side_cor <- ggplot(bs_tbl, aes(x = same_side, y = cor)) +
    geom_boxplot()
  
  return(p_bs_side_cor)
}

get_bs_cors <- function(tbl, cols=quo(A:Y)){
  select(tbl, !!cols) %>%
    as.matrix() %>%
    t() %>%
    cor(use = 'pairwise.complete.obs') %>%
    set_colnames(str_c('V', 1:dim(.)[2])) %>%
    as_tibble(.name_repair = 'unique') %>%
    mutate(pos1 = 1:dim(.)[1]) %>%
    gather(key = 'pos2', value = 'cor', -pos1) %>%
    mutate(pos2 = as.integer(str_sub(pos2, 2)),
           side1 = mod(pos1, 2),
           side2 = mod(pos2, 2),
           same_side = side1 == side2) %>%
    filter(pos2 < pos1) %>%
    return()
}

# PCA of residues in AH/BS structures
plot_sec_struct_pca <- function(tbl){
  ah_pca <- filter(tbl, !is.na(alpha_helix)) %>%
    positional_profile_PCA()
  
  bs_pca <- filter(tbl, !is.na(beta_sheet)) %>%
    positional_profile_PCA()
  
  plots <- list()
  plots$ah_pcs <- plot_all_pcs(ah_pca$variants, colour_var = 'sig_count')
  plots$bs_pcs <- plot_all_pcs(bs_pca$variants, colour_var = 'sig_count')
  plots$bs_orientation <- plot_beta_sheet_orientation(bs_pca$variants, cols=quo(PC1:PC20))
  plots <- c(plots, plot_alpha_helix_dist_plots(ah_pca$variants, cols=quo(PC1:PC20)))
  
  return(plots)
}

# Correlation between AA frequencies and Enrichment
plot_sec_strct_freq_enrichment_correlation <- function(tbl, overall_freqs, ah_per_pos_freqs, bs_per_pos_freqs,
                                                       max_ah_length=20, max_bs_length=6){
  overall_enrich <- select(tbl, study, pos, wt, A:Y, beta_sheet, alpha_helix) %>%
    mutate(pos_class = ifelse(is.na(beta_sheet),
                              ifelse(is.na(alpha_helix),
                                     'background',
                                     'alpha_helix'),
                              'beta_sheet')) %>%
    group_by(pos_class) %>%
    summarise_at(vars(A:Y), .funs = ~ mean(., na.rm=TRUE)) %>%
    tibble_to_matrix(., -pos_class, row_names = 'pos_class') %>% 
    t() %>% 
    as_tibble(rownames = 'aa') 
  
  background_enrichment <- overall_enrich$background # Save to normalise others
  
  overall_enrich <- mutate(overall_enrich, alpha_helix_rel = alpha_helix - background, # ER is already on a Log2 scale so minus not divide
                           beta_sheet_rel = beta_sheet - background) %>%
    select(aa, `alpha helix`=alpha_helix_rel, `beta sheet`=beta_sheet_rel) %>%
    gather(key = 'ss', value = 'enrichment', -aa) %>%
    left_join(.,
              select(overall_freqs, aa, `alpha helix`=alpha_helix_rel, `beta sheet`=beta_sheet_rel) %>%
                gather(key='ss', value='log2_freq', -aa),
              by=c("aa", "ss"))
  
  p_overall_freq_enrich <- ggplot(overall_enrich, aes(x=enrichment, y=log2_freq, colour=ss, label=aa)) +
    geom_text() +
    geom_smooth(method = 'lm') +
    stat_cor()
  
  # Per position alpha helix
  per_position_ah <- select(tbl, study, pos, wt, A:Y, alpha_helix, alpha_helix_position) %>%
    group_by(alpha_helix_position) %>%
    summarise_at(vars(A:Y), .funs = ~ mean(., na.rm=TRUE)) %>%
    filter(alpha_helix_position <= max_ah_length) %>%
    tibble_to_matrix(., -alpha_helix_position) %>%
    apply(1, function(x){x - background_enrichment}) %>%
    t() %>%
    as_tibble() %>%
    mutate(alpha_helix_position = 1:max_ah_length) %>%
    gather(key = 'aa', value = 'enrichment', -alpha_helix_position) %>%
    left_join(., ah_per_pos_freqs, by = c('alpha_helix_position', 'aa'))
  
  p_per_pos_ah_freq_enrich <- ggplot(per_position_ah, aes(x=enrichment, y=log2_freq, label=aa)) +
    facet_wrap(~alpha_helix_position) +
    geom_text() +
    geom_smooth(method = 'lm') +
    stat_cor()
  
  # Per position beta_sheet
  per_position_bs <- select(tbl, study, pos, wt, A:Y, beta_sheet, beta_sheet_position) %>%
    group_by(beta_sheet_position) %>%
    summarise_at(vars(A:Y), .funs = ~ mean(., na.rm=TRUE)) %>%
    filter(beta_sheet_position <= max_bs_length) %>%
    tibble_to_matrix(., -beta_sheet_position) %>%
    apply(1, function(x){x - background_enrichment}) %>%
    t() %>%
    as_tibble() %>%
    mutate(beta_sheet_position = 1:max_bs_length) %>%
    gather(key = 'aa', value = 'enrichment', -beta_sheet_position) %>%
    left_join(., bs_per_pos_freqs, by = c('beta_sheet_position', 'aa'))
  
  p_per_pos_bs_freq_enrich <- ggplot(per_position_bs, aes(x=enrichment, y=log2_freq, label=aa)) +
    facet_wrap(~beta_sheet_position) +
    geom_text() +
    geom_smooth(method = 'lm') +
    stat_cor()
  
  return(list(overall_enrichment_frequency_cor = p_overall_freq_enrich,
              ah_per_pos_enrichment_frequency_cor = p_per_pos_ah_freq_enrich,
              bs_per_pos_enrichment_frequency_cor = p_per_pos_bs_freq_enrich))
}

# Label secondary structure on variant profile matrices
label_secondary_structure <- function(tbl, ss_col='sec_struct', min_ah_length=6, min_bs_length=4){
  tbl <- group_by(tbl, study) %>%
    mutate(region = split_protein_regions(pos)) %>%
    group_by(study, region) %>%
    mutate(beta_sheet = find_secondary_structures(.data[[ss_col]], target='E', min_length = min_bs_length,
                                                  prefix = str_c(first(study), '_', first(region), '_')),
           alpha_helix = find_secondary_structures(.data[[ss_col]], target='H', min_length = min_ah_length,
                                                   prefix = str_c(first(study), '_', first(region), '_'))) %>%
    group_by(beta_sheet) %>%
    mutate(beta_sheet_position = 1:n()) %>%
    group_by(alpha_helix) %>%
    mutate(alpha_helix_position = 1:n()) %>%
    ungroup()
  tbl$alpha_helix_position[is.na(tbl$alpha_helix)] <- NA
  tbl$beta_sheet_position[is.na(tbl$beta_sheet)] <- NA
  
  return(tbl)
}

# Find simple secondary structure runs
find_secondary_structures <- function(x, target='E', min_length=4, prefix=NULL){
  runs <- rle(x)
  target_runs <- runs$values == target & runs$lengths >= min_length
  
  final_runs <- rep(NA, length(target_runs))
  final_runs[target_runs] <- 1:sum(target_runs)
  final_runs <- rep(final_runs, times = runs$lengths)
  
  if (!is.null(prefix)){
    final_runs <- str_c(prefix, final_runs)
  }
  
  return(final_runs)
}

# Find continuous regions with data
split_protein_regions <- function(pos){
  reg_ends <- which(c(diff(pos) > 1, TRUE))
  reg_labels <- rep(1:length(reg_ends), times=c(reg_ends[1], diff(reg_ends)))
  return(reg_labels)
}

