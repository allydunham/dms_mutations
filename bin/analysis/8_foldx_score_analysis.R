#!/usr/bin/env Rscript 
# Perform complimentary analysis on FoldX scores to that done on DMS data

source('src/config.R')
source('src/analysis/foldx_score_analysis.R')
source('src/analysis/position_profile_clustering.R')

foldx <- readRDS('data/rdata/human_foldx_tiny.RDS')

plots <- list()

#### Analyse by average at position ####
# foldx_pos_avg <- group_by(foldx, uniprot_id, position, wt) %>%
#   summarise_at(vars(total_energy:entropy_complex), mean, na.rm=TRUE) %>%
#   ungroup()
# saveRDS(foldx_pos_avg, 'data/rdata/human_foldx_tiny_pos_avg.RDS')
foldx_pos_avg <- readRDS('data/rdata/human_foldx_tiny_pos_avg.RDS')

# drop unused factors then center and scale (pssobly add water_bridge/kon back if structures where they're meaningful are used later)
foldx_pos_avg_scaled <- select(foldx_pos_avg, -sloop_entropy, -mloop_entropy, -entropy_complex, -electrostatic_kon, -water_bridge) %>% 
  tibble_to_matrix(., total_energy:energy_ionisation) %>% 
  scale(center = FALSE, scale = TRUE) %>% 
  as_tibble() %>%
  bind_cols(select(foldx_pos_avg, uniprot_id:modification), .)

pos_avg_settings <- list(h=2.5, k=NULL, max_k=6)
pos_avg_hclust <- group_by(foldx_pos_avg_scaled, wt) %>%
  do(hclust = make_hclust_clusters(., total_energy:energy_ionisation, h = pos_avg_settings$h, max_k = pos_avg_settings$max_k))
#sapply(pos_avg_hclust$hclust, function(x){plot(x$hclust); abline(h = 12)})

pos_avg_hclust_tbl <- map_dfr(pos_avg_hclust$hclust, .f = ~ .[[1]]) %>%
  mutate(cluster = str_c(wt, '_', cluster))

pos_avg_hclust_analysis <- analyse_clusters(pos_avg_hclust_tbl,
                                            str_c('hclust (h = ', pos_avg_settings$h, 
                                                  ', k = ', pos_avg_settings$k, 
                                                  ', max_k = ', pos_avg_settings$max_k, ')'),
                                            total_energy:energy_ionisation,
                                            transform_ddg = function(x){x / max(abs(x))})

plots$pos_avg_clusters <- pos_avg_hclust_analysis$plots
########

#### Analyse by total energy of each substitution ####
foldx_total_energy <- select(foldx, uniprot_id:mut, total_energy) %>%
  spread(key = mut, value = total_energy)

tot_eng_settings <- list(h=75, k=NULL, max_k=6)
tot_eng_hclust <- group_by(foldx_total_energy, wt) %>%
  do(hclust = make_hclust_clusters(., A:Y, h = tot_eng_settings$h, max_k = tot_eng_settings$max_k))
#sapply(tot_eng_hclust$hclust, function(x){plot(x$hclust); abline(h = tot_eng_settings$h)})

tot_eng_hclust_tbl <- map_dfr(tot_eng_hclust$hclust, .f = ~ .[[1]]) %>%
  mutate(cluster = str_c(wt, '_', cluster))

tot_eng_hclust_analysis <- analyse_clusters(tot_eng_hclust_tbl,
                                            str_c('hclust (h = ', tot_eng_settings$h, 
                                                  ', k = ', tot_eng_settings$k, 
                                                  ', max_k = ', tot_eng_settings$max_k, ')'),
                                            A:Y,
                                            transform_ddg = function(x){atan(0.5 * x)})
tot_eng_hclust_analysis$plots$mean_profile$plot <- tot_eng_hclust_analysis$plots$mean_profile$plot + 
  theme(axis.text.x = element_text(angle = 0))
  scale_fill_gradientn(colors = c('red', 'white', 'blue', 'yellow', 'green'),
                       values = rescale(c(min(tot_eng_hclust_analysis$mean_profiles_long$trans_ddg),
                                          0 ,
                                          -1 * min(tot_eng_hclust_analysis$mean_profiles_long$trans_ddg),
                                          30,
                                          max(tot_eng_hclust_analysis$mean_profiles_long$trans_ddg))))

plots$total_energy <- tot_eng_hclust_analysis$plots
########

# Save Plots
save_plot_list(plots, root='figures/8_foldx_score_analysis/')
