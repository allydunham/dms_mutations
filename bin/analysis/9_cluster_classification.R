#!/usr/bin/env Rscript 
# Attempt to predict ER clusters of positions via SIFT/FoldX scores

source('src/config.R')
source('src/analysis/cluster_classification.R')

library(caret)

dm_hclust <- readRDS('data/rdata/deep_mut_hclust_clusters.RDS')

# Narrow down data to useful columns
dm_tbl <- select(dm_hclust$tbl, cluster, study, gene_name, position=pos, wt, A:Y, all_atom_rel, ss, phi, psi) %>%
  left_join(., select(dm_hclust$foldx, cluster, study, position=pos, wt, total_energy:entropy_complex),
            by = c('cluster', 'study', 'position', 'wt')) %>%
  select(-total_energy, -sloop_entropy, -mloop_entropy, -entropy_complex, -electrostatic_kon, -water_bridge) # drop unused foldx terms 

# Classify 
train_control <- trainControl(method = 'cv', number = 10, savePredictions = TRUE)
model <- train(cluster ~ ., data = filter(dm_tbl, wt == 'A') %>% select(cluster, all_atom_rel, phi:energy_ionisation) %>% drop_na(),
               method = 'rf')
