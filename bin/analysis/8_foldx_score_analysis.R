#!/usr/bin/env Rscript 
# Perform complimentary analysis on FoldX scores to that done on DMS data

source('src/config.R')
source('src/analysis/foldx_score_analysis.R')
source('src/analysis/position_profile_clustering.R')

foldx <- readRDS('data/rdata/human_foldx_reduced.RDS')

foldx_ptms <- readRDS('data/rdata/human_foldx_ptms.RDS')

#### Analyse by average at position ####
# foldx_pos_avg <- group_by(foldx, uniprot_id, position, wt, modification) %>%
#   summarise_at(vars(total_energy:entropy_complex), mean, na.rm=TRUE) %>%
#   ungroup()
# saveRDS(foldx_pos_avg, 'data/rdata/human_foldx_reduced_pos_avg.RDS')
foldx_pos_avg <- readRDS('data/rdata/human_foldx_reduced_pos_avg.RDS')

########