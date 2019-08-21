#!/usr/bin/env Rscript 
# Attempt to predict ER clusters of positions via SIFT/FoldX scores

source('src/config.R')
source('src/analysis/cluster_classification.R')

library(caret)

