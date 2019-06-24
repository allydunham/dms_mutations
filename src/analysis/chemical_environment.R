#!/usr/bin/env Rscript 
# Functions to analyse chemical environments with respect to deep mutational scanning data

# PCA of profiles
chem_env_pca <- function(tbl, var='nearest_10'){
  
  pca <- pull(tbl, !!var) %>%
    do.call(rbind, .) %>%
    set_colnames(sort(Biostrings::AA_STANDARD)) %>%
    prcomp()
  
  out_tbl <- bind_cols(select(tbl, -!!var), as_tibble(pca$x))
  return(list(profiles=out_tbl, pca=pca))
}
