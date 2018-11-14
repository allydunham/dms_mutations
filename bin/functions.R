#!/usr/bin/env Rscript 
# Script containing general functions for deep mutagenesis analysis

# generate mut id for variants
gen_mut_id <- function(acc, ref, alt, pos){
  if (is.na(ref)){
    ref <- 'X'
  }
  return(paste0(acc, '_', pos, '_', alt))
}

# Find the next empty row (All NA or bottom) in a tbl
find_next_empty_row <- function(start_row, tbl){
  r <- start_row
  final_row <- dim(tbl)[1]
  repeat{
    if (all(is.na(tbl[r,]))){
      return(r)
    } else if (r == final_row){
      return(r)
    } else {
      r <- r + 1
    }
  }
}

# Wrapper to pass correct background and selection counts to fitness function, based on format of Melnikov 2014 data
# Expects sel to be a data.frame with cols for position, ref_aa and all alt_aa's in one selection/drug/library category
# these are given as exp_name in the SX_DRUG_LX format of melnikov
melnikov_fitness <- function(sel, exp_name, bkg){
  # Extract meta info on experiment
  meta <- as.list(strsplit(exp_name, '_')[[1]])
  names(meta) <- c('selection_round', 'drug', 'library')
  
  # Select correct background reads for library
  bkg <- bkg[[paste0('Bkg', str_sub(meta$library, -1))]]
  
  # Format bkg and sel as matrices
  ref_aas <- bkg$ref_aa
  gene_length <- length(ref_aas)
  sel <- as.matrix(select(sel, -position, -ref_aa))
  bkg <- as.matrix(select(bkg, -position, -ref_aa))
  
  # Find WT positions and set to NA (as WT measurement is poor)
  wt_inds <- cbind(1:gene_length, match(ref_aas, colnames(sel)))
  sel[wt_inds] <- NA
  bkg[wt_inds] <- NA
  
  e_scores <- e_score(sel, bkg)
  
  # Currently don't convert E score to fitness since hard to determine number of WT seqs in library (counted many times over)
  # Could use mean of very similar AAs (e.g. high blosum pairs) as WT equivalent
  #fitness <- e_scores/mean(e_scores[wt_inds])
  fitness <- e_scores

  fitness %<>% as.tibble() %>%
    mutate(position = 1:gene_length,
           ref_aa = ref_aas) %>%
    gather(key = 'alt_aa', value = 'e_score', -ref_aa, -position)
  return(fitness)
}

# Calculate E-score equivalent to Enrich1 
# Currently no pseudocount etc. (simple implementation without error checking for prelim analysis)
e_score <- function(sel, bkg){
  bkg[bkg == 0] <- NA
  
  freq_sel <- sel/sum(sel, na.rm = TRUE)
  freq_bkg <- bkg/sum(bkg, na.rm = TRUE)
  
  return(freq_sel/freq_bkg)
}


