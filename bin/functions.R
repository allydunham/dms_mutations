#!/usr/bin/env Rscript 
# Script containing general functions for deep mutagenesis analysis

# generate mut id for variants
gen_mut_id <- function(acc, ref, alt, pos){
  return(paste(acc, pos, alt, sep = '_'))
}

# Find the next row of all NA values in a tbl
find_next_na_row <- function(start_row, tbl){
  r <- start_row
  final_row <- dim(tbl)[1]
  repeat{
    if (all(is.na(tbl[r,]))){
      return(r)
    } else if (r == final_row){
      warning('No more all NA rows, returning final row')
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
  sel <- as.matrix(select(sel, -position, -ref_aa))
  bkg <- as.matrix(select(bkg, -position, -ref_aa))
}

# Calculate fitness as increase using Enrichment ratio
# Expects bkg and sel as matrices with rows as loci and ref_aas as a char vector
# Currently no pseudocount etc.
er_fitness <- function(sel, bkg, ref_aas){
  ### Determine Approx WT reads
  wt_inds <- cbind(1:length(ref_aas), match(ref_aas, colnames(sel)))
  
}


