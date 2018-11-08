#!/usr/bin/env Rscript 
# Script containing general functions for deep mutagenesis analysis

# generate mut id for variants
gen_mut_id <- function(acc, ref, alt, pos){
  return(paste(acc, pos, alt, sep = '_'))
}
