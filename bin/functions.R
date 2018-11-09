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
