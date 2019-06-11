#!/usr/bin/env Rscript 
# Script containing small utility functions used throughout the project

# Convert AA encoding between 1 letter/3 letter
AA_CODE_1_TO_3 <- AMINO_ACID_CODE
AA_CODE_3_TO_1 <- structure(names(AMINO_ACID_CODE), names=AMINO_ACID_CODE)

aa_1_to_3 <- function(x){
  return(unname(AA_CODE_1_TO_3[x]))
}

aa_3_to_1 <- function(x){
  return(unname(AA_CODE_3_TO_1[x]))
}