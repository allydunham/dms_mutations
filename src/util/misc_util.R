#!/usr/bin/env Rscript 
# Script containing small utility functions used throughout the project

#### raw computation ####
# returns a list with a vector of unique row labels and a vector of bools as from duplicated(mat)
# combined these can map rowwise results from de-duplicated matrix back to original as int labels correspond to
# row numbers in de-duped matrix
enumerate_unique_rows <- function(mat){
  dupe_rows <- duplicated(mat)
  
  unique_inds <- rep(0, length(dupe_rows))
  ind <- 1
  for (i in 1:length(dupe_rows)){
    if (!dupe_rows[i]){
      # enumerate non-dupe rows with unused index
      unique_inds[i] <- ind
      ind <- ind + 1
      
    } else if (unique_inds[i] == 0) {
      # If its a dupe row, set indeces to first occurance (which will have been previously set)
      dupe_row_nums <- which(apply(mat, 1, identical, y=mat[i,]))
      unique_inds[dupe_row_nums] <- unique_inds[dupe_row_nums[1]]
    }
  }
  return(list(indeces=unique_inds, duplicate=dupe_rows))
}

# Given a matrix return ind
find_matching_rows <- function(target, mat){
  apply(mat, 1, identical, y=target)
}

########

#### Tibbles ####
# Convert a subset of a tibble to a matrix with labeled rows
tibble_to_matrix <- function(tbl, columns, row_names=NULL){
  columns <- enquo(columns)
  if (!is.null(row_names)){
    if (length(row_names) == 1 & is.character(row_names)){
      row_names <- tbl[[row_names]]
    }
  }
  tbl <- select(tbl, !!columns) %>%
    as.matrix()
  if (!is.null(row_names)){
    tbl <- set_rownames(tbl, row_names)
  }
  return(tbl)
}

# Count the number of unique values in each column of a tibble
col_unique_counts <- function(tbl){
  return(apply(tbl, 2, function(x){length(unique(x))}))
}

########

#### Proteins ####
# Convert AA encoding between 1 letter/3 letter
AA_CODE_1_TO_3 <- AMINO_ACID_CODE
AA_CODE_3_TO_1 <- structure(names(AMINO_ACID_CODE), names=AMINO_ACID_CODE)

aa_1_to_3 <- function(x){
  return(unname(AA_CODE_1_TO_3[x]))
}

aa_3_to_1 <- function(x){
  return(unname(AA_CODE_3_TO_1[x]))
}
########

#### Plotting ####
# Plot all PCs
plot_all_pcs <- function(tbl, max_pc=20, colour_var='wt', nrow=2, ncol=5, width = NULL, height = NULL, geom=geom_point()){
  if (nrow * ncol * 2 != max_pc){
    stop('Invalid row/col numbers: nrow * ncol * 2 != max_pc')
  }
  return(
    labeled_ggplot(
      ggarrange(
        plotlist = lapply(
          seq(1, max_pc-1, 2),
          function(x){
            ggplot(tbl, aes_string(x=str_c('PC', x), y=str_c('PC', x + 1), colour=colour_var)) + 
              geom
          }),
        nrow = nrow, ncol = ncol, common.legend = TRUE, legend = 'right'
      ),
      width = ifelse(is.null(width), ncol * 4, width), height = ifelse(is.null(height), nrow * 4, height)
    )
  )
}

# Plot scatter plots against a list of factors (supplied as symbols)
plot_factors <- function(tbl, x, y, factors){
  x = enquo(x)
  y = enquo(y)
  plist <- sapply(factors,
                  function(fac){ggplot(tbl, aes(x=!!x, y=!!y, colour=!!fac)) + geom_point()},
                  simplify = FALSE)
  return(plist)
}
########