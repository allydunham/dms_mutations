#!/usr/bin/env Rscript 
# Script containing small utility functions used throughout the project

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

# Convert AA encoding between 1 letter/3 letter
AA_CODE_1_TO_3 <- AMINO_ACID_CODE
AA_CODE_3_TO_1 <- structure(names(AMINO_ACID_CODE), names=AMINO_ACID_CODE)

aa_1_to_3 <- function(x){
  return(unname(AA_CODE_1_TO_3[x]))
}

aa_3_to_1 <- function(x){
  return(unname(AA_CODE_3_TO_1[x]))
}

# Plot all PCs
plot_all_pcs <- function(tbl, max_pc=20, colour_var='wt', nrow=2, ncol=5, width = NULL, height = NULL){
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
              geom_point()
          }),
        nrow = nrow, ncol = ncol, common.legend = TRUE, legend = 'right'
      ),
      width = ifelse(is.null(width), ncol * 4, width), height = ifelse(is.null(height), nrow * 4, height)
    )
  )
}
