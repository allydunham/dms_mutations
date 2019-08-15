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

alphabetise_matrix <- function(x, by=c('rows', 'columns', 'both')){
  by <- match.arg(by)
  
  if (by %in% c('rows', 'both')){
    x <- x[sort(rownames(x)),]
  }
  
  if (by %in% c('columns', 'both')){
    x <- x[,sort(colnames(x))]
  }
  
  return(x)
}
########

#### Tibbles ####
# Convert a subset of a tibble to a matrix with labeled rows
tibble_to_matrix <- function(tbl, ..., row_names=NULL){
  cols <- enquos(...)
  
  if (!is.null(row_names)){
    if (length(row_names) == 1 & is.character(row_names)){
      row_names <- tbl[[row_names]]
    }
  }
  
  tbl <- select(tbl, !!! cols) %>%
    as.matrix()
  if (!is.null(row_names)){
    tbl <- set_rownames(tbl, row_names)
  }
  return(tbl)
}

# Transpose tibble
transpose_tibble <- function(tbl, col, name_col = 'rows'){
  col <- enquo(col)
  
  tibble_to_matrix(tbl, -!!col, row_names = pull(tbl, !!col)) %>%
    t() %>%
    as_tibble(rownames = name_col) %>%
    return()
}

# Count the number of unique values in each column of a tibble
col_unique_counts <- function(tbl){
  return(apply(tbl, 2, function(x){length(unique(x))}))
}

# Calc PCA based on selected columns from a tibble
tibble_pca <- function(tbl, ...){
  cols <- enquos(...)
  
  pca <- tibble_to_matrix(tbl, !!! cols) %>%
    prcomp()
  
  return(pca)
}

# Calculate correlation of tibble columns
tibble_correlation <- function(tbl, ..., filter_diag=FALSE){
  prof_cols <- enquos(...)
  
  cor_mat <- tibble_to_matrix(tbl, !!! prof_cols) %>%
    cor()
  
  group_order <- rownames(cor_mat)[hclust(dist(cor_mat))$order]
  
  cors <- as_tibble(cor_mat, rownames = 'cat1') %>%
    gather(key = 'cat2', value = 'cor', -cat1) %>%
    mutate(cat1 = factor(cat1, levels = group_order),
           cat2 = factor(cat2, levels = group_order))
  
  if (filter_diag){
    cors[cors$cat1 == cors$cat2, 'cor'] <- NA
  }
  
  return(cors)
}

# Determine order of long factor cols based on a third col (i.e. assume it is a gathered version of a relationship matrix)
add_factor_order <- function(tbl, col1, col2, value, sym=TRUE){
  col1 <- enquo(col1)
  col2 <- enquo(col2)
  value <- enquo(value)
  
  mat <- select(tbl, !!col1, !!col2, !!value) %>%
    spread(key = !!col2, value = !!value) %>%
    tibble_to_matrix(-!!col1, row_names = pull(., !!col1))
  
  if (sym){
    order1 <- rownames(mat)[hclust(dist(mat))$order]
    order2 <- order1
  } else {
    mat2 <- t(mat)
    
    order1 <- rownames(mat)[hclust(dist(mat))$order]
    order2 <- rownames(mat2)[hclust(dist(mat2))$order]
  }
  
  return(mutate(tbl,
                !!col1 := factor(!!col1, levels=order1),
                !!col2 := factor(!!col2, levels=order2))
  )
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

# Guess a vaguely sensible number of rows and columns for plotting multiple plots
get_good_rows_cols <- function(n_plots){
  for (i in 5:1){
    if (n_plots %% i == 0){
      return(c(n_plots/i, i))
    }
  }
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

#### Misc ####
# Convert log values into axis labels
# force arg forces labels to be in 'std' (..., 0.1, 1, 10, ...) or 'exp' (..., 10^-1, 1, 10^1, ...) rather than choosen 'inteligently'
# If force != 'std'/'exp' exponents below exp_notation threshold are displayed in full
# Scientific switches between 10^x and 1E0 format
make_log_labeler <- function(base=10, force='none', exp_notation_threshold=3, scientific=FALSE){
  f <- function(x){
    if ((!force == 'exp' & max(abs(x), na.rm = TRUE) < 3) | force == 'std'){
      return(sapply(base^x, format, scientific=FALSE, trim=TRUE))
    } else if (scientific) {
      return(format(base^x, scientific = TRUE, trim = TRUE))
    } else {
      return(sapply(x, function(n){if (!is.na(n) & n == 0) {1} else {bquote(.(base)^.(n))}}))
    }
  }
  return(f)
}
########
