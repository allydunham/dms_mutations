#!/usr/bin/env Rscript 
# Script containing general functions for deep mutagenesis analysis

#### Formatting and Storage ####
# generate mut id for variants
gen_mut_id <- function(acc, ref, alt, pos){
  if (is.na(ref)){
    ref <- 'X'
  }
  return(paste0(acc, '_', pos, '_', alt))
}

## DeepMut data 'Class'
# variant_data = data frame with variants (requires: protein variant string (comma separated list of variants in hgvs protein format),
# score & raw_score, supports other orig columns too)
# gene_name = name string
# alt_name = alternative gene name
# accessions = named chr vector of accessions
# study = named chr vector of paper metadata
# species = species string
# aa_seq = reference AA seq string
# authour = paper authours
# year = paper year (minimal info to locate, best to include more (such as url) in misc)
# transform = method used to transform score
# misc = named list of other meta data to attach (e.g. other protein accessions, other study data.
# Various special misc values are expected and are written in a logical positon - alt_name, doi, pmid, url, title
# anything else with _id is treated as a gene id
DeepMut <- function(variant_data, ...){
  # (minimal) Error Checking
  if (!all(c('variants', 'score', 'raw_score') %in% colnames(variant_data))) {
    stop('variant_data does not contain the required columns (protein_variants, score, raw_score)')
  }
  
  
  defaults <- list(gene_name=NA, domain=NA, species=NA, aa_seq=NA, transform='None', uniprot_id=NA, authour=NA, year=NA)

  deep_mut <- list(variant_data=variant_data, ...)
  deep_mut <- c(deep_mut, defaults[!names(defaults) %in% names(deep_mut)])

  class(deep_mut) <- c('DeepMut', class(deep_mut))
  return(deep_mut)
}

# Class for grounping similar DeepMut datasets
# lst - named list of DeepMut Objects
DeepMutSet <- function(lst){
  class(lst) <- c('DeepMutSet', class(lst))
  return(lst)
}

# Generic for writing deep mut file
write_deep_mut <- function(x, ...){
  UseMethod('write_deep_mut', x)
}

# Default Write deep mut method
write_deep_mut.default <- function(x){
  stop(str_c('Error: to write a dm file(s) "x" must be of class "DeepMut" or "DeepMutSet" but class(x) = ', class(x)))
}

# write 'DeepMut' classed objects to a consistent file type (termed dm file for now)
write_deep_mut.DeepMut <- function(x, outfile, create_dir=TRUE){
  keys <- names(x)
  keys <- keys[!keys == 'variant_data']
  
  # Check if dir exists and create if necessary
  target_dir <- str_split(outfile, '/')[[1]]
  target_dir <- str_c(target_dir[-length(target_dir)], collapse = '/')
  if (!dir.exists(target_dir)){
    if (create_dir){
      dir.create(target_dir, recursive = TRUE)
    } else {
      stop('Target directory does not exist')
    }
  }
  
  ## Write meta data
  write_lines(c('#deep_mut_file_version:2.1.0'), outfile)
  
  gene_keys <- c('gene_name', 'domain', 'species', 'alt_name')
  acc_keys <- keys[grepl('_id', keys)]
  study_keys <- c('authour', 'year', 'title', 'pmid', 'url', 'doi')
  seq_keys <- keys[grepl('_seq', keys)]
  misc_keys <- keys[!keys %in% c(gene_keys, acc_keys, study_keys, seq_keys, 'deep_mut_file_version')]
  
  ordered_keys <- c(gene_keys, acc_keys, study_keys, misc_keys, seq_keys)
  ordered_keys <- ordered_keys[ordered_keys %in% keys]
  
  # Write meta data
  for (k in ordered_keys){
    if (length(x[[k]]) > 1){
      # Write list/vectors
      write_lines(str_c('#*', k, ':', length(x[[k]])), outfile, append = TRUE)
      write_lines(str_c('#*', as.character(x[[k]])), outfile, append = TRUE)
      
    } else if (grepl('_seq', k)){
      # write long seq strings to multiple lines
      if (!is.na(x[[k]])){
        s <- str_split(x[[k]], '')[[1]]
        l <- ceiling(length(s)/80)
        
        split_s <- sapply(1:l, function(i){
          t <- s[((i - 1) * 80 + 1):(i*80)];
          return(str_c(t[!is.na(t)], collapse = ''))
        })
        
        # Header has number of lines to follow, all lines denoted '#+'
        write_lines(str_c('#+', k, ':', length(split_s)), outfile, append = TRUE)
        write_lines(str_c('#+', split_s), outfile, append = TRUE)
      } else {
        write_lines(str_c('#', k, ':NA', outfile, append = TRUE))
      }
    } else {
      # Write generic (e.g. short string)
      write_lines(str_c('#', k, ':', ifelse(is.na(x[[k]]), 'NA', x[[k]])), outfile, append = TRUE)
    }
    
  }
  
  ## Write variant table (header line (& final metadata line) denoted by '?')
  write_lines(str_c(c(str_c('?',colnames(x$variant_data)[1]),
                      colnames(x$variant_data)[-1]),
                   collapse = '\t'),
              outfile, append = TRUE)
  write_tsv(x$variant_data, outfile, append = TRUE, col_names = FALSE)
}

write_deep_mut.DeepMutSet <- function(x, root_dir, create_dir=TRUE){
  if (!dir.exists(root_dir)){
    if (create_dir){
      dir.create(root_dir, recursive = TRUE)
    } else {
      stop('Directory does not exist')
    }
  }
  for (set in names(x)){
    dm_path <- gsub(' ', '',
                    str_c(root_dir, '/', str_replace_na(x[[set]]$uniprot_id, replacement = ''), '_',
                          str_replace_na(x[[set]]$gene_name, replacement = ''), '.', set, '.dm'))
    write_deep_mut(x[[set]], dm_path)
  }
}

# Read deep mutagenesis data from a 'dm' file and return a DeepMut object
read_deep_mut <- function(filepath){
  tbl <- rename_all(read_tsv(filepath, comment = '#', col_names = TRUE), .funs = ~ gsub('\\?', '', .))
  dm <- DeepMut(tbl)
  
  fi <- file(filepath, 'r')
  while (TRUE){
    ln <- readLines(fi, n = 1)
    # Check for end of file
    if (length(ln) == 0){
      close(fi)
      stop('Reached end of file without encountering header line (marked "?")')
    } 
    
    first_char <- str_sub(ln, end=1)
    # Check line is a meta line
    if (first_char == '#'){
      ln <- str_sub(ln, 2)
      second_char <- str_sub(ln, end=1)
      
      pair <- str_split(ln, ':')[[1]]
      if (length(pair) > 2){
        # Reform any values that contain ':' chars
        pair <- c(pair[1], str_c(pair[-1], collapse = ':'))
      }
      
      if (second_char == '*'){
        # Read lists
        pair[1] <- str_sub(pair[1], 2)
        if (pair[2] == 'NA'){
          dm[[pair[1]]] <- NA
        } else {
          s <- readLines(fi, n = as.numeric(pair[2]))
          dm[[pair[1]]] <- gsub('\\#\\*', '', s)
        }
        
      } else if (second_char == '+'){
        # Read multiline strings
        pair[1] <- str_sub(pair[1], 2)
        if (pair[2] == 'NA'){
          dm[[pair[1]]] <- NA
        } else {
          s <- readLines(fi, n = as.numeric(pair[2]))
          dm[[pair[1]]] <- gsub('\\#\\+', '', str_c(s, collapse = ''))
        }
      } else {
        # Read normal entries
        if (pair[2] == 'NA'){
          dm[[pair[1]]] <- NA
        } else {
          dm[[pair[1]]] <- pair[2]
        }
      }
      
      # Convert numeric values
      if (all(grepl('^\\-?[0-9]*(\\.[0-9]*)?$', dm[[pair[1]]]))){
        dm[[pair[1]]] <- as.numeric(dm[[pair[1]]])
      }
    } else if (first_char == '?'){
      # Detect header line and stop parsing
      close(fi)
      break
    } else {
      # Otherwise file is not formatted properly
      close(fi)
      stop('Reached end of meta lines without encountering header line (marked "?")')
    }
  }
  return(dm)
}

# Get set of unique variants from deep mut data
get_variants <- function(x, ...){
  UseMethod('get_variants', x)
}

get_variants.default <- function(x){
  stop(str_c('Error: to fetch variants "x" must be of class "DeepMut" or "DeepMutSet" but class(x) = ', class(x)))
}

get_variants.DeepMut <- function(x){
  muts <- str_replace_all(x$variant_data$variants, 'p.', '')
  muts <- str_split(muts, ',')
  muts <- unique(unlist(muts))
  return(muts[!is.na(muts)])
}

get_variants.DeepMutSet <- function(x){
  muts <- lapply(x, get_variants)
  muts <- unique(unlist(muts))
  return(muts)
}

# Fetch meta data about deep mut data
get_meta <- function(x, var){
  if ('DeepMutSet' %in% class(x)){
    return(sapply(x, get_meta, var=var))
  } else if ('DeepMut' %in% class(x)){
    return(x[[var]])
  }
}

# Fetch deep mut dataset size
get_size <- function(x){
  if ('DeepMutSet' %in% class(x)){
    return(sapply(x, get_size))
  } else if ('DeepMut' %in% class(x)){
    return(dim(x$variant_data)[1])
  }
}

#### Analysis ####
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
  
  # Apply simple pseudocount of minimum non zero
  pseudo <- min(sel[sel>0], na.rm = TRUE)
  sel <- sel + pseudo
  bkg <- bkg + pseudo
  
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


