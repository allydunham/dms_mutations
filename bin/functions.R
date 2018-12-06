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
# ref_seq = reference AA seq string
# authour = paper authours
# year = paper year (minimal info to locate, best to include more (such as url) in misc)
# transform = method used to transform score
# misc = named list of other meta data to attach (e.g. other protein accessions, other study data.
# Various special misc values are expected and are written in a logical positon - alt_name, doi, pubmed_id, url, title
# anything else with _id is treated as a gene id
DeepMut <- function(variant_data, gene_name=NA, domain=NA, species=NA, ref_seq=NA, transform='None',
                    uniprot_id=NA, authour=NA, year=NA, misc=NULL){
  # (minimal) Error Checking
  if (!all(c('variants', 'score', 'raw_score') %in% colnames(variant_data))) {
    stop('variant_data does not contain the required columns (dna_variants, protein_variants, score, raw_score)')
  }
  
  # Construct List of required fields
  deep_mut <- list(variant_data=variant_data,
                   gene_name=gene_name,
                   domain=domain,
                   species=species,
                   uniprot_id=uniprot_id,
                   ref_seq=ref_seq,
                   authour=authour,
                   year=year,
                   transform=transform)
  
  if (!is_null(misc)){
    deep_mut <- c(deep_mut, misc)
  }
  
  class(deep_mut) <- 'DeepMut'
  return(deep_mut)
}

# write 'DeepMut' classed objects to a consistent file type (termed dm file for now)
write_deep_mut <- function(x, outfile){
  if (!'DeepMut' %in% class(x)){
    stop('x must be an object of class DeepMut (see DeepMut() function)')
  }
  keys <- names(x)
  
  ## Write meta data
  write_lines(c('#deep_mut_file_version:1.1'), outfile)
  
  gene_keys <- c('gene_name', 'domain', 'species', 'alt_name')
  acc_keys <- keys[grepl('_id', keys) & !keys == 'pubmed_id']
  study_keys <- c('authour', 'year', 'title', 'pubmed_id', 'url', 'doi')
  misc_keys <- keys[!keys %in% c('gene_name', 'domain', 'species', 'alt_name', 'authour', 'year', 'title',
                                 'pubmed_id', 'url', 'doi', 'transform', 'ref_seq', 'variant_data') &
                      !grepl('_id', keys)]
  
  ordered_keys <- c(gene_keys, acc_keys, study_keys, 'transform', misc_keys)
  ordered_keys <- ordered_keys[ordered_keys %in% keys]
  
  for (k in ordered_keys){
    write_lines(str_c('#', k, ':', ifelse(is.na(x[[k]]), 'NA', x[[k]])), outfile, append = TRUE)
  }
  
  ## Write ref sequence
  if (!(is.na(x$ref_seq) & is.null(x$ref_seq))){
    seq <- str_split(x$ref_seq, '')[[1]]
    l <- ceiling(length(seq)/80)
    
    # Following sequence lines denoted as '#+'
    split_seq <- sapply(1:l, function(i){
      t <- seq[((i - 1) * 80 + 1):(i*80)];
      return(str_c(t[!is.na(t)], collapse = ''))
    })
    
    # Header with number of seq lines to follow
    write_lines(str_c('#ref_seq:', length(split_seq)), outfile, append = TRUE)
    # Seq lines, led by #+
    write_lines(str_c('#+', split_seq), outfile, append = TRUE)
  } else {
    write_lines('#ref_seq:NA', outfile, append = TRUE)
  }
  
  ## Write variant table (header line (& final metadata line) denoted by '?')
  write_lines(str_c(c(str_c('?',colnames(x$variant_data)[1]),
                      colnames(x$variant_data)[-1]),
                   collapse = '\t'),
              outfile, append = TRUE)
  write_tsv(x$variant_data, outfile, append = TRUE, col_names = FALSE)
}

# Read deep mutagenesis data from a 'dm' file and return a DeepMut object
read_deep_mut <- function(filepath){
  tbl <- rename_all(read_tsv(filepath, comment = '#', col_names = TRUE), funs(gsub('\\?', '', .)))
  dm <- DeepMut(tbl)
  
  fi <- file(filepath, 'r')
  ln <- readLines(fi, n = 1)
  while (TRUE){
    # Check for end of file
    if (length(ln) == 0){
      close(fi)
      stop('Reached end of file without encountering header line (marked "?")')
    } 
    
    first_char <- str_sub(ln, end=1)
    # Check line is a meta line
    if (first_char == '#'){
      ln <- str_sub(ln, 2)
      
      # Check for seq line before leader
      if (str_sub(ln, end=1) == '+'){
        close(fi)
        stop('Error: Reached sequence lines before finding a "#ref_seq" line')
      }
      
      # process meta pair
      pair <- str_split(ln, ':')[[1]]
      if (length(pair) > 2){
        # Reform any values that contain ':' chars
        pair <- c(pair[1], str_c(pair[-1], collapse = ':'))
      }
      
      if (pair[1] == 'ref_seq'){
        # Read in seq lines
        seq <- readLines(fi, n = pair[2])
        dm$ref_seq <- gsub('\\#\\+', '', str_c(seq, collapse = ''))
        
      } else if (!pair[1] == 'deep_mut_file_version'){
        # Read normal meta pairs (ignore file version), cheking if they look like a number
        if (grepl('^\\-?[0-9]*(\\.[0-9]*)?$', pair[2])){
          dm[[pair[1]]] <- as.numeric(pair[2])
        } else {
          dm[[pair[1]]] <- pair[2]
        }
      }
      ln <- readLines(fi, n = 1)
      
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

# tt <- deep_mut_data$hietpas_2011_hsp90
# tt$variants <- paste0('X', tt$position, tt$alt_aa)
# tt$score <- tt$selection_coefficient
# tt$raw_score <- tt$selection_coefficient
# tt <- tt[,c('variants', 'score', 'raw_score')]
# dm <- DeepMut(tt, gene_name = 'hsc82', alt_name = 'hsp90', accessions = c(uniprot_id='P02829', ensembl_id='AJS72977'),
#               species = 'Saccharomyces cerevisiae', study = c(year='2011', authour='Hietpas et al.'))

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


