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
# domain = domain string
# accessions = named chr vector of accessions
# study = named chr vector of paper metadata
# species = species string
# ref_seq = reference AA seq string
# len = gene length (overwritten by ref_seq len if given)
# transform = method used to transform score
DeepMut <- function(variant_data, gene_name=NULL, alt_name = NULL, domain=NULL, accessions=NULL,
                    study=NULL, species=NULL, ref_seq=NULL, len=NULL, transform='None'){
  
  # (minimal) Error Checking
  if (!all(c('variants', 'score', 'raw_score') %in% colnames(variant_data))) {
    stop('variant_data does not contain the required columns (dna_variants, protein_variants, score, raw_score)')
  }
  if (!is_null(ref_seq)){
    len <- nchar(ref_seq)
  }

  # Construct List
  deep_mut <- list(varaiant_data = variant_data,
                   gene_name = gene_name,
                   alt_name = alt_name,
                   domain = domain,
                   accessions = accessions,
                   study = study,
                   species = species,
                   ref_seq = ref_seq,
                   len = len)
  
  class(deep_mut) <- 'DeepMut'
  return(deep_mut)
}

# write 'DeepMut' classed objects to a consistent file type
write_deep_mut <- function(x, outfile){
  if (!'DeepMut' %in% class(x)){
    stop('x must be an object of class DeepMut (see DeepMut() function)')
  }
  
  # Write file meta data
  write_lines(c('#deep_mut_file_version:1.0'), outfile)
  
  # Write gene data
  write_lines(str_c('#gene_name:', x$gene_name), outfile, append = TRUE)
  
  if (!is_null(x$alt_name)){
    write_lines(str_c('#alt_name:', x$alt_name), outfile, append = TRUE)
  }
  
  write_lines(str_c('#', names(x$accessions), ':', x$accessions), outfile, append = TRUE)
  write_lines(str_c('#domain:', x$domain), outfile, append = TRUE)
  write_lines(str_c('#species:', x$species), outfile, append = TRUE)
  write_lines(str_c('#amino_acid_length:', x$len), outfile, append = TRUE)
  
  # Write study information
  write_lines(str_c('#study_', names(x$study), ':', x$study), outfile, append = TRUE)
  
  # Write ref sequence
  # Header line
  write_lines('#ref_seq:', outfile, append = TRUE)
  seq <- str_split(x$ref_seq, '')[[1]]
  l <- ceiling(length(seq)/80)

  # Following sequence lines denoted as '#+'
  split_seq <- sapply(1:l, function(i){
    t <- seq[((i - 1) * 80 + 1):(i*80)];
    return(str_c(t[!is.na(t)], collapse = ''))
    })
  
  write_lines(str_c('#+', split_seq), outfile, append = TRUE)
  
  # Write variant table (header line (& final metadata line) denoted by '#?')
  write_lines(str_c(c(str_c('#?',colnames(x$varaiant_data)[1]),
                      colnames(x$varaiant_data)[-1]),
                   collapse = '\t'),
              outfile, append = TRUE)
  write_tsv(x$varaiant_data, outfile, append = TRUE, col_names = FALSE)
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


