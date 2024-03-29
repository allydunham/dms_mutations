#!/usr/bin/env Rscript 
# Functions used to load and process data on variants from external tools (SIFT, FoldX, etc.)

#### General ####
import_dm_predictions_dataset <- function(dm_file, per_codon=FALSE){
  root_dir <- dirname(dm_file)
  
  # Read DeepMut File and add normalised scores
  dm <- read_deep_mut(dm_file)
  norm_factor <- max(abs(dm$variant_data$score), na.rm = TRUE)
  dm$variant_data <- mutate(dm$variant_data, norm_score = score / norm_factor)
  
  dataset_name <- make_dm_dataset_name(dm)
  message(dataset_name)
  
  # Import prediction data
  sift_file <- str_c(root_dir, '/', dm$gene_name, '.SIFTprediction')
  if (file.exists(sift_file)){
    sift <- read_sift_predictions(sift_file)
  } else {
    sift <- NA
  }
  
  foldx <- read_dm_foldx_predictions(dm$pdb_id, root_dir)
  
  env <- select_envision_data(root_dir, dm$gene_name, dm$uniprot_id)
  
  pph_file <- str_c(root_dir, '/pph_', dm$gene_name,'.predictions')
  if (file.exists(pph_file)){
    pph <- read_polyphen2_predictions(pph_file)
  } else {
    pph <- NA
  }
  
  evcoup_file <- str_c(root_dir, '/ev/mutate/ev_dataset_predicted.csv')
  if (file.exists(evcoup_file)){
    evcoup <- read_evcouplings_predictions(evcoup_file)
  } else {
    evcoup <- NA
  }
  
  naccess <- read_all_accesibility(dm$pdb_id, root_dir)
  
  secondary_structure <- read_secondary_structure(root_dir, dm)
  
  chem_env <- read_all_chemical_environment(dm$pdb_id, root_dir)
  
  phi_psi <- read_all_phi_psi(dm$pdb_id, root_dir)
  
  # If dataset is given per codon take mean per codon
  if (per_codon){
    dm <- adjust_dm_per_codon(dm, norm_factor)
    
    if (!identical(foldx, NA)){
      foldx <- lapply(foldx, adjust_foldx_per_codon)
    }
    if (!identical(evcoup, NA)){
      evcoup <- adjust_evcouplings_per_codon(evcoup)
    }
  }
  
  predictions <- list(dm=dm,
                      sift=sift,
                      envision=env,
                      polyphen2=pph,
                      foldx=foldx,
                      evcouplings=evcoup,
                      surface_accesibility=naccess,
                      secondary_structure=secondary_structure,
                      chem_env=chem_env,
                      backbone_angles=phi_psi,
                      norm_factor=norm_factor,
                      manual_threshold=unname(MANUAL_THRESHOLDS[dataset_name]))
  
  if (any(grepl('\\,', dm$variant_data$variants))){
    cls <- 'multi_variant'
  } else {
    cls <- 'single_variant'
  }
  class(predictions) <- c(cls, class(predictions))

  predictions <- combine_dms_predictions(predictions)
  return(predictions)
}

# Take mean of repeated variants enrichment score if there are repeats
adjust_dm_per_codon <- function(dm, norm_factor){
  dm$variant_data <- dm$variant_data %>%
    group_by(variants) %>%
    summarise(score = mean(score, na.rm=TRUE),
              raw_score = mean(raw_score, na.rm=TRUE)) %>%
    mutate(norm_score = score / norm_factor,
           position = as.integer(str_sub(variants, start = 4, end = -2))) %>%
    arrange(position) %>%
    select(-position)
  
  return(dm)
}

#### Combine variant data ####
combine_dms_predictions <- function(x, ...){
  UseMethod('combine_dms_predictions', x)
}

# If only single variant predictions are present
combine_dms_predictions.single_variant <- function(x){
  x$single_variants <- x$dm$variant_data %>%
    mutate(variants = str_replace_all(variants, 'p\\.', ''))
  
  if (!identical(x$evcoup, NA)){
    x$single_variants <- left_join(x$single_variants, select(x$evcoup, variants=mutant, evcoup_epistatic=prediction_epistatic,
                                                             evcoup_independent=prediction_independent), by='variants')
  }
  
  for (i in names(x$foldx)){
    re <- structure(c('sd', 'total_energy'), names=c(str_c('foldx_', i, '_sd'), str_c('foldx_', i, '_ddG')))
    x$single_variants %<>% left_join(., select(x$foldx[[i]], variants, !!re), by='variants')
  }
  
  if (!identical(x$polyphen2, NA)){
    x$single_variants <- left_join(x$single_variants, select(x$polyphen2, variants, pph2_prediction=prediction, pph2_class, pph2_prob,
                                                             pph2_FPR, pph2_TPR, pph2_FDR), by='variants')
  }
  if (!identical(x$env, NA)){
    x$single_variants <- left_join(x$single_variants, select(x$env, variants=Variant, envision_prediction=Envision_predictions,
                                                             log2_envision_prediction), by='variants')
  }
  if (!identical(x$sift, NA)){
    x$single_variants <- left_join(x$single_variants, select(x$sift, variants=variant, sift_prediction,
                                                             sift_score, sift_median), by='variants')
  }
  return(x)
}

# If any multi-variant predictions are present
combine_dms_predictions.multi_variant <- function(x){
  x$multi_variants <- x$dm$variant_data %>%
    mutate(variants = str_replace_all(variants, 'p\\.', ''),
           count = factor(sapply(variants, function(x){dim(str_split(x, ',', simplify = TRUE))[2]})))
  
  
  if (!identical(x$evcoup, NA)){
    x$multi_variants <- left_join(x$multi_variants, select(x$evcoup, variants=mutant, evcoup_epistatic=prediction_epistatic,
                                                       evcoup_independent=prediction_independent), by='variants')
  }
  
  for (i in names(x$foldx)){
    re <- structure(c('sd', 'total_energy'), names=c(str_c('foldx_', i, '_sd'), str_c('foldx_', i, '_ddG')))
    x$multi_variants %<>% left_join(., select(x$foldx[[i]], variants, !!re), by='variants')
  }
  
  # Create single variants table
  x$single_variants <- filter(x$multi_variants, count==1) %>%
    select(-count)
  
  # Take average of effect for each single variant that hasn't been directly measured
  max_vars <- dim(str_split(x$dm$variant_data$variants, ',', simplify = TRUE))[2]
  mean_single_variants <- x$dm$variant_data %>%
    mutate(variants = str_replace_all(variants, 'p\\.', '')) %>%
    separate('variants', str_c('variant', 1:max_vars, sep='_'), sep=',', extra='drop', fill='right') %>%
    gather(key = 'num', value = 'variants', contains('variant_')) %>%
    select(-num) %>%
    drop_na(variants) %>%
    group_by(variants) %>%
    summarise(sd=sd(score, na.rm = TRUE),
              score=mean(score, na.rm = TRUE), # currently just take mean, maybe use better metric?
              raw_score=mean(raw_score, na.rm=TRUE),
              norm_score=mean(norm_score, na.rm=TRUE),
              n=n()) %>%
    filter(!variants %in% x$single_variants$variants)
  
  x$single_variants <- bind_rows(x$single_variants, mean_single_variants)
  
  if (!identical(x$polyphen2, NA)){
    x$single_variants <- left_join(x$single_variants, select(x$polyphen2, variants, pph2_prediction=prediction, pph2_class, pph2_prob,
                                                             pph2_FPR, pph2_TPR, pph2_FDR), by='variants')
  }
  if (!identical(x$env, NA)){
    x$single_variants <- left_join(x$single_variants, select(x$env, variants=Variant, envision_prediction=Envision_predictions,
                                                             log2_envision_prediction), by='variants')
  }
  if (!identical(x$sift, NA)){
    x$single_variants <- left_join(x$single_variants, select(x$sift, variants=variant, sift_prediction,
                                                             sift_score, sift_median), by='variants')
  }
  return(x)
}

########

#### SIFT ####
read_sift_predictions <- function(filepath){
  return(read_tsv(filepath, col_names = c('variant', 'sift_prediction', 'sift_score', 'sift_median', 'num_seq', 'align_count')))
}

########

#### FoldX ####

# Read all foldx entries for a DM study, assuming my directory setup 
read_dm_foldx_predictions <- function(pdb_ids, root){
  if (length(pdb_ids) == 0){
    return(NA)
  }
  
  foldx <- list()
  for (pdb in str_split(pdb_ids, ':')){
    pdb_id <- pdb[1]
    pdb_offset <- ifelse(is.na(pdb[3]), 0, as.integer(pdb[3]))
    
    muts_file <- str_c(root, '/', pdb_id, '/individual_list_', pdb_id, '.txt')
    ddg_file <- str_c(root, '/', pdb_id, '/Average_', pdb_id, '_Repair.fxout')
    if (file.exists(muts_file) & file.exists(ddg_file)){
      foldx[[pdb_id]] <- read_foldx_ddg(ddg_file, muts_file, offset=pdb_offset)
    }
  }
  
  if (length(foldx) > 0){
    return(foldx)
  } else {
    return(NA)
  }
}

# read FoldX ddG predictions from FoldX output and a file giving the input variants 
read_foldx_ddg <- function(ddg_path, muts_path, offset=0){
  muts <- read_lines(muts_path) %>%
    str_replace(., ';', '') %>%
    lapply(., format_pdb_variants, pdb_offset=offset) %>%
    unlist()
  
  ddg <- read_tsv(ddg_path, skip = 8,
                  col_types = cols(.default=col_double(), Pdb=col_character())) %>%
    rename_all(~ (str_to_lower(str_replace_all(., ' ', '_')))) %>%
    mutate(pdb=muts) %>%
    rename(variants=pdb)
  
  return(ddg)
}

# Format variant IDs output by FoldX
format_pdb_variants <- function(x, pdb_offset=0){
  x <- str_split(x, ',')[[1]]
  str_sub(x, 2, 2) <- ''
  str_sub(x, 2, -2) <- as.character(as.integer(str_sub(x, 2, -2)) + pdb_offset)
  return(str_c(x, collapse = ','))
}

# Take means if variants given multiple times (e.g. per codon dms)
adjust_foldx_per_codon <- function(x){
  return(
    group_by(x, variants) %>%
      summarise_all(.funs = mean) %>%
      mutate(position = as.integer(str_sub(variants, start = 2, end = -2))) %>%
      arrange(position) %>%
      select(-position)
  )
}

########

#### Envision ####
ENVISION_COLS <- cols(.default = col_double(), id2 = col_character(), AA1 = col_character(),
                      AA2 = col_character(), position = col_integer(), Uniprot = col_character(),
                      WT_Mut = col_character(), Variant = col_character(), AA1_polarity = col_character(),
                      AA2_polarity = col_character(), SecondaryStructure=col_character(),
                      phi_psi_angles=col_character())

ENVISION_ROOT <- 'data/envision_data'

# Choose appropriate dataset to load, based on presence of batch and manually selected results (prioritising manual)
# Assumes my setup (file {uniprot_id}_envisionData.csv in envision_root or {gene_name}_envision_vars.csv in the batch root dir)
select_envision_data <- function(batch_dir, gene_name, uniprot_id){
  env_path_manual <- str_c(ENVISION_ROOT, '/', uniprot_id, '_envisionData.csv')
  env_path_batch <- str_c(batch_dir, '/', gene_name, '_envision_vars.csv')
  
  env <- NA
  if (file.exists(env_path_manual)){
    env <- read_envision_predictions(env_path_manual)
  } else if (file.exists(env_path_batch)){
    env <- read_envision_predictions(env_path_batch)
  }
  return(env)
}

# Import Envision predictions
read_envision_predictions <- function(filepath, variant_filter=NULL){
  env <- read_csv(filepath, col_types = ENVISION_COLS)

  if (!is.null(variant_filter)){
    env <- filter(env, Variant %in% variant_filter)
  }
  
  env <- mutate(env, log2_envision_prediction = log2(Envision_predictions))
  return(env)
}

########

#### Polyphen2 ####
PPH_COLS = cols(o_pos=col_integer(),pos=col_integer(),pph2_prob=col_double(),pph2_FPR=col_double(),
                pph2_TPR=col_double(),pph2_FDR=col_double(),dScore=col_double(),Score1=col_double(),Score2=col_double(),
                MSAv=col_integer(),Nobs=col_integer(),Nstruct=col_integer(),Nfilt=col_integer(),PDB_pos=col_integer(),
                ident=col_double(),length=col_integer(),NormASA=col_double(),dVol=col_integer(),dProp=col_double(),
                `B-fact`=col_double(),`H-bonds`=col_double(),AveNHet=col_double(),MinDHet=col_double(),AveNInt=col_double(),
                MinDInt=col_double(),AveNSit=col_double(),MinDSit=col_double(),Transv=col_integer(),CodPos=col_integer(),
                CpG=col_integer(),MinDJxn=col_integer(),IdPmax=col_double(),IdPSNP=col_double(),
                IdQmin=col_double(), .default = col_character())

read_polyphen2_predictions <- function(filepath){
  pph <- read_tsv(filepath, col_types = PPH_COLS, na = '') %>%
    mutate(variants=str_c(aa1, pos, aa2))
  names(pph) <- str_replace(names(pph), '#', '')
  return(pph)
}

########

#### EVCouplings ####
read_evcouplings_predictions <- function(filepath){
  return(read_csv(filepath, col_types = cols(mutant = col_character(), .default = col_double())))
}

adjust_evcouplings_per_codon <- function(x){
  return(
    evcoup <- group_by(x, mutant) %>%
      summarise_all(.funs = mean) %>%
      mutate(position = as.integer(str_sub(mutant, start = 2, end = -2))) %>%
      arrange(position) %>%
      select(-position)
  )
}

########

#### naccess ####
read_all_accesibility <- function(pdb_ids, root){
  if (length(pdb_ids) == 0){
    return(NA)
  }
  
  naccess <- list()
  for (pdb in str_split(pdb_ids, ':')){
    pdb_id <- pdb[1]
    pdb_chain <- pdb[2]
    
    filepath <- str_c(root, '/', pdb_id, '/', pdb_id, '_Repair.rsa')
    if (file.exists(filepath)){
      naccess[[pdb_id]] <- read_naccess_rsa(filepath, pdb_chain=pdb_chain)
    }
  }
  
  if (length(naccess) == 0){
    return(NA)
  }
  
  # Combine accesiblity per residues into final table (take average where multiple structures cover residue)
  combined_accesibility <- sapply(names(naccess), function(x){naccess[[x]]$residues}, simplify = FALSE) %>%
    bind_rows() %>%
    select(-res3, -chain) %>%
    group_by(res1, pos) %>%
    summarise_all(.funs = list(~ mean(., na.rm=TRUE))) %>%
    arrange(pos) %>%
    ungroup()
    
  naccess$combined <- combined_accesibility
  return(naccess)
}
  
# Import a given naccess per residue output. chain filters which chain to keep
read_naccess_rsa <- function(filepath, pdb_chain=NULL){
  fi <- read_lines(filepath)
  fi <- grep('^REM', fi, invert = TRUE, value = TRUE)
  
  # Select table lines and read
  tbl_str <- str_replace(grep('^RES', fi, value = TRUE),'RES ', '')
  str_sub(tbl_str, 6, 5) <- ' ' # hack inserting space between chain and position in 4digit positions, which naccess doesn't do 
  acc <- read_table2(tbl_str, col_names = c('res3', 'chain', 'pos', 'all_atom_abs', 'all_atom_rel',
                                            'side_chain_abs', 'side_chain_rel', 'backbone_abs', 'backbone_rel',
                                            'non_polar_abs', 'non_polar_rel', 'polar_abs', 'polar_rel')) %>%
    mutate(res3 = str_to_title(res3)) %>%
    add_column(res1 = structure(names(Biostrings::AMINO_ACID_CODE), names = Biostrings::AMINO_ACID_CODE)[.$res3], .after = 'res3')
  
  if (!is.null(pdb_chain)){
    acc <- filter(acc, chain == pdb_chain)
  }
    
  # Parse summary lines
  chain <- str_split(grep('^CHAIN', fi, value = TRUE), '\\s+', simplify = TRUE)[,-c(1,2)]
  if (is.null(dim(chain))){
    chain <- set_names(chain, c('chain', 'all_atoms', 'side_chain', 'backbone', 'non_polar', 'polar'))
  } else {
    chain <- set_colnames(chain, c('chain', 'all_atoms', 'side_chain', 'backbone', 'non_polar', 'polar')) %>%
      as_tibble()
  }
    
  
  total <- str_split(grep('^TOTAL', fi, value = TRUE), '\\s+', simplify = TRUE)[,-1] %>%
    set_names(c('all_atoms', 'side_chain', 'backbone', 'non_polar', 'polar'))

  return(list(residues = acc, chains = chain, total = total))
}


########

#### Secondary Structure ####
read_secondary_structure <- function(root, dm){
  seq_preds <- read_ss8(str_c(root, '/', dm$gene_name, '.fa.ss8'))
  
  if (identical(NA, dm$pdb_id)){
    return(seq_preds)
  }
  
  struc_preds <- sapply(str_split(dm$pdb_id, ':', simplify = TRUE)[,1],
                        function(x){read_ss8(str_c(root, '/', x, '/', x, '.ss8'))},
                        simplify = FALSE) %>%
    bind_rows(.id = 'pdb') %>%
    group_by(pos, aa) %>%
    summarise_at(.vars = vars(-pdb), .funs = first)
  
  overall_preds <- seq_preds %>%
    filter(!pos %in% struc_preds$pos) %>%
    bind_rows(., struc_preds) %>%
    arrange(pos)
  
  return(overall_preds)
}

read_ss8 <- function(path){
  return(
    read_tsv(path) %>%
      rename(pos = `#`) %>%
      rename_all(.funs = str_to_lower)
  )
}

########

#### Chemical environment ####
read_all_chemical_environment <- function(pdb_ids, root_dir, ext='.chem_env'){
  if (identical(pdb_ids, NA)){
    return(NA)
  }
  
  pdb_ids <- str_split(pdb_ids, ':', simplify = TRUE)[,1]
  tbls <- sapply(pdb_ids, function(x){read_chemical_environment(str_c(root_dir, '/', x, '/', x, ext))}, simplify = FALSE)
  
  if (length(tbls) > 1){
    # Combining makes several assumptions to make sense:
    # - all residues listed are from the same protein, grouped in one way, with the same position numbering
    combine_wide <- Reduce(function(x, y){full_join(tbls[[x]], tbls[[y]], by = c('position', 'aa'), suffix=str_c('_', c(x,y)))}, names(tbls)) %>%
      select(position, aa, everything()) %>%
      mutate(num_pdbs = rowSums(select(., starts_with('pdb_id_')) %>% mutate_all(~ !is.na(.)))) %>%
      select(-starts_with('pdb_id_'))
  } else {
    combine_wide <- tbls[[1]] %>%
      select(position, aa, everything()) %>%
      rename_at(vars(-position, -aa), ~ str_c(., '_', names(tbls)[1]))
  }
  
  combine_long <- bind_rows(tbls) %>%
    arrange(position)

  return(list(combine_long=combine_long, combine_wide=combine_wide, split=tbls))
}

read_chemical_environment <- function(path){
  read_tsv(path) %>%
    mutate_at(.vars = vars(-pdb_id, -chain, -group, -position, -aa), .funs = ~ str_split(., ',')) %>%
    mutate_at(.vars = vars(-pdb_id, -chain, -group, -position, -aa), .funs = ~ lapply(., as.integer)) %>%
    return()
}

########

#### Backbone angles ####
read_all_phi_psi <- function(pdb_ids, root_dir, ext='.bb_angles'){
  if (identical(pdb_ids, NA)){
    return(NA)
  }
  
  pdb_ids <- str_split(pdb_ids, ':', simplify = TRUE)[,1]
  sapply(pdb_ids, function(x){read_tsv(str_c(root_dir, '/', x, '/', x, ext))}, simplify = FALSE) %>%
    bind_rows() %>%
    return()
}
########
