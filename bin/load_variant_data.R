#!/usr/bin/env Rscript 
# Script to load processed deep mutagenesis data for analysis

source('bin/libraries.R')
source('bin/dm_functions.R')
source('bin/variant_functions.R')

#### Load Data ####
deep_variant_data <- list()

deep_datasets <- dir('data/standardised')

for (dataset in deep_datasets){
  print(dataset)
  root <- str_c('data/standardised/', dataset)
  dm_files <- grep('*.dm', dir(root), value = TRUE)
  
  for (dm_file in dm_files){
    dataset_name <- str_c(c(str_split(dataset, '_', simplify = TRUE)[,c(1,2)], str_split(dm_file, '_', simplify = TRUE)[,2]), collapse='_') %>%
      str_replace(., '\\.dm', '') %>% 
      str_replace(., '\\.', '_')
    
    # Read DeepMut File
    dm <- read_deep_mut(str_c(root, '/', dm_file[1]))
    
    # Read SIFT Scores
    sift_file <- str_c(dm$gene_name, '.SIFTprediction')
    if (sift_file %in% dir(root)){
      sift <- read_tsv(str_c(root, '/', sift_file), col_names = c('variant', 'sift_prediction', 'sift_score', 'sift_median', 'num_seq', 'align_count'))
    } else {
      sift <- NA
    }
    
    # Read FoldX Scores
    foldx <- list()
    if (length(dm$pdb_id) > 0){
      pdb_ids <- dm$pdb_id
      for (pdb in str_split(pdb_ids, ':')){
        pdb_id <- pdb[1]
        pdb_offset <- ifelse(is.na(pdb[3]), 0, as.integer(pdb[3]))
        fx_file <- str_c(root, '/', pdb_id, '/individual_list_', pdb_id, '.txt')
        if (file.exists(fx_file)){
          muts <- read_lines(fx_file) %>%
            str_replace(., ';', '') %>%
            lapply(., format_pdb_variants, pdb_offset=pdb_offset) %>%
            unlist()
          
          ddg <- read_tsv(str_c(root, '/', pdb_id, '/Average_', pdb_id, '_Repair.fxout'), skip = 8,
                          col_types = cols(.default=col_double(), Pdb=col_character())) %>%
            rename_all(funs(str_to_lower(str_replace_all(., ' ', '_')))) %>%
            mutate(pdb=muts) %>%
            rename(variants=pdb)
          
          foldx[[pdb_id]] <- ddg        
        }
      }
    }
    if (length(foldx) == 0){
      foldx <- NA
    }
    
    # Read Envision Scores
    env_path_batch <- str_c('data/envision_data/', dm$uniprot_id, '_envisionData.csv')
    env_path_specific <- str_c(root, '/', dm$gene_name, '_envision_vars.csv')
    if (file.exists(env_path_batch)){
      env <- read_csv(env_path_batch, col_types = cols(.default = col_double(), id2 = col_character(), AA1 = col_character(),
                                                       AA2 = col_character(), position = col_integer(), Uniprot = col_character(),
                                                       WT_Mut = col_character(), Variant = col_character(), AA1_polarity = col_character(),
                                                       AA2_polarity = col_character(), SecondaryStructure=col_character(),
                                                       phi_psi_angles=col_character())) %>%
        filter(Variant %in% get_variants(dm)) %>%
        mutate(log2_envision_prediction = log2(Envision_predictions))
    } else if (file.exists(env_path_specific)){
      env <- read_csv(env_path_specific, col_types = cols(.default = col_double(), id2 = col_character(), AA1 = col_character(),
                                                          AA2 = col_character(), position = col_integer(), Uniprot = col_character(),
                                                          WT_Mut = col_character(), Variant = col_character(), AA1_polarity = col_character(),
                                                          AA2_polarity = col_character(), SecondaryStructure=col_character(),
                                                          phi_psi_angles=col_character())) %>%
        mutate(log2_envision_prediction = log2(Envision_predictions))
    } else {
      env <- NA
    }
  
    # Read EVCouplings Scores
    evcoup_path <- str_c(root, '/ev/mutate/ev_dataset_predicted.csv')
    if ('ev' %in% dir(root)){
      evcoup <- read_csv(evcoup_path)
    } else {
      evcoup <- NA
    }
    
    # Read PolyPhen2 Scores
    pph_file <- str_c('pph_', dm$gene_name,'.predictions')
    if (pph_file %in% dir(root)){
      pph <- read_tsv(str_c(root, '/', pph_file), col_types = cols(effect=col_character())) %>%
        mutate(variants=str_c(aa1, pos, aa2))
      names(pph) <- str_replace(names(pph), '#', '')
    } else {
      pph <- NA
    }
    
    # Create Combined Score Table
    if (any(grepl('\\,', dm$variant_data$variants))){
      # Datasets with multiple mutations
      cls <- 'multi_variant'
      
      multi_variants <- dm$variant_data %>%
        mutate(variants=str_replace_all(variants, 'p\\.', ''),
               count = factor(sapply(variants, function(x){dim(str_split(x, ',', simplify = TRUE))[2]})))
        
      
      if (!all(is.na(evcoup))){
        multi_variants <- left_join(multi_variants, select(evcoup, variants=mutant, evcoup_epistatic=prediction_epistatic,
                                                           evcoup_independent=prediction_independent), by='variants')
      }
      
      for (i in names(foldx)){
        re <- structure(c('sd', 'total_energy'), names=c(str_c('foldx_', i, '_sd'), str_c('foldx_', i, '_ddG')))
        multi_variants %<>% left_join(., select(foldx[[i]], variants, !!re), by='variants')
      }
      
      # Take average of effect for each single variant
      max_vars <- dim(str_split(dm$variant_data$variants, ',', simplify = TRUE))[2]
      single_variants <- dm$variant_data %>%
        mutate(variants=str_replace_all(variants, 'p\\.', '')) %>%
        separate('variants', str_c('variant', 1:max_vars, sep='_'), sep=',', extra='drop', fill='right') %>%
        gather(key = 'num', value = 'variants', contains('variant_')) %>%
        select(-num) %>%
        group_by(variants) %>%
        summarise(sd=sd(score, na.rm = TRUE),
                  score=mean(score, na.rm = TRUE),
                  raw_score=mean(raw_score, na.rm=TRUE),
                  n=n())# currently just take mean, maybe use better metric?
      
      if (!all(is.na(pph))){
        single_variants <- left_join(single_variants, select(pph, variants, pph2_prediction=prediction, pph2_class, pph2_prob,
                                                             pph2_FPR, pph2_TPR, pph2_FDR), by='variants')
      }
      
      if (!all(is.na(env))){
        single_variants <- left_join(single_variants, select(env, variants=Variant, envision_prediction=Envision_predictions,
                                                             log2_envision_prediction), by='variants')
      }
      
      if (!all(is.na(sift))){
        single_variants <- left_join(single_variants, select(sift, variants=variant, sift_prediction,
                                                             sift_score, sift_median), by='variants')
      }
        
    } else {
      # Datasets with single variants only
      cls <- 'single_variant'
      multi_variants <- NULL
      
      single_variants <- dm$variant_data %>%
        mutate(variants=str_replace_all(variants, 'p\\.', ''))
        
      if (!all(is.na(pph))){
        single_variants <- left_join(single_variants, select(pph, variants, pph2_prediction=prediction, pph2_class,
                                                             pph2_prob, pph2_FPR, pph2_TPR, pph2_FDR), by='variants')
      }
      if (!all(is.na(env))){
        single_variants <- left_join(single_variants, select(env, variants=Variant, envision_prediction=Envision_predictions,
                                                             log2_envision_prediction), by='variants')
      }
      if (!all(is.na(sift))){
        single_variants <- left_join(single_variants, select(sift, variants=variant, sift_prediction,
                                                             sift_score, sift_median), by='variants')
      }
      if (!all(is.na(evcoup))){
        single_variants <- left_join(single_variants, select(evcoup, variants=mutant, evcoup_epistatic=prediction_epistatic,
                                                             evcoup_independent=prediction_independent), by='variants')
      }
      
      for (i in names(foldx)){
        re <- structure(c('sd', 'total_energy'), names=c(str_c('foldx_', i, '_sd'), str_c('foldx_', i, '_ddG')))
        single_variants %<>% left_join(., select(foldx[[i]], variants, !!re), by='variants')
      }
    }
    
    # Generate Data Object
    deep_variant_data[[dataset]] <- list(dm=dm,
                                         sift=sift,
                                         envision=env,
                                         polyphen2=pph,
                                         foldx=foldx,
                                         evcouplings=evcoup,
                                         single_variants=single_variants,
                                         multi_variants=multi_variants)
    
    class(deep_variant_data[[dataset]]) <- c(cls, class(deep_variant_data[[dataset]]))
  }
}

# Save generated dataset
write_rds(deep_variant_data, 'data/variant_data.RDS')

