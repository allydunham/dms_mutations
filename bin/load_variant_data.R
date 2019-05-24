#!/usr/bin/env Rscript 
# Script to load processed deep mutagenesis data for analysis

source('bin/config.R')

#### Load Data ####
deep_variant_data <- list()

deep_datasets <- dir('data/standardised')

per_codon_datasets <- c('hietpas_2011_hsp90', 'weile_2017_sumo1', 'weile_2017_ube2i', 'findlay_2018_brca1', 'firnberg_2014_tem1')

pph_col_classes = cols(o_pos=col_integer(),pos=col_integer(),pph2_prob=col_double(),pph2_FPR=col_double(),
                       pph2_TPR=col_double(),pph2_FDR=col_double(),dScore=col_double(),Score1=col_double(),Score2=col_double(),
                       MSAv=col_integer(),Nobs=col_integer(),Nstruct=col_integer(),Nfilt=col_integer(),PDB_pos=col_integer(),
                       ident=col_double(),length=col_integer(),NormASA=col_double(),dVol=col_integer(),dProp=col_double(),
                       `B-fact`=col_double(),`H-bonds`=col_double(),AveNHet=col_double(),MinDHet=col_double(),AveNInt=col_double(),
                       MinDInt=col_double(),AveNSit=col_double(),MinDSit=col_double(),Transv=col_integer(),CodPos=col_integer(),
                       CpG=col_integer(),MinDJxn=col_integer(),IdPmax=col_double(),IdPSNP=col_double(),
                       IdQmin=col_double(), .default = col_character())

# TODO Refactor into smaller functions
for (dataset in deep_datasets){
  print(dataset)
  root <- str_c('data/standardised/', dataset)
  dm_files <- grep('*.dm', dir(root), value = TRUE)
  
  for (dm_file in dm_files){
    # Read DeepMut File
    dm <- read_deep_mut(str_c(root, '/', dm_file))
    norm_factor <- max(abs(dm$variant_data$score), na.rm = TRUE)
    dm$variant_data <- mutate(dm$variant_data, norm_score = score / norm_factor)
    
    #Generate dataset name
    dm_batch <- tryCatch(str_split(str_remove(dm_file, '\\.dm'), '\\.', simplify = TRUE)[1,2],
                                   error=function(cond){return(NULL)})
    dataset_name <- str_c(str_split(dm$authour, ' ', simplify = TRUE)[1,1], '_',
                          dm$year, '_',
                          dm$gene_name, ifelse(is.null(dm_batch), '', '.'),
                          dm_batch) %>%
      str_to_lower() %>%
      str_replace_all(., '-', '_')
    
    
    
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
            rename_all(~ (str_to_lower(str_replace_all(., ' ', '_')))) %>%
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
      evcoup <- read_csv(evcoup_path, col_types = cols(mutant = col_character(), .default = col_double()))
    } else {
      evcoup <- NA
    }
    
    # Read PolyPhen2 Scores
    pph_file <- str_c('pph_', dm$gene_name,'.predictions')
    if (pph_file %in% dir(root)){
      pph <- read_tsv(str_c(root, '/', pph_file), col_types = pph_col_classes, na = '') %>%
        mutate(variants=str_c(aa1, pos, aa2))
      names(pph) <- str_replace(names(pph), '#', '')
    } else {
      pph <- NA
    }
    
    # Adjust scores that are given per codon not per AA
    if (dataset %in% per_codon_datasets){
      dm$variant_data <- dm$variant_data %>%
        group_by(variants) %>%
        summarise(score = mean(score, na.rm=TRUE),
                  raw_score = mean(raw_score, na.rm=TRUE)) %>%
        mutate(norm_score = score / norm_factor,
               position = as.integer(str_sub(variants, start = 4, end = -2))) %>%
        arrange(position) %>%
        select(-position)
      
      if (!identical(foldx, NA)){
        foldx <- lapply(foldx, function(x){
          group_by(x, variants) %>%
            summarise_all(.funs = mean) %>%
            mutate(position = as.integer(str_sub(variants, start = 2, end = -2))) %>%
            arrange(position) %>%
            select(-position)
          })
      }
      if (!identical(evcoup, NA)){
        evcoup <- group_by(evcoup, mutant) %>%
          summarise_all(.funs = mean) %>%
          mutate(position = as.integer(str_sub(mutant, start = 2, end = -2))) %>%
          arrange(position) %>%
          select(-position)
      }
    }
    
    # Create Combined Score Table
    if (any(grepl('\\,', dm$variant_data$variants))){
      # Datasets with multiple mutations
      cls <- 'multi_variant'
      
      multi_variants <- dm$variant_data %>%
        mutate(variants = str_replace_all(variants, 'p\\.', ''),
               count = factor(sapply(variants, function(x){dim(str_split(x, ',', simplify = TRUE))[2]})))
        
      
      if (!identical(evcoup, NA)){
        multi_variants <- left_join(multi_variants, select(evcoup, variants=mutant, evcoup_epistatic=prediction_epistatic,
                                                           evcoup_independent=prediction_independent), by='variants')
      }
      
      for (i in names(foldx)){
        re <- structure(c('sd', 'total_energy'), names=c(str_c('foldx_', i, '_sd'), str_c('foldx_', i, '_ddG')))
        multi_variants %<>% left_join(., select(foldx[[i]], variants, !!re), by='variants')
      }
      
      # Create single variants table
      single_variants <- filter(multi_variants, count==1) %>%
        select(-count)
      
      # Take average of effect for each single variant that hasn't been directly measured
      max_vars <- dim(str_split(dm$variant_data$variants, ',', simplify = TRUE))[2]
      mean_single_variants <- dm$variant_data %>%
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
        filter(!variants %in% single_variants$variants)
      
      single_variants <- bind_rows(single_variants, mean_single_variants)
        
    } else {
      # Datasets with single variants only
      cls <- 'single_variant'
      multi_variants <- NULL
      
      single_variants <- dm$variant_data %>%
        mutate(variants = str_replace_all(variants, 'p\\.', ''))
        
      if (!identical(evcoup, NA)){
        single_variants <- left_join(single_variants, select(evcoup, variants=mutant, evcoup_epistatic=prediction_epistatic,
                                                             evcoup_independent=prediction_independent), by='variants')
      }
      
      for (i in names(foldx)){
        re <- structure(c('sd', 'total_energy'), names=c(str_c('foldx_', i, '_sd'), str_c('foldx_', i, '_ddG')))
        single_variants %<>% left_join(., select(foldx[[i]], variants, !!re), by='variants')
      }
    }
    
    # Add common values to single variants
    if (!identical(pph, NA)){
      single_variants <- left_join(single_variants, select(pph, variants, pph2_prediction=prediction, pph2_class, pph2_prob,
                                                           pph2_FPR, pph2_TPR, pph2_FDR), by='variants')
    }
    if (!identical(env, NA)){
      single_variants <- left_join(single_variants, select(env, variants=Variant, envision_prediction=Envision_predictions,
                                                           log2_envision_prediction), by='variants')
    }
    if (!identical(sift, NA)){
      single_variants <- left_join(single_variants, select(sift, variants=variant, sift_prediction,
                                                           sift_score, sift_median), by='variants')
    }
    
    # Generate Data Object
    deep_variant_data[[dataset_name]] <- list(dm=dm,
                                         sift=sift,
                                         envision=env,
                                         polyphen2=pph,
                                         foldx=foldx,
                                         evcouplings=evcoup,
                                         single_variants=single_variants,
                                         multi_variants=multi_variants,
                                         norm_factor=norm_factor,
                                         manual_threshold=unname(MANUAL_THRESHOLDS[dataset_name]))
    
    class(deep_variant_data[[dataset_name]]) <- c(cls, class(deep_variant_data[[dataset_name]]))
  }
}

# Save generated dataset
write_rds(deep_variant_data, 'data/variant_data.RDS')

