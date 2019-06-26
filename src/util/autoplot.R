#!/usr/bin/env Rscript 
# Convienience functions and objects for automated plot creation

# ggplot object with attached saving parameters (additional named arguments passed to ggsave)
labeled_ggplot <- function(p, units='in', ...){
  l <- list(plot = p, units = units, ...)
  class(l) <- 'labeled_ggplot'
  return(l)
}

# smart save plots based on type (arguments passed through)
smart_save <- function(p, ...){
  UseMethod('smart_save', p)
}

smart_save.ggplot <- function(p, filename, override=FALSE, ...){
  ggsave(filename, plot = p, ...)
}

smart_save.ggmatrix <- function(p, filename, override=FALSE, ...){
  ggsave(filename, plot = p, ...)
}

smart_save.labeled_ggplot <- function(p, filename, override=FALSE, ...){
  if (override){
    params <- c(list(...), p, formals(ggsave))
  } else {
    params <- c(p, list(...), formals(ggsave))
  }
   
  ggsave(filename = filename, plot = params$plot, device = params$device, path = params$path, scale = params$scale,
         width = params$width, height = params$height, units = params$units, dpi = params$dpi, limitsize = params$limitsize)
}

# Save named nested lists of plots
save_plot_list <- function(plotlist, root, listname='', overwrite=FALSE, level=0, report_level=1){
  if (str_sub(root, -1) == '/'){
    root <- str_sub(root, end=-2)
  }
  
  if (level > 0 & level <= report_level & nchar(listname) > 0){
    cat(str_c(rep('--', level-1), collapse=''), 'Writing ', listname, sep = '')
  } else if (level == report_level + 1){
    cat(' -', listname)
  }
  
  dir.create(root, showWarnings = FALSE, recursive = TRUE)
  
  for (name in names(plotlist)){
    if (class(plotlist[[name]])[1] == 'list'){
      # Recursively process sublists
      save_plot_list(plotlist[[name]], root=str_c(root, '/', name), listname=name, overwrite=FALSE, level=level+1)
      
    } else {
      # Assume everything else is a plot, and try to save
      plot_path <- str_c(root, '/', name, '.pdf')
      # Only write figures that don't exist, delete figs to regenerate
      if (overwrite | !file.exists(plot_path)){
        smart_save(plotlist[[name]], plot_path, override=FALSE, width = 7, height = 5)
      }
    }
  }
  
  if (level > 0 & level <= report_level & nchar(listname) > 0){
    cat('\n')
  }
}