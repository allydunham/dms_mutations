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

smart_save.labeled_ggplot <- function(p, filename, override=FALSE, ...){
  if (override){
    params <- c(list(...), p, formals(ggsave))
  } else {
    params <- c(p, list(...), formals(ggsave))
  }
   
  ggsave(filename = filename, plot = params$plot, device = params$device, path = params$path, scale = params$scale,
         width = params$width, height = params$height, units = params$units, dpi = params$dpi, limitsize = params$limitsize)
}