
library(colormap)
library(dplyr)
library(ggseg)
library(gsubfn)
library(tidyverse)

check.in.ggseg <- function(label) {
  return(label %in% brain_labels(dk) | label %in% brain_labels(aseg))
}

adni.labels.to.ggseg <- function(labels, filter.failures = T,
                                 warn = T) {
  
  # cortical
  cortical <- labels %>%
    str_replace('^CTX_', '') %>%
    str_replace('_VOLUME$', '') %>%
    str_replace('_SUVR$', '') %>%
    tolower()
  
  # subcortical
  subcortical <- labels %>%
    str_replace('_SUVR$', '') %>%
    str_replace('_VOLUME$', '') %>%
    str_replace_all('_', '-') %>%
    tolower() %>%
    str_replace('dc$', 'DC')
  subcortical <- gsubfn('^.|-.', toupper, subcortical)
  
  # assemble output
  output <- ifelse(str_detect(labels, '^CTX'), cortical, subcortical)
  
  # check success
  failures <- ! check.in.ggseg(output)
  
  # filter
  if (filter.failures) {
    output <- ifelse(failures, labels, output)
  }
  
  # warning
  if (warn) {
    issues <- output[failures]
    if (length(issues) > 0) {
      s <- paste(issues, collapse=', ')
      print(sprintf('WARNING! labels not matching ggseg: %s', s))
    }
  }
  
  return(output)
}

plot.cortex<- function(values, regions, vmin = NULL, vmax = NULL,
                       name = NULL, legend = F, cm = 'inferno',
                       cm.reverse = F) {
  if (is.null(vmin)) vmin <- min(values)
  if (is.null(vmax)) vmax <- max(values)
  
  df <- data.frame(label=regions, value=values)
  
  # cortical
  atlas <- data.frame(label = brain_labels(dk)) %>%
    left_join(df, by='label') %>%
    drop_na()
  
  plot.cortical <- ggplot(atlas) +
    geom_brain(atlas = dk,
               color='black',
               size=.5,
               aes(fill = value)) +
    theme_void() +
    scale_fill_colormap(colormap=cm, limits=c(vmin, vmax), oob = scales::squish, reverse = cm.reverse)
  
  if (! legend) plot.cortical <- plot.cortical + theme(legend.position = 'none')
  if (! is.null(name)) plot.cortical <- plot.cortical + ggtitle(name)
  
  return(plot.cortical)
}

plot.subcortex <- function(values, regions, vmin = NULL, vmax = NULL,
                           name = NULL, legend = F, cm = 'inferno',
                           cm.reverse = F) {
  if (is.null(vmin)) vmin <- min(values)
  if (is.null(vmax)) vmax <- max(values)
  
  df <- data.frame(label=regions, value=values)
  
  atlas <- data.frame(label = brain_labels(aseg)) %>%
    left_join(df, by='label') %>%
    drop_na()
  
  plot.subcortical <- ggplot(atlas) +
    geom_brain(atlas = aseg,
               color='black',
               size=.5,
               aes(fill = value)) +
    theme_void() +
    scale_fill_colormap(colormap=cm, limits=c(vmin, vmax), oob = scales::squish, reverse = cm.reverse)
  
  if (! legend) plot.subcortical <- plot.subcortical + theme(legend.position = 'none')
  if (! is.null(name)) plot.subcortical <- plot.subcortical + ggtitle(name)
  
  return(plot.subcortical)
}

plot.colorbar <- function(mini, maxi, cm,
                          cm.reverse = F,
                          density = 10000,
                          ticks = 4 ,
                          width = 12,
                          height = 1,
                          savepath = NULL,
                          text.size = 16) {
  breaks = round(seq(mini, maxi, length.out = ticks), 2)
  
  data <- data.frame(x = seq(mini, maxi, length.out=density),
                     y = rep(1, density),
                     intensity = seq(mini, maxi, length.out=density))
  
  ggplot(data, aes(x=x, y=y, fill=intensity, color=intensity)) + 
    geom_raster() +
    scale_fill_colormap(colormap = cm, reverse = cm.reverse) + 
    scale_color_colormap(colormap = cm, reverse = cm.reverse) + 
    scale_x_continuous(expand=c(0,0), breaks=breaks) +
    scale_y_continuous(expand=c(0,0)) +
    xlab('') + 
    ylab ('') +
    theme(legend.position = 'none',
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          text = element_text(size=text.size))
  
  if (! is.null(savepath)) {
    ggsave(savepath, width = width, height = height)
  }
}