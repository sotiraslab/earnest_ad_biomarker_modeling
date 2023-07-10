# === imports ======

sh <- suppressPackageStartupMessages

sh(library(ADNIMERGE))
sh(library(colormap))
sh(library(ggseg))
sh(library(gridExtra))
sh(library(lubridate))
sh(library(neuroCombat))
sh(library(stringr))
sh(library(this.path))
sh(library(tidyverse))

# === Set working directory ======

setwd(this.dir())

# === Set paths ======

PATH.ADNI <- '../../data/derivatives/adni_base_table.csv'
PATH.OASIS <- '../../data/derivatives/oasis_base_table.csv'

# === read data =======

adni <- read.csv(PATH.ADNI)
oasis <- read.csv(PATH.OASIS)

# === volume harmonization =======

gm.rois <- colnames(adni)[str_detect(colnames(adni), '_VOLUME$')]
gm.rois <- gm.rois[! str_detect(gm.rois, 'META_TEMPORAL|BRAAK')]

adni.mat <- t(as.matrix(adni[, gm.rois]))
oasis.mat <- t(as.matrix(oasis[, gm.rois]))

# design matrix
adni$CDRFactor <- factor(adni$CDRGlobal)
oasis$CDRFactor <- factor(oasis$CDR)
covars <- c('Age', 'Sex', 'HasE4', 'CDRFactor')
combine.data <- rbind(
  adni[, covars],
  oasis[, covars]
)

dat <- cbind(adni.mat, oasis.mat)
batch <- c(rep(1, ncol(adni.mat)), rep(2, ncol(oasis.mat)))
mod <- model.matrix(~Age+Sex+HasE4+CDRFactor, data=combine.data)

harmonized <- neuroCombat(dat=dat, batch=batch, mod=mod)
harmonized.data <- harmonized$dat.combat

adni.harmonized <- adni
adni.harmonized[, gm.rois] <- as.data.frame(t(harmonized.data[, batch == 1]))

oasis.harmonized <- oasis
oasis.harmonized[, gm.rois] <- as.data.frame(t(harmonized.data[, batch == 2]))

# === create plotting functions ========

labels.to.ggseg <- function(labels) {
  
  # cortical
  cortical <- tolower(labels) %>%
    str_replace_all('_volume$', '')
  cortical <- paste('lh_', cortical, sep='')
  is.cortical <- cortical %in% brain_labels(dk)
  
  # subcortical
  subcortical <- tolower(labels) %>%
    str_replace_all('_volume$', '') %>%
    str_to_title() %>%
    str_replace_all('Ventraldc', 'VentralDC') %>%
    str_replace_all('Thalamus_proper', 'Thalamus-Proper')
  subcortical <- paste('Left-', subcortical, sep = '')
  is.subcortical <- subcortical %in% brain_labels(aseg)
  
  return (c(cortical[is.cortical], subcortical[is.subcortical]))
}

plot.unilateral.ggseg <- function(values, regions, vmin = NULL, vmax = NULL,
                                  name = 'My Plot', legend = T, cm = 'inferno') {
  if (is.null(vmin)) vmin <- min(values)
  if (is.null(vmax)) vmax <- max(values)
  
  df <- data.frame(label=regions, value=values)
  
  # cortical
  atlas <- data.frame(label = brain_labels(dk)) %>%
    left_join(df, by='label') %>%
    drop_na()
  
  plot.cortical <- ggplot(atlas) +
    geom_brain(atlas = dk,
               hemi = 'left',
               color='black',
               size=.5,
               aes(fill = value)) +
    theme_void() +
    theme(legend.position = 'none') +
    ggtitle(name)
  
  if (cm == 'DIVERGE') {
    plot.cortical <- plot.cortical + 
      scale_fill_gradient2(low='blue', high='red', mid='white', limits=c(vmin, vmax), oob = scales::squish)
  } else {
    plot.cortical <- plot.cortical + 
      scale_fill_colormap(colormap='inferno', limits=c(vmin, vmax), oob = scales::squish)
  }
  
  # subcortical
  atlas <- data.frame(label = brain_labels(aseg)) %>%
    left_join(df, by='label') %>%
    drop_na()
  
  plot.subcortical <- ggplot(atlas) +
    geom_brain(atlas = aseg,
               hemi = 'left',
               color='black',
               size=.5,
               aes(fill = value)) +
    theme_void()
  
  if (cm == 'DIVERGE') {
    plot.subcortical <- plot.subcortical + 
      scale_fill_gradient2(low='blue', high='red', mid='white', limits=c(vmin, vmax), oob = scales::squish)
  } else {
    plot.subcortical <- plot.subcortical + 
      scale_fill_colormap(colormap='inferno', limits=c(vmin, vmax), oob = scales::squish)
  }

  if (! legend) plot.subcortical <- plot.subcortical + theme(legend.position = 'none')
  
  final <- grid.arrange(plot.cortical, plot.subcortical, ncol=2)
  
  return(final)
}

plot.helper <- function(adni.mask, oasis.mask, cols = gm.rois) {
  
  adni.before <- colMeans(adni[adni.mask, cols])
  oasis.before <- colMeans(oasis[oasis.mask, cols])
  diff.before <- adni.before - oasis.before
  
  adni.after <- colMeans(adni.harmonized[adni.mask, cols])
  oasis.after <- colMeans(oasis.harmonized[oasis.mask, cols])
  diff.after <- adni.after - oasis.after
  
  vmax <- max(adni.before, oasis.before, adni.after, oasis.after)
  vmin <- min(adni.before, oasis.before, adni.after, oasis.after)
  dmax <- max(diff.before, diff.after)
  dmin <- min(diff.before, diff.after)
  
  regions <- labels.to.ggseg(cols)
  a <- plot.unilateral.ggseg(adni.before, regions, vmin=vmin, vmax=vmax, name='ADNI (original)')
  b <- plot.unilateral.ggseg(oasis.before, regions, vmin=vmin, vmax=vmax, name='OASIS (original)')
  c <- plot.unilateral.ggseg(diff.before, regions, vmin=dmin, vmax=dmax, name='Difference (original)', cm='DIVERGE')
  
  d <- plot.unilateral.ggseg(adni.after, regions, vmin=vmin, vmax=vmax, name='ADNI (ComBat)')
  e <- plot.unilateral.ggseg(oasis.after, regions, vmin=vmin, vmax=vmax, name='OASIS (ComBat)')
  f <- plot.unilateral.ggseg(diff.after, regions, vmin=dmin, vmax=dmax, name='Difference (ComBat)', cm='DIVERGE')
  
  final <- grid.arrange(a, d, b, e, c, f, nrow = 3, ncol = 2)
  
  return (final)
}

# === plots =====

p <- plot.helper(adni$CDRFactor == 0 & adni$AmyloidPositive == 0,
                 oasis$CDRFactor == 0 & oasis$AmyloidPositive == 0)
ggsave('gm_harmonization_cn.png', p, width = 8, height = 6, units = 'in')

p <- plot.helper(adni$CDRFactor == 0.5, oasis$CDRFactor == 0.5)
ggsave('gm_harmonization_cdr0.5.png', p, width = 8, height = 6, units = 'in')

p <- plot.helper(adni$CDRFactor == 1.0, oasis$CDRFactor == 1.0)
ggsave('gm_harmonization_cdr1.0.png', p, width = 8, height = 6, units = 'in')

# === save data =======

write.csv(adni.harmonized, '../../data/derivatives/adni_harmonized.csv', row.names = F, quote = F, na = '')
write.csv(oasis.harmonized, '../../data/derivatives/oasis_harmonized.csv', row.names = F, quote = F, na = '')