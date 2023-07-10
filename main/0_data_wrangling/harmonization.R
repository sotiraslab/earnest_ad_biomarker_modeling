# === imports ======

sh <- suppressPackageStartupMessages

sh(library(colormap))
sh(library(ggseg))
sh(library(gridExtra))
sh(library(gsubfn))
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

adni$CDRFactor <- factor(adni$CDR)
oasis$CDRFactor <- factor(oasis$CDR)

# === harmonization function =======

harmonize <- function(data.a, data.b, columns, covariates) {
  a <- t(as.matrix(data.a[, columns]))
  b <- t(as.matrix(data.b[, columns]))
  
  combine.data <- rbind(data.a[, covariates], data.b[, covariates])
  
  fml <- as.formula(paste('~', paste(covariates, collapse='+')))
  
  dat <- cbind(a, b)
  batch <- c(rep(1, ncol(a)), rep(2, ncol(b)))
  mod <- model.matrix(fml, data=combine.data)
  
  harmonized <- neuroCombat(dat=dat, batch=batch, mod=mod)
  harmonized.data <- harmonized$dat.combat
  
  output <- list('harmonized.data.a' = as.data.frame(t(harmonized.data[, batch == 1])),
                 'harmonized.data.b' = as.data.frame(t(harmonized.data[, batch == 2])),
                 'columns' = columns,
                 'covariates' = covariates)
  
  return (output)
}

# === harmonize =======

av45.rois <- colnames(adni)[str_detect(colnames(adni), '^AV45')]
ftp.rois <- colnames(adni)[str_detect(colnames(adni), '^FTP')]
gm.rois <- colnames(adni)[str_detect(colnames(adni), '_VOLUME')]

covariates <- c('Age', 'Sex', 'HasE4', 'CDRFactor')

av45.harmonize <- harmonize(adni, oasis, av45.rois, covariates = covariates)
ftp.harmonize <- harmonize(adni, oasis, ftp.rois, covariates = covariates)
gm.harmonize <- harmonize(adni, oasis, gm.rois, covariates = covariates)

adni.harmonized <- adni
adni.harmonized[, av45.rois] <- av45.harmonize$harmonized.data.a
adni.harmonized[, ftp.rois] <- ftp.harmonize$harmonized.data.a
adni.harmonized[, gm.rois] <- gm.harmonize$harmonized.data.a

oasis.harmonized <- oasis
oasis.harmonized[, av45.rois] <- av45.harmonize$harmonized.data.b
oasis.harmonized[, ftp.rois] <- ftp.harmonize$harmonized.data.b
oasis.harmonized[, gm.rois] <- gm.harmonize$harmonized.data.b

# === create plotting functions ========


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


plot.ggseg <- function(values, regions, vmin = NULL, vmax = NULL,
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
               position = position_brain(hemi ~ side),
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
      scale_fill_colormap(colormap=cm, limits=c(vmin, vmax), oob = scales::squish)
  }
  
  # subcortical
  atlas <- data.frame(label = brain_labels(aseg)) %>%
    left_join(df, by='label') %>%
    drop_na()
  
  plot.subcortical <- ggplot(atlas) +
    geom_brain(atlas = aseg,
               color='black',
               size=.5,
               aes(fill = value)) +
    theme_void()
  
  if (cm == 'DIVERGE') {
    plot.subcortical <- plot.subcortical + 
      scale_fill_gradient2(low='blue', high='red', mid='white', limits=c(vmin, vmax), oob = scales::squish)
  } else {
    plot.subcortical <- plot.subcortical + 
      scale_fill_colormap(colormap=cm, limits=c(vmin, vmax), oob = scales::squish)
  }

  if (! legend) plot.subcortical <- plot.subcortical + theme(legend.position = 'none')
  
  final <- grid.arrange(plot.cortical, plot.subcortical, ncol=2)
  
  return(final)
}

plot.helper <- function(adni.mask, oasis.mask, modality) {

  if (modality == 'av45') {
    cols <- av45.rois
    regions <- adni.labels.to.ggseg(gsub('AV45_', '', av45.rois))
  } else if (modality == 'ftp') {
    cols <- ftp.rois
    regions <- adni.labels.to.ggseg(gsub('FTP_', '', ftp.rois))
  } else if (modality == 'gm') {
    cols <- gm.rois
    regions <- adni.labels.to.ggseg(gm.rois)
  } else {
    stop('modality must be "av45", "ftp", or "gm"')
  }

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

  a <- plot.ggseg(adni.before, regions, vmin=vmin, vmax=vmax, name='ADNI (original)')
  b <- plot.ggseg(oasis.before, regions, vmin=vmin, vmax=vmax, name='OASIS (original)')
  c <- plot.ggseg(diff.before, regions, vmin=dmin, vmax=dmax, name='Difference (original)', cm='DIVERGE')

  d <- plot.ggseg(adni.after, regions, vmin=vmin, vmax=vmax, name='ADNI (ComBat)')
  e <- plot.ggseg(oasis.after, regions, vmin=vmin, vmax=vmax, name='OASIS (ComBat)')
  f <- plot.ggseg(diff.after, regions, vmin=dmin, vmax=dmax, name='Difference (ComBat)', cm='DIVERGE')

  final <- grid.arrange(a, d, b, e, c, f, nrow = 3, ncol = 2)

  return (final)
}

ttest.df <- function(a, b, cols = NULL, zero.nonsig = T) {
  if (is.null(cols)) cols <- colnames(a)
  ts <- abs(sapply(cols, function (x) t.test(a[x], b[x])$statistic))
  ps <- sapply(cols, function (x) t.test(a[x], b[x])$p.value)
  if (zero.nonsig) {
    ts[ps > 0.05] <- NA
  }
  df <- data.frame(t=ts, p=ps)
  return (df)
}

combat.evaluation <- function(modality, adni.mask, oasis.mask) {
  
  if (modality == 'av45') {
    cols <- av45.rois
    regions <- adni.labels.to.ggseg(gsub('AV45_', '', av45.rois))
  } else if (modality == 'ftp') {
    cols <- ftp.rois
    regions <- adni.labels.to.ggseg(gsub('FTP_', '', ftp.rois))
  } else if (modality == 'gm') {
    cols <- gm.rois
    regions <- adni.labels.to.ggseg(gm.rois)
  } else {
    stop('modality must be "av45", "ftp", or "gm"')
  }
  
  adni.before <- adni[adni.mask, cols]
  oasis.before <- oasis[oasis.mask, cols]
  diff.before <- colMeans(adni.before) - colMeans(oasis.before)

  adni.after <- adni.harmonized[adni.mask, cols]
  oasis.after <- oasis.harmonized[oasis.mask, cols]
  diff.after <- colMeans(adni.after) - colMeans(oasis.after)
  
  t.before <- ttest.df(adni.before, oasis.before, zero.nonsig = T)$t
  t.after <- ttest.df(adni.after, oasis.after, zero.nonsig = T)$t
  
  dlim <- max(abs(diff.before), abs(diff.after))
  tmax <- max(t.before, t.after, na.rm = T)
  tmin <- 0
  
  a <- plot.ggseg(diff.before, regions, vmin=-dlim, vmax=dlim, name='Difference (original)', cm='DIVERGE')
  b <- plot.ggseg(diff.after, regions, vmin=-dlim, vmax=dlim, name='Difference (ComBat)', cm='DIVERGE')
  c <- plot.ggseg(t.before, regions, vmin=tmin, vmax=tmax, name='TTest (before)', cm='viridis')
  d <- plot.ggseg(t.after, regions, vmin=tmin, vmax=tmax, name='TTest (ComBat)', cm='viridis')
  
  final <- grid.arrange(a, b, c, d, nrow = 2, ncol = 2)
  
  return (final)
}

# === plots =====

for (m in c('av45', 'ftp', 'gm')) {
  p <- combat.evaluation(m, adni$CDRFactor == 0 & adni$AmyloidPositive == 0, oasis$CDRFactor == 0 & oasis$AmyloidPositive == 0)
  ggsave(sprintf('%s_harmonize_cn.png', m), p, width=8, height=4, units='in')
  
  p <- combat.evaluation(m, adni$CDRFactor == 0 & adni$AmyloidPositive == 1, oasis$CDRFactor == 0 & oasis$AmyloidPositive == 1)
  ggsave(sprintf('%s_harmonize_preclinical.png', m), p, width=8, height=4, units='in')
  
  p <- combat.evaluation(m, adni$CDRFactor == 0.5,  oasis$CDRFactor == 0.5)
  ggsave(sprintf('%s_harmonize_mci.png', m), p, width=8, height=4, units='in')
  
  p <- combat.evaluation(m, adni$CDRFactor == 1.0,  oasis$CDRFactor == 1.0)
  ggsave(sprintf('%s_harmonize_ad.png', m), p, width=8, height=4, units='in')
}

# === save data =======

write.csv(adni.harmonized, '../../data/derivatives/adni_harmonized.csv', row.names = F, quote = F, na = '')
write.csv(oasis.harmonized, '../../data/derivatives/oasis_harmonized.csv', row.names = F, quote = F, na = '')