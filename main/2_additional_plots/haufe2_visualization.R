
# --- Imports -----
sh <- suppressPackageStartupMessages

sh(library(showtext))
sh(library(stringr))
sh(library(this.path))
sh(library(tidyverse))

# --- Working directory -----
setwd(this.dir())

# --- Source plot script -----
source('../../scripts/ggseg_plots.R')

# --- Load SVM weights ------

INFOLDER.LONG <- '../../outputs/exp1_svms_global_cognition'
INFOLDER.SHORT <- '../../outputs/exp1_svms_global_cognition_short/'

if (file.exists(file.path(INFOLDER.LONG, 'results.csv'))) {
  INFOLDER <- INFOLDER.LONG
} else {
  print('!!!!!!!!!!!!!!')
  print('WARNING: Using short results for SVM weights.')
  print('!!!!!!!!!!!!!!')
  INFOLDER <- INFOLDER.SHORT
}

amy.path <- file.path(INFOLDER, 'svm_weights_haufe_AmyloidSVM.csv')
tau.path <- file.path(INFOLDER, 'svm_weights_haufe_TauSVM.csv')
gm.path <- file.path(INFOLDER, 'svm_weights_haufe_GMSVM.csv')
atn.path <- file.path(INFOLDER, 'svm_weights_haufe_ATNSVM.csv')

OUTPUT <- '../../outputs/additional_plots/'
dir.create(OUTPUT, showWarnings = F)

# --- Add font -------

font_add('arial', '../../fonts/arial.ttf')

# --- Figures from combined SVM -------
df <- read.csv(atn.path)
all.cols <- colnames(df)

amy.cols <- all.cols %>%
  str_subset('^AV45_')

tau.cols <- all.cols %>%
  str_subset('^FTP_')

gm.cols <- all.cols %>%
  str_subset('_VOLUME$')

region.labels <- amy.cols %>%
  str_replace('AV45_', '') %>%
  adni.labels.to.ggseg()

all.rois <- c(amy.cols, tau.cols, gm.cols)

mean.weights <- colMeans(abs(df))
maxi <- max(mean.weights)
mini <- 0

# amyloid
plot.cortex(values = mean.weights[amy.cols],
            regions = region.labels,
            vmin = mini,
            vmax = maxi,
            name = 'Amyloid') +
  theme(text = element_text(family = 'arial'))
ggsave(file.path(OUTPUT, 'combined_haufe_weights_amyloid.png'), width = 6, height = 1.8, units='in') +
  theme(text = element_text(family = 'arial'))

plot.subcortex(values = mean.weights[amy.cols],
               regions = region.labels,
               vmin = mini,
               vmax = maxi,
               name = '',
               legend=F) +
  theme(text = element_text(family = 'arial'))
ggsave(file.path(OUTPUT, 'combined_haufe_weights_amyloid_subcortical.png'), width = 3, height = 1.8, units='in') +
  theme(text = element_text(family = 'arial'))

# tau
plot.cortex(values = mean.weights[tau.cols],
            regions = region.labels,
            vmin = mini,
            vmax = maxi,
            name = 'Tau') +
  theme(text = element_text(family = 'arial'))
ggsave(file.path(OUTPUT, 'combined_haufe_weights_tau.png'), width = 6, height = 1.8, units='in')

plot.subcortex(values = mean.weights[tau.cols],
               regions = region.labels,
               vmin = mini,
               vmax = maxi,
               name = '',
               legend=F) +
  theme(text = element_text(family = 'arial'))
ggsave(file.path(OUTPUT, 'combined_haufe_weights_tau_subcortical.png'), width = 3, height = 1.8, units='in') +
  theme(text = element_text(family = 'arial'))

# GM
plot.cortex(values = mean.weights[gm.cols],
            regions = region.labels,
            vmin = mini,
            vmax = maxi,
            name = 'Neurodegeneration') +
  theme(text = element_text(family = 'arial'))
ggsave(file.path(OUTPUT, 'combined_haufe_weights_gm.png'), width = 6, height = 1.8, units='in')

plot.subcortex(values = mean.weights[gm.cols],
               regions = region.labels,
               vmin = mini,
               vmax = maxi,
               name = '',
               legend=F) +
  theme(text = element_text(family = 'arial'))
ggsave(file.path(OUTPUT, 'combined_haufe_weights_gm_subcortical.png'), width = 3, height = 1.8, units='in') +
  theme(text = element_text(family = 'arial'))

# --- Figures from unimodal SVM -------

plot.unimodal <- function(biomarker) {
  paths <- list(amyloid=amy.path,
                tau=tau.path,
                gm=gm.path)
  cols <- list(amyloid=amy.cols,
              tau=tau.cols,
              gm=gm.cols)
  selected.path <- paths[[biomarker]]
  selected.cols <- cols[[biomarker]]
  if (is.null(selected.path)) {
    stop(sprintf('Cannot recognize biomarker "%s"', biomarker))
  }
  
  df <- read.csv(selected.path) %>%
    select(all_of(selected.cols))
  mean.weights <- colMeans(abs(df))
  mini <- 0
  maxi <- max(mean.weights)
  
  plot.cortex(values = mean.weights,
              regions = region.labels,
              name = biomarker,
              vmin = mini,
              vmax = maxi) +
    theme(text = element_text(family = 'arial'))
  filename <- sprintf('unimodal_haufe_weights_%s.png', biomarker)
  ggsave(file.path(OUTPUT, filename), width = 6, height = 1.8, units='in')
  
  plot.subcortex(values = mean.weights,
                 regions = region.labels,
                 name = '',
                 legend=F,
                 vmin = mini,
                 vmax = maxi) +
    theme(text = element_text(family = 'arial'))
  filename <- sprintf('unimodal_haufe_weights_%s_subcortical.png', biomarker)
  ggsave(file.path(OUTPUT, filename), width = 6, height = 1.8, units='in')
}

plot.unimodal('amyloid')
plot.unimodal('tau')
plot.unimodal('gm')