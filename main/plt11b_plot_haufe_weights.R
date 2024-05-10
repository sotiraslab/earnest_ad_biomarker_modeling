
# --- Imports -----
sh <- suppressPackageStartupMessages

sh(library(stringr))
sh(library(svglite))
sh(library(this.path))
sh(library(tidyverse))

# --- Working directory -----
setwd(this.dir())

# --- Source plot script -----
source('../rscripts/ggseg_plots.R')

# --- Load SVM weights ------

HAUFE_FOLDER <- 'outputs/haufe_weights/'

amy.path <- file.path(HAUFE_FOLDER, 'svm_weights_haufe_AmyloidSVM.csv')
tau.path <- file.path(HAUFE_FOLDER, 'svm_weights_haufe_TauSVM.csv')
gm.path <- file.path(HAUFE_FOLDER, 'svm_weights_haufe_GMSVM.csv')
atn.path <- file.path(HAUFE_FOLDER, 'svm_weights_haufe_ATNSVM.csv')

OUTPUT <- file.path('figures', 'haufe_brains')
dir.create(OUTPUT, showWarnings = F)

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
            name = 'Amyloid')
ggsave(file.path(OUTPUT, 'combined_haufe_weights_amyloid.svg'), width = 6, height = 1.8, units='in')

plot.subcortex(values = mean.weights[amy.cols],
               regions = region.labels,
               vmin = mini,
               vmax = maxi,
               name = '',
               legend=F)
ggsave(file.path(OUTPUT, 'combined_haufe_weights_amyloid_subcortical.svg'), width = 3, height = 1.8, units='in')

# tau
plot.cortex(values = mean.weights[tau.cols],
            regions = region.labels,
            vmin = mini,
            vmax = maxi,
            name = 'Tau')
ggsave(file.path(OUTPUT, 'combined_haufe_weights_tau.svg'), width = 6, height = 1.8, units='in')

plot.subcortex(values = mean.weights[tau.cols],
               regions = region.labels,
               vmin = mini,
               vmax = maxi,
               name = '',
               legend=F)
ggsave(file.path(OUTPUT, 'combined_haufe_weights_tau_subcortical.svg'), width = 3, height = 1.8, units='in')

# GM
plot.cortex(values = mean.weights[gm.cols],
            regions = region.labels,
            vmin = mini,
            vmax = maxi,
            name = 'Neurodegeneration')
ggsave(file.path(OUTPUT, 'combined_haufe_weights_gm.svg'), width = 6, height = 1.8, units='in')

plot.subcortex(values = mean.weights[gm.cols],
               regions = region.labels,
               vmin = mini,
               vmax = maxi,
               name = '',
               legend=F)
ggsave(file.path(OUTPUT, 'combined_haufe_weights_gm_subcortical.svg'), width = 3, height = 1.8, units='in')

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
              vmax = maxi)
  filename <- sprintf('unimodal_haufe_weights_%s.svg', biomarker)
  ggsave(file.path(OUTPUT, filename), width = 6, height = 1.8, units='in')
  
  plot.subcortex(values = mean.weights,
                 regions = region.labels,
                 name = '',
                 legend=F,
                 vmin = mini,
                 vmax = maxi)
  filename <- sprintf('unimodal_haufe_weights_%s_subcortical.svg', biomarker)
  ggsave(file.path(OUTPUT, filename), width = 6, height = 1.8, units='in')
}

plot.unimodal('amyloid')
plot.unimodal('tau')
plot.unimodal('gm')