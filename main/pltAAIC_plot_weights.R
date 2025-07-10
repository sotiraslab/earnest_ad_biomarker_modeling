
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

HAUFE_FOLDER <- 'outputs/haufe_weights_aaic/'

biomarker.path <- file.path(HAUFE_FOLDER, 'svm_weights_haufe_ATNbiomarker.csv')
roi.path <- file.path(HAUFE_FOLDER, 'svm_weights_haufe_AmyloidROI.csv')

OUTPUT <- file.path('figures', 'haufe_brains_aaic')
dir.create(OUTPUT, showWarnings = F)

# --- ATN biomarkers barplot -----

df <- read.csv(biomarker.path)

p.data <- df %>%
  pivot_longer(everything(), names_to = 'Feature', values_to = 'Importance') %>%
  group_by(Feature) %>%
  summarise(Mean = mean(abs(Importance)),
            SD = sd(abs(Importance))) %>%
  ungroup() %>%
  recode(Feature,
         AmyloidComposite = 'Amyloid SUVR',
         META_TEMPORAL_TAU = 'Tau SUVR',
         META_TEMPORAL_VOL = 'Temporal Volume',
         Age = 'Age',
         SexBinary = 'Sex',
         HasE4Binary = 'APOEE4')
