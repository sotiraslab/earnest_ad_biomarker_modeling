
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
roi.path <- file.path(HAUFE_FOLDER, 'svm_weights_haufe_ATNROI.csv')

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
  mutate(
    Feature = recode(
      Feature, 
      AmyloidComposite = 'Amyloid SUVR',
      META_TEMPORAL_TAU = 'Tau SUVR',
      META_TEMPORAL_VOL = 'Temporal Volume',
      Age = 'Age',
      SexBinary = 'Sex',
      HasE4Binary = 'APOEE4'),
    Feature = factor(Feature, levels=c('Amyloid SUVR', 'Tau SUVR', 'Temporal Volume', 'Age', "Sex", 'APOEE4'))
         )

ggplot(p.data, aes(x=Feature, y=Mean, fill=Feature)) +
  geom_bar(stat='identity') + 
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width=0.2) +
  theme_bw() +
  scale_y_continuous(expand=expansion(c(0, .1), 0)) +
  scale_fill_manual(
    values = c(
      'Amyloid SUVR'='#fe6100',
      'Tau SUVR'='#dc267f',
      'Temporal Volume'='#648fff'
    )
  ) +
  theme(legend.position = 'none', text=element_text(size=15)) +
  ylab('Model Importance')

ggsave(file.path(OUTPUT, 'atn_biomarkers_importance.svg'), width = 8, height = 4, units='in')

# --- ATN regional brains -----

df <- read.csv(roi.path)

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

# colormap
cm = 'viridis' 

# amyloid
plot.cortex(values = mean.weights[amy.cols],
            regions = region.labels,
            vmin = mini,
            vmax = maxi,
            cm = cm)
ggsave(file.path(OUTPUT, 'combined_haufe_weights_amyloid.svg'), width = 6.5, height = 1.2, units='in')

plot.subcortex(values = mean.weights[amy.cols],
               regions = region.labels,
               vmin = mini,
               vmax = maxi,
               legend=F,
               cm = cm)
ggsave(file.path(OUTPUT, 'combined_haufe_weights_amyloid_subcortical.svg'), width = 3, height = 1.2, units='in')

# tau
plot.cortex(values = mean.weights[tau.cols],
            regions = region.labels,
            vmin = mini,
            vmax = maxi,
            cm = cm)
ggsave(file.path(OUTPUT, 'combined_haufe_weights_tau.svg'), width = 6.5, height = 1.2, units='in')

plot.subcortex(values = mean.weights[tau.cols],
               regions = region.labels,
               vmin = mini,
               vmax = maxi,
               cm = cm)
ggsave(file.path(OUTPUT, 'combined_haufe_weights_tau_subcortical.svg'), width = 3, height = 1.2, units='in')


# GM
plot.cortex(values = mean.weights[gm.cols],
            regions = region.labels,
            vmin = mini,
            vmax = maxi,
            cm = cm)
ggsave(file.path(OUTPUT, 'combined_haufe_weights_gm.svg'), width = 6.5, height = 1.2, units='in')

plot.subcortex(values = mean.weights[gm.cols],
               regions = region.labels,
               vmin = mini,
               vmax = maxi,
               cm = cm)
ggsave(file.path(OUTPUT, 'combined_haufe_weights_gm_subcortical.svg'), width = 3, height = 1.8, units='in')

# colorbar
plot.colorbar(mini, maxi, cm = cm, text.size = 20, 
              savepath = file.path(OUTPUT, 'combined_haufe_weights_colorbar.svg'))