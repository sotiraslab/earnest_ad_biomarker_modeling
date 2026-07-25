
# --- Imports -----
sh <- suppressPackageStartupMessages

sh(library(stringr))
sh(library(this.path))
sh(library(tidyverse))

# --- Working directory -----
setwd(this.dir())

# --- Load data -----

df <- read.csv('datasets/maindata.csv')
predictions <- read.csv('outputs/individual_predictions.csv')
mmse <- read.csv('datasets/mmse_lmem_fit_data.csv')

last.mmse <- mmse %>%
  group_by(ID) %>%
  slice_tail(n=1) %>%
  ungroup() %>%
  filter(ID %in% df$RID)

# --- Helper function ------

p.annotation <- function(p, text=T, stars=T) {
  stars.labels <- c('***', "**", "*", "")
  text.labels <- c('<0.001', "<0.01", "<0.05", "n.s.")
  
  labels <- rep('', 4)
  if (text) {
    labels <-str_c(labels, text.labels)
  }
  if (stars) {
    labels <- str_c(labels, stars.labels)
  }
  
  
  annot <- cut(p,
               breaks = c(0, 0.001, 0.01, 0.05, Inf),
               labels = labels,
               include.lowest = T)
  annot <- as.character(annot)
  return (annot)
}

# --- Compare LMEM residual to prediction accuracy -------

odir <- file.path('figures', 'residual_vs_lmemSE')
dir.create(odir, showWarnings = F)

cols <- colnames(predictions)
model_keys <- cols[str_detect(cols, '_residual')]
# model_keys <- str_replace(model_keys, '_residual', '')

reg.stats <- matrix(data = NA, nrow=length(model_keys), ncol = 4)

for (i in 1:length(model_keys)) {
  x <- predictions$DeltaMMSEBootSE
  y <- predictions[[model_keys[i]]]
  title <- model_keys[[i]]
  title <- str_replace_all(title, '_residual', '')
  title <- str_replace_all(title, '\\.', ' ')
  
  pdata <- data.frame(x=x, y=y)
  p <- ggplot(pdata, aes(x=x, y=y)) +
    geom_point() +
    geom_smooth(method='lm') +
    theme_bw() +
    ggtitle(title) +
    xlab('Standard error of LMEM fit') +
    ylab('Prediction residual for ΔMMSE')
  # print(p)
  
  m <- lm(y ~ x, data = pdata)
  m.sum <- summary(m)
  f <- m.sum$fstatistic
  p <- pf(f[1], f[2], f[3], lower.tail = F)
  attributes(p) <- NULL
  
  reg.stats[i, 1] <- title
  reg.stats[i, 2] <- round(m.sum$r.squared, 3)
  reg.stats[i, 3] <- p.annotation(p)
  reg.stats[i, 4] <- round(m.sum$coefficients[2, 1], 3)
}

reg.stats <- as.data.frame(reg.stats)
colnames(reg.stats) <- c('Model', 'R2', 'p', 'beta')

# --- Compare length of followup to model accuracy -------

# followup <- last.mmse %>%
#   select(ID, DELTA)
# colnames(followup) <- c('RID', 'DELTA')
# 
# predictions.with.followup <- left_join(predictions, followup, by='RID')
# 
# for (i in 1:length(model_keys)) {
#   x <- predictions.with.followup$DELTA
#   y <- predictions[[model_keys[i]]]
#   
#   pdata <- data.frame(x = x, y=y)
#   p <- ggplot(pdata, aes(x=x, y=y)) +
#     geom_point() +
#     geom_smooth(method='lm')
#   print(p)
#   
#   m <- lm(y ~ x, data = pdata)
#   print(summary(m))
# }