
library(anticlust)
library(argparse)
library(this.path)

# Inputs
inpath <- '/Users/earnestt1234/Desktop/training.csv'
outpath <- '/Users/earnestt1234/Desktop/splits.csv'
K <- 2
N <- 15
id.columns <- c('Subject', 'Session')
continuous.columns <- c('Age', 'SummarySUVRAmyloid', 'SummarySUVRTau')
categorical.columns <- c('Dataset', 'CDRBinned')
verbose <- TRUE
seed <- 42

# ---- MAIN -----

set.seed(seed)

df <- read.csv(path)

splits <- matrix(data=NA, nrow=nrow(df), ncol=N)
continuous.input <- df[, continuous.columns]
categorical.input <- df[, categorical.columns]

for (i in 1:N) {
  v <- anticlustering(
    x = continuous.input,
    k = K,
    categories = categorical.input
    )
  splits[, i] <- v
}

splits <- as.data.frame(splits)

if (! is.null(id.columns)) {
  ids <- df[, id.columns]
  splits <- cbind(ids, splits)
}

write.table(splits, outpath, row.names = F, col.names = F, sep = ',')