# === imports ======

sh <- suppressPackageStartupMessages

sh(library(ADNIMERGE))
sh(library(lubridate))
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
