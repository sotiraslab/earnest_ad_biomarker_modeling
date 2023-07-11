# === imports ======

sh <- suppressPackageStartupMessages

sh(library(lubridate))
sh(library(stringr))
sh(library(this.path))
sh(library(tidyverse))

# === Set working directory ======

setwd(this.dir())

# === Set paths ======

PATH.ADNI <- '../../data/derivatives/adni_harmonized.csv'
PATH.OASIS <- '../../data/derivatives/oasis_harmonized.csv'

# === read data =======

adni <- read.csv(PATH.ADNI)
oasis <- read.csv(PATH.OASIS)

av45.cols <- colnames(adni)[str_detect(colnames(adni), '^AV45')]
ftp.cols <- colnames(adni)[str_detect(colnames(adni), '^FTP')]
gm.cols <- colnames(adni)[str_detect(colnames(adni), '_VOLUME')]

adni.av45 <- adni[, av45.cols]
adni.ftp <- adni[, ftp.cols]
adni.gm <- adni[, gm.cols]

oasis.av45 <- oasis[, av45.cols]
oasis.ftp <- oasis[, ftp.cols]
oasis.gm <- oasis[, gm.cols]

# === ROI definitions ======

comp.amyloid.regs <- c('CAUDALMIDDLEFRONTAL',
                       'LATERALORBITOFRONTAL',
                       'MEDIALORBITOFRONTAL',
                       'PARSOPERCULARIS',
                       'PARSORBITALIS',
                       'PARSTRIANGULARIS',
                       'ROSTRALMIDDLEFRONTAL',
                       'SUPERIORFRONTAL',
                       'FRONTALPOLE',
                       'CAUDALANTERIORCINGULATE',
                       'ISTHMUSCINGULATE',
                       'POSTERIORCINGULATE',
                       'ROSTRALANTERIORCINGULATE',
                       'INFERIORPARIETAL',
                       'PRECUNEUS',
                       'SUPERIORPARIETAL',
                       'SUPRAMARGINAL',
                       'INFERIORTEMPORAL',
                       'MIDDLETEMPORAL',
                       'SUPERIORTEMPORAL')

braak1.regs <- c('ENTORHINAL')

braak34.regs <- c('PARAHIPPOCAMPAL',
                  'FUSIFORM',
                  'LINGUAL',
                  'AMYGDALA',
                  'MIDDLETEMPORAL',
                  'CAUDALANTERIORCINGULATE',
                  'ROSTRALANTERIORCINGULATE',
                  'POSTERIORCINGULATE',
                  'ISTHMUSCINGULATE',
                  'INSULA',
                  'INFERIORTEMPORAL',
                  'TEMPORALPOLE')

braak56.regs <- c('SUPERIORFRONTAL',
                  'LATERALORBITOFRONTAL',
                  'MEDIALORBITOFRONTAL',
                  'FRONTALPOLE',
                  'CAUDALMIDDLEFRONTAL',
                  'ROSTRALMIDDLEFRONTAL',
                  'PARSOPERCULARIS',
                  'PARSORBITALIS',
                  'PARSTRIANGULARIS',
                  'LATERALOCCIPITAL',
                  'SUPRAMARGINAL',
                  'INFERIORPARIETAL',
                  'SUPERIORTEMPORAL',
                  'SUPERIORPARIETAL',
                  'PRECUNEUS',
                  'BANKSSTS',
                  'TRANSVERSETEMPORAL',
                  'PERICALCARINE',
                  'POSTCENTRAL',
                  'CUNEUS',
                  'PRECENTRAL',
                  'PARACENTRAL')

mtt.regs <- c('AMYGDALA',
              'ENTORHINAL',
              'FUSIFORM',
              'INFERIORTEMPORAL',
              'MIDDLETEMPORAL')

# ==== compute ======

volume.weighted.mean <- function(pet.data, vol.data, search.columns) {
  
  pattern <- paste(search.columns, collapse='|')
  pet.cols <- colnames(pet.data)[str_detect(colnames(pet.data), pattern)]
  vol.cols <- colnames(vol.data)[str_detect(colnames(vol.data), pattern)]
  
  volumes.norm <- vol.data / rowSums(vol.data)
  pet.norm <- pet.data * volumes.norm
  result <- rowSums(pet.norm)

  return(result)
}

adni$AMYLOID_COMPOSITE <- volume.weighted.mean(adni.av45, adni.gm, comp.amyloid.regs)
adni$META_TEMPORAL_SUVR <- volume.weighted.mean(adni.ftp, adni.gm, mtt.regs)
adni$BRAAK1_SUVR <- volume.weighted.mean(adni.ftp, adni.gm, braak1.regs)
adni$BRAAK34_SUVR <- volume.weighted.mean(adni.ftp, adni.gm, braak34.regs)
adni$BRAAK56_SUVR <- volume.weighted.mean(adni.ftp, adni.gm, braak56.regs)

oasis$AMYLOID_COMPOSITE <- volume.weighted.mean(oasis.av45, oasis.gm, comp.amyloid.regs)
oasis$META_TEMPORAL_SUVR <- volume.weighted.mean(oasis.ftp, oasis.gm, mtt.regs)
oasis$BRAAK1_SUVR <- volume.weighted.mean(oasis.ftp, oasis.gm, braak1.regs)
oasis$BRAAK34_SUVR <- volume.weighted.mean(oasis.ftp, oasis.gm, braak34.regs)
oasis$BRAAK56_SUVR <- volume.weighted.mean(oasis.ftp, oasis.gm, braak56.regs)

# mtt volume
mtt.cols <- gm.cols[str_detect(gm.cols, paste(mtt.regs, collapse='|'))]
adni$META_TEMPORAL_VOLUME <- rowSums(adni.gm[, mtt.cols]) 
adni$HIPPOCAMPUS_VOLUME <- rowSums(adni.gm[, c('LEFT_HIPPOCAMPUS_VOLUME', 'RIGHT_HIPPOCAMPUS_VOLUME')])

oasis$META_TEMPORAL_VOLUME <- rowSums(oasis.gm[, mtt.cols])
oasis$HIPPOCAMPUS_VOLUME <- rowSums(oasis.gm[, c('LEFT_HIPPOCAMPUS_VOLUME', 'RIGHT_HIPPOCAMPUS_VOLUME')])

# === save =======

write.csv(adni, '../../data/derivatives/adni_harmonized_augmented.csv', row.names = F, quote = F, na = '')
write.csv(oasis, '../../data/derivatives/oasis_harmonized_augmented.csv', row.names = F, quote = F, na = '')
