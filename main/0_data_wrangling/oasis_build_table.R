
# === Imports ============

sh <- suppressPackageStartupMessages

sh(library(ggplot2))
sh(library(stringr))
sh(library(this.path))
sh(library(tidyverse))

# === Set Working Directory ============

setwd(this.dir())

# === Files needed ===========


PATH.CENTILOID <- '../../data/derivatives/oasis_manual_centiloid.csv'
PATH.CLINICAL <- '../../data/rawdata/OASIS3_data_files/scans/UDSb4-Form_B4__Global_Staging__CDR__Standard_and_Supplemental/resources/csv/files/OASIS3_UDSb4_cdr.csv'
PATH.DEMO <- '../../data/rawdata/OASIS3_data_files/scans/demo-demographics/resources/csv/files/OASIS3_demographics.csv'
PATH.NEUROPSYCH <- '../../data/rawdata/OASIS3_data_files/scans/pychometrics-Form_C1__Cognitive_Assessments/resources/csv/files/OASIS3_UDSc1_cognitive_assessments.csv'
PATH.PACC.SCRIPT <- '../../scripts/pacc.R'

PATH.AV45 <- '../../data/rawdata/oasis_amyloid.csv'
PATH.AV45.ROIS <- '../../data/derivatives/av45_regions.csv'
PATH.FTP <- '../../data/rawdata/oasis_flortaucipir.csv'
PATH.FTP.ROIS <- '../../data/derivatives/ftp_regions.csv'
PATH.GM<- '../../data/rawdata/oasis_freesurfer.csv'
PATH.GM.ROIS <- '../../data/derivatives/gm_regions.csv'

# === Set variables ======

THRESHOLD.IMAGING.DAYS = 365
THRESHOLD.COGNITIVE.DAYS = 365

# === Load tau ROI data =============

# This contains all the ROI SUVR values extracted from
# PUP.  Both PVC and non-PVC.  This at least gives
# the very maximum number of subjects that can be used
# for NMF (434).  But need to filter out cases that are
# amyloid negative.

big.df <- read.csv(PATH.FTP)

# extract only the subjects
adrc_session_to_number <- function(col) {
  extracted <- str_extract(col, 'd\\d+')
  no.leading.d <- substr(extracted, 2, nchar(extracted))
  number <- as.numeric(no.leading.d)
  
  return (number)
}

df <- big.df[, c("PUP_PUPTIMECOURSEDATA.ID", 'FSId', 'MRId'), drop=0]
colnames(df) <- c('TauID', 'FSID', 'MRIID')
df$Subject <- str_extract(df$TauID, 'OAS\\d+')
df$SessionTau <- adrc_session_to_number(df$TauID)

# === Add Amyloid ===========

# based on the amyloid cortical values
# see check_centiloid_conversion.R for a little more explanation

amyloid.df <- read.csv(PATH.CENTILOID)

amyloid.df.proc <- amyloid.df %>%
  select(PUP_PUPTIMECOURSEDATA.ID,
         tracer,
         ManualCentiloid) %>%
  rename(AmyloidTracer = tracer,
         Centiloid = ManualCentiloid,
         AmyloidID = PUP_PUPTIMECOURSEDATA.ID) %>%
  mutate(Subject = str_extract(AmyloidID, 'OAS\\d+'),
         SessionAmyloid = adrc_session_to_number(AmyloidID),
         AmyloidPositive = ifelse(AmyloidTracer == 'AV45',
                                  Centiloid >= 20.6,
                                  Centiloid >= 16.4)) %>%
  filter(AmyloidTracer == 'AV45')

link.av45 <- left_join(df, amyloid.df.proc, by='Subject') %>%
  mutate(DiffTauAmyloid = SessionTau - SessionAmyloid) %>%
  group_by(Subject) %>%
  slice_min(abs(DiffTauAmyloid), with_ties = F) %>%
  filter(abs(DiffTauAmyloid) < THRESHOLD.IMAGING.DAYS)

# df.merge <- left_join(df, amyloid.df.proc, by='Subject') %>%
#   mutate(TauAmyloidDiff = Session - SessionAmyloid) %>%
#   group_by(Subject) %>%
#   summarise(Indexer=amyloid.selector(Indexer, TauAmyloidDiff, AmyloidPositive)) %>%
#   ungroup() %>%
#   left_join(select(amyloid.df.proc, -Subject), by='Indexer')
# 
# # assign status
# df <- left_join(df, df.merge, by='Subject') %>%
#   mutate(TauAmyloidDiff = Session - SessionAmyloid)

df <- link.av45 %>%
  mutate(MeanImagingDate = (SessionTau + SessionAmyloid) / 2)

# === Add CDR ==========

clinical.df <- read.csv(PATH.CLINICAL)

clinical.df.proc <- clinical.df %>%
  select(OASIS_session_label, CDRTOT, CDRSUM) %>%
  mutate(Subject=str_extract(OASIS_session_label, 'OAS\\d+'),
         SessionCDR=adrc_session_to_number(OASIS_session_label)) %>%
  rename(CDR=CDRTOT, CDRSumBoxes=CDRSUM) %>%
  select(-OASIS_session_label)

df <- left_join(df, clinical.df.proc, by='Subject') %>%
  mutate(DiffMeanImagingCDR = abs(MeanImagingDate - SessionCDR)) %>%
  group_by(Subject) %>%
  slice_min(DiffMeanImagingCDR, with_ties = F) %>%
  mutate(CDR = ifelse(DiffMeanImagingCDR > THRESHOLD.COGNITIVE.DAYS, NA, CDR),
         CDRSumBoxes = ifelse(DiffMeanImagingCDR > THRESHOLD.COGNITIVE.DAYS, NA, CDRSumBoxes),
         CDRBinned=cut(CDR, breaks=c(0, .5, 1, Inf), right=F)) %>%
  drop_na(CDR)

levels(df$CDRBinned) <- c("0.0", "0.5", "1.0+")
df$CDRBinned <- as.character(df$CDRBinned)

df$Dementia <- ifelse(df$CDR >= 0.5 & ! is.na(df$CDR), 
                      'Yes',
                      'No')
df[is.na(df$CDR), 'Dementia'] <- 'Unknown'
df$Control <- ifelse(! df$AmyloidPositive & df$Dementia == 'No', 1, 0)

# === Add Demographics ==========

demo <- read.csv(PATH.DEMO)

demo.proc <- demo %>%
  select(OASISID, AgeatEntry, AgeatDeath, GENDER, APOE) %>%
  rename(Subject=OASISID, Sex=GENDER) %>%
  mutate(HasE4 = grepl('4', APOE),
         Sex=ifelse(Sex == 1, 'Male', 'Female'))

df <- left_join(df, demo.proc, by='Subject')
df$Age <- df$AgeatEntry + (df$MeanImagingDate / 365.25)

# === Add MMSE ==========

mmse <- read.csv(PATH.CLINICAL)

mmse <- mmse %>%
  select(OASISID, MMSE, days_to_visit) %>%
  rename(SessionMMSE=days_to_visit,
         Subject=OASISID) %>%
  select(Subject, MMSE, SessionMMSE)

df <- left_join(df, mmse, by='Subject') %>%
  mutate(DiffMeanImagingMMSE = abs(MeanImagingDate - SessionMMSE)) %>%
  group_by(Subject) %>%
  slice_min(DiffMeanImagingMMSE, with_ties = F) %>%
  mutate(MMSE = ifelse(DiffMeanImagingMMSE > THRESHOLD.COGNITIVE.DAYS, NA, MMSE))

# === Add neuropsych ==========

nps <- read.csv(PATH.NEUROPSYCH)

nps <- nps %>%
  select(OASISID, days_to_visit, srtfree, MEMUNITS, digsym, ANIMALS, tmb) %>%
  rename(SessionNeuroPsych = days_to_visit,
         Subject=OASISID)

df <- left_join(df, nps, by='Subject') %>%
  mutate(DiffMeanImagingNPS = abs(MeanImagingDate - SessionNeuroPsych)) %>%
  group_by(Subject) %>%
  slice_min(DiffMeanImagingNPS, with_ties = F)

bad.nps <- (df$DiffMeanImagingNPS > THRESHOLD.COGNITIVE.DAYS) | (is.na(df$DiffMeanImagingNPS))
df[bad.nps, c('srtfree', 'MEMUNITS', 'digsym', 'ANIMALS', 'tmb')] <- NA

# === Calculate PACC ======

source(PATH.PACC.SCRIPT)

df$PACC <- compute.pacc(df,
                        pacc.columns = c('srtfree', 'MEMUNITS', 'digsym', 'MMSE'),
                        cn.mask = df$Control == 1,
                        higher.better = c(T, T, T, T))

# === Configure subcortical regions

# MAKE SURE THIS MATCHES THE ORDER AS ADNI!
# SUBCORTICAL = c('amygdala',
#                 'caudate',
#                 'hippocampus',
#                 'pallidum',
#                 'putamen',
#                 'thalamus',
#                 'ventraldc')

SUBCORTICAL = c('AMYGDALA',
                'CAUD',
                'HIPPOCAMPUS',
                'PALLIDUM',
                'PUTAMEN',
                'THALAMUS',
                'VENTRALDC')
SUBCORTICAL_PAT <- paste(SUBCORTICAL, collapse='|')

# === Add regional AV45 =========

rois <- read.csv(PATH.AV45)

# corpus callosum & cerebllum are omitted to match ADNI
cort.cols <- colnames(rois)[grepl("PET_fSUVR_(L|R)_CTX_.*", colnames(rois), perl = T) &
                              ! grepl("CRPCLM|CBLL", colnames(rois), perl = T)]
subcort.cols <- colnames(rois)[grepl(SUBCORTICAL_PAT, colnames(rois), perl = T) &
                                 grepl('PET_fSUVR_(L|R)', colnames(rois), perl = T) &
                                 ! grepl('CTX|WM', colnames(rois), perl = T)]
cols <- c(cort.cols, subcort.cols)

# ROIS from OASIS & ADNI have different names but should be ordered the same
# good point to verify this!
adni.rois <- read.csv(PATH.AV45.ROIS)
converter <-  data.frame(OASIS=cols, ADNI=adni.rois$Region)

# merge
rois$AmyloidID <- rois$PUP_PUPTIMECOURSEDATA.ID
joiner <- rois[, c("AmyloidID", cols)]
colnames(joiner) <- c("AmyloidID", converter$ADNI)
df <- left_join(df, joiner, by='AmyloidID')

# === Add regional tau =========

rois <- read.csv(PATH.FTP)

# corpus callosum & cerebllum are omitted to match ADNI
cort.cols <- colnames(rois)[grepl("PET_fSUVR_(L|R)_CTX_.*", colnames(rois), perl = T) &
                            ! grepl("CRPCLM|CBLL", colnames(rois), perl = T)]
subcort.cols <- colnames(rois)[grepl(SUBCORTICAL_PAT, colnames(rois), perl = T) &
                               grepl('PET_fSUVR_(L|R)', colnames(rois), perl = T) &
                               ! grepl('CTX|WM', colnames(rois), perl = T)]
cols <- c(cort.cols, subcort.cols)

# ROIS from OASIS & ADNI have different names but should be ordered the same
# good point to verify this!
adni.rois <- read.csv(PATH.FTP.ROIS)
converter <-  data.frame(OASIS=cols, ADNI=adni.rois$Region)

# merge
rois$TauID <- rois$PUP_PUPTIMECOURSEDATA.ID
joiner <- rois[, c("TauID", cols)]
colnames(joiner) <- c("TauID", converter$ADNI)
df <- left_join(df, joiner, by='TauID')

# === Add regional GM volume =========

fs <- read.csv(PATH.GM)

rois <- fs %>%
  select(matches('^(lh|rh)_.*_volume$'),
         matches(SUBCORTICAL_PAT, ignore.case = T) &
           contains('volume') &
           ! contains('WM') &
           ! contains('TOTAL'))

# check name conversion
adni.rois <- read.csv(PATH.GM.ROIS)
converter <-  data.frame(OASIS=colnames(rois), ADNI=adni.rois$Region)
colnames(rois) <- converter$ADNI

merger <- rois %>%
  mutate(FSID = fs$FS_FSDATA.ID,
         ICV = fs$IntraCranialVol)

df <- left_join(df, merger, by = 'FSID') %>%
  mutate(across(ends_with('_VOLUME'), function (x) (x * 1000 / ICV)))

# === remove NAs =========

all.roi.cols <- colnames(df)[str_detect(colnames(df), '_SUVR$|_VOLUME$')]
na.cols <- c('Age', 'Sex', 'PACC', 'HasE4', 'CDRBinned', all.roi.cols)

df.withna <- df
df <- df %>%
  drop_na(all_of(na.cols))

# === Add engineered features =========

# # these should closely match ADNI
# # so that models estimated in ADNI can be applied in OASIS
# # without renaming the features
# 
# # Note that ADNI mostly uses volume-weighted uptakes
# # for PET.  This is apparently calulated using
# # unilateral ROIs, which is repeated here.
# 
# base <- df[, c('AmyloidID', 'TauID', 'FSID'), drop=F]
# bilateral.cols <- c(gsub('TOT', 'L', cort.cols), gsub('TOT', 'R', cort.cols),
#                     gsub('TOT', 'L', subcort.cols), gsub('TOT', 'R', subcort.cols))
# 
# # unilateral av45
# rois.av45 <- read.csv(PATH.AV45)
# rois.av45 <- rois.av45 %>%
#   mutate(AmyloidID = PUP_PUPTIMECOURSEDATA.ID) %>%
#   select(AmyloidID,
#          matches('PET_fSUVR_(L|R)_CTX_.*') & ! matches('CBLL') & ! matches('CRPCLM'),
#          matches(SUBCORTICAL_PAT) & ! contains('_TOT_') & ! contains('WM') & ! contains('CTX'))
# rois.av45 <- left_join(base, rois.av45, by='AmyloidID') %>%
#   select(-AmyloidID, -TauID, -FSID)
# check.av45 <- data.frame(OASIS=colnames(rois.av45), bilateral.cols)
# colnames(rois.av45) <- bilateral.cols
# 
# # unilateral tau
# rois.tau <- read.csv(PATH.FTP)
# rois.tau <- rois.tau %>%
#   mutate(TauID = PUP_PUPTIMECOURSEDATA.ID) %>%
#   select(TauID,
#          matches('PET_fSUVR_(L|R)_CTX_.*') & ! matches('CBLL') & ! matches('CRPCLM'),
#          matches(SUBCORTICAL_PAT) & ! contains('_TOT_') & ! contains('WM') & ! contains('CTX'))
# rois.tau <- left_join(base, rois.tau, by='TauID') %>%
#   select(-AmyloidID, -TauID, -FSID)
# check.tau <- data.frame(OASIS=colnames(rois.tau), bilateral.cols)
# colnames(rois.tau) <- bilateral.cols
# 
# # unilateral GM
# rois.gm <- read.csv(PATH.GM)
# rois.gm <- rois.gm %>%
#   mutate(FSID = FS_FSDATA.ID) %>%
#   select(FSID,
#          contains('lh_') & contains('volume') & ! contains('WM') & ! contains('TOTAL'),
#          contains('rh_') & contains('volume') & ! contains('WM') & ! contains('TOTAL'),
#          matches(SUBCORTICAL_PAT) & contains('volume') & ! contains('WM') & ! contains('TOTAL'))
# rois.gm <- left_join(base, rois.gm, by='FSID') %>%
#   select(-AmyloidID, -TauID, -FSID)
# check.gm <- data.frame(OASIS=colnames(rois.gm), bilateral.cols)
# colnames(rois.gm) <- bilateral.cols
# 
# # create volume weighting function
# volume.weighted.mean <- function(pet.rois, volumes, columns) {
#   pet.rois <- rois.tau
#   volumes <- rois.gm
#   columns <- c("PET_fSUVR_L_CTX_ENTORHINAL", "PET_fSUVR_R_CTX_ENTORHINAL")
#   
#   pet <- pet.rois[, columns]
#   volumes <- volumes[, columns]
#   volumes.norm <- volumes / rowSums(volumes)
#   pet.norm <- pet * volumes.norm
#   result <- rowSums(pet.norm)
#   
#   return(result)
# }
# 
# comp.amyloid.regs <- c('CAUDMIDFRN',
#                        'LATORBFRN',
#                        'MEDORBFRN',
#                        'PARSOPCLRS',
#                        'PARSORBLS',
#                        'PARSTRNGLS',
#                        'ROSMIDFRN',
#                        'SUPERFRN',
#                        'FRNPOLE',
#                        'CAUDANTCNG',
#                        'ISTHMUSCNG',
#                        'POSTCNG',
#                        'ROSANTCNG',
#                        'INFERPRTL',
#                        'PRECUNEUS',
#                        'SUPERPRTL',
#                        'SUPRAMRGNL',
#                        'INFERTMP',
#                        'MIDTMP',
#                        'SUPERTMP')
# comp.amyloid.cols <- bilateral.cols[str_detect(bilateral.cols, paste(comp.amyloid.regs, collapse='|'))]
# 
# braak1.regs <- c('ENTORHINAL')
# braak1.cols <- bilateral.cols[str_detect(bilateral.cols, paste(braak1.regs, collapse='|'))]
# 
# braak34.regs <- c('PARAHPCMPL',
#                   'FUSIFORM',
#                   'LINGUAL',
#                   'AMYGDALA',
#                   'MIDTMP',
#                   'CAUDANTCNG',
#                   'ROSANTCNG',
#                   'POSTCNG',
#                   'ISTHMUSCNG',
#                   'INSULA',
#                   'INFERTMP',
#                   'TMPPOLE')
# braak34.cols <- bilateral.cols[str_detect(bilateral.cols, paste(braak34.regs, collapse='|'))]
# 
# braak56.regs <- c('SUPERFRN',
#                   'LATORBFRN',
#                   'MEDORBFRN',
#                   'FRNPOLE',
#                   'CADMIDFRN',
#                   'ROSMIDFRN',
#                   'PARSOPCLRS',
#                   'PARSORBLS',
#                   'PARSTRNGLS',
#                   'LATOCC',
#                   'SUPRAMRGNL',
#                   'INFERPRTL',
#                   'SUPERTMP',
#                   'SUPERPRTL',
#                   'PRECUNEUS',
#                   'SSTSBANK',
#                   'TRANSTMP',
#                   'PERICLCRN',
#                   'POSTCNTRL',
#                   'CUNEUS',
#                   'PRECENTRL',
#                   'PARACNTRL')
# braak56.cols <- bilateral.cols[str_detect(bilateral.cols, paste(braak56.regs, collapse='|'))]
# 
# mtt.regs <- c('AMYGDALA',
#               'ENTORHINAL',
#               'FUSIFORM',
#               'INFERTMP',
#               'MIDTMP')
# mtt.cols <- bilateral.cols[str_detect(bilateral.cols, paste(mtt.regs, collapse='|'))]
# 
# # Amyloid
# #  - Centiloid
# #  - SUMMARYSUVR_WHOLECEREBNORM
# #
# # Tau
# #  - META_TEMPORAL_SUVR
# #  - BRAAK1_SUVR
# #  - BRAAK34_SUVR
# #  - BRAAK56_SUVR
# # 
# # GM
# #  - HIPPOCAMPUS_VOLUME
# #  - META_TEMPORAL_VOLUME
# 
# 
# # apply!
# df$SUMMARYSUVR_WHOLECEREBNORM <- volume.weighted.mean(rois.av45, rois.gm, comp.amyloid.cols)
# df$META_TEMPORAL_SUVR <- volume.weighted.mean(rois.tau, rois.gm, mtt.cols)
# df$BRAAK1_SUVR <- volume.weighted.mean(rois.tau, rois.gm, braak1.cols)
# df$BRAAK34_SUVR <- volume.weighted.mean(rois.tau, rois.gm, braak34.cols)
# df$BRAAK56_SUVR <- volume.weighted.mean(rois.tau, rois.gm, braak56.cols)
# 
# df$META_TEMPORAL_VOLUME <- (rowSums(rois.gm[, mtt.cols]) * 1000) / df$ICV

# === Save ==========

write.csv(df, '../../data/derivatives/oasis_base_table.csv', na='', quote=F, row.names = F)
  