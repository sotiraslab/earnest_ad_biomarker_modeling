
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
         CDRSumBoxes = ifelse(DiffMeanImagingCDR > THRESHOLD.COGNITIVE.DAYS, NA, CDRSumBoxes)) %>%
  drop_na(CDR)

df$Dementia <- ifelse(df$CDR >= 0.5 & ! is.na(df$CDR), 
                      'Yes',
                      'No')
df[is.na(df$CDR), 'Dementia'] <- 'Unknown'
df$Control <- ifelse(! df$AmyloidPositive & df$Dementia == 'No', 1, 0)

# === Add Demographics ==========

demo <- read.csv(PATH.DEMO)

demo.proc <- demo %>%
  select(OASISID, AgeatEntry, AgeatDeath, GENDER, APOE) %>%
  rename(Subject=OASISID) %>%
  mutate(HasE4 = grepl('4', APOE))

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

df$PACC.Original <- compute.pacc(df,
                                 pacc.columns = c('srtfree', 'MEMUNITS', 'digsym', 'MMSE'),
                                 cn.mask <- df$Control,
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
cort.cols <- colnames(rois)[grepl("PET_fSUVR_TOT_CTX_.*", colnames(rois), perl = T) &
                              ! grepl("CRPCLM|CBLL", colnames(rois), perl = T)]
subcort.cols <- colnames(rois)[grepl(SUBCORTICAL_PAT, colnames(rois), perl = T) &
                                 grepl('PET_fSUVR_TOT', colnames(rois), perl = T) &
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
cort.cols <- colnames(rois)[grepl("PET_fSUVR_TOT_CTX_.*", colnames(rois), perl = T) &
                            ! grepl("CRPCLM|CBLL", colnames(rois), perl = T)]
subcort.cols <- colnames(rois)[grepl(SUBCORTICAL_PAT, colnames(rois), perl = T) &
                               grepl('PET_fSUVR_TOT', colnames(rois), perl = T) &
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
  select(Subject,
         matches('^(lh|rh)_.*_volume$'),
         matches(SUBCORTICAL_PAT, ignore.case = T) &
           contains('volume') &
           ! contains('WM') &
           ! contains('TOTAL'))

lh <- select(rois, contains('lh_'), contains('Left'))
rh <- select(rois, contains('rh_'), contains('Right'))
rois.bilateral <- (lh + rh)

# check name conversion
adni.rois <- read.csv(PATH.GM.ROIS)
converter <-  data.frame(OASIS=colnames(rois.bilateral), ADNI=adni.rois$Region)
colnames(rois.bilateral) <- converter$ADNI

merger <- rois.bilateral %>%
  mutate(FSID = fs$FS_FSDATA.ID,
         ICV = fs$IntraCranialVol)

df <- left_join(df, merger, by = 'FSID') %>%
  mutate(across(ends_with('_VOLUME'), function (x) (x * 1000 / ICV)))

# === Save ==========

write.csv(df, '../../data/derivatives/oasis_base_table.csv', na='', quote=F, row.names = F)
  