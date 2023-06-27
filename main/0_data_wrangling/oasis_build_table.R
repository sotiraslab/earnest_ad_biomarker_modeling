
# === Imports ============

sh <- suppressPackageStartupMessages

sh(library(ggplot2))
sh(library(stringr))
sh(library(this.path))
sh(library(tidyverse))

# === Set Working Directory ============

setwd(this.dir())

# === Files needed ===========

PATH.FTP <- '../../data/rawdata/oasis_flortaucipir.csv'
PATH.CENTILOID <- '../../data/derivatives/oasis_manual_centiloid.csv'
PATH.CLINICAL <- '../../data/rawdata/OASIS3_data_files/scans/UDSb4-Form_B4__Global_Staging__CDR__Standard_and_Supplemental/resources/csv/files/OASIS3_UDSb4_cdr.csv'
PATH.DEMO <- '../../data/rawdata/OASIS3_data_files/scans/demo-demographics/resources/csv/files/OASIS3_demographics.csv'
PATH.NEUROPSYCH <- '../../data/rawdata/OASIS3_data_files/scans/pychometrics-Form_C1__Cognitive_Assessments/resources/csv/files/OASIS3_UDSc1_cognitive_assessments.csv'
PATH.PACC.SCRIPT <- '../../scripts/pacc.R'

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
         PUPLongID = PUP_PUPTIMECOURSEDATA.ID) %>%
  mutate(Subject = str_extract(PUPLongID, 'OAS\\d+'),
         SessionAmyloid = adrc_session_to_number(PUPLongID),
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

# === Add regional tau =========

tau.rois <- read.csv(PATH.FTP)

# # === Select amyloid positive df ==========
# 
# # corpus callosum & cerebllum are omitted to match ADNI
# pat.inc <- "PET_fSUVR_TOT_CTX_.*"
# pat.exc <- "CRPCLM|CBLL"
# cols <- colnames(tau.rois)[grepl(pat.inc, colnames(tau.rois), perl = T) &
#                            ! grepl(pat.exc, colnames(tau.rois), perl = T)]
# 
# # ROIS from ADRC have different names than ADNI/ggseg
# # but are ordered the same
# adni.rois <- read.csv('../../derivatives/adni/nmf_regions.csv')
# converter <-  data.frame(ADRC=cols, ADNI=adni.rois$Feature)
# 
# tau.rois$Subject <- str_extract(tau.rois$PUP_PUPTIMECOURSEDATA.ID, 'OAS\\d+')
# 
# # add total cortical mean for a course tau index
# joiner <- tau.rois[, c("Subject", cols, 'PET_fSUVR_TOT_CORTMEAN')]
# colnames(joiner) <- c("Subject", converter$ADNI, 'TotalCtxTauMean')
# 
# df <- left_join(df, joiner, by='Subject')

# # === Add groups =============
# 
# df$Group <- NA
# 
# df$Group <- ifelse(is.na(df$Group) & df$Subject %in% df.amyloidpos$Subject,
#                    'TrainingSet',
#                    NA)
# 
# cdr0 <- (df$CDR == 0 & ! is.na(df$CDR))
# amyneg <- (! is.na(df$AmyloidPositive) & df$AmyloidPositive == F)
# both <- cdr0 & amyneg
# df$Group <- ifelse(is.na(df$Group) & both,
#                    'ControlSet',
#                    df$Group)
# 
# df$Group <- ifelse(is.na(df$Group),
#                    'Other',
#                    df$Group)
# 
# # === Calculate PACC ======
# 
# source(PATH.PACC.SCRIPT)
# 
# df$PACC.Original <- compute.pacc(df,
#                                  pacc.columns = c('srtfree', 'MEMUNITS', 'digsym', 'MMSE'),
#                                  cn.mask <- df$Group == 'ControlSet',
#                                  higher.better = c(T, T, T, T))
# 
# # === Generate subject list =========
# 
# # this code creates the subject-ids file for OASIS3
# # NOTE: overwriting the subject lists in the subject_ids folder
# #       may result in a different set of subjects being analyzed!
# 
# # sub.ids <- select(df, Subject, Session, Group)
# # write.csv(sub.ids, 'oasis_subjects.csv', row.names = F)
# 
# 
# # === Apply Subject list =========
# 
# original.subs <- read.csv(PATH.OASIS.SUBS)
# current.subs <- select(df, -Group)
# df.original.subs <- left_join(original.subs, current.subs, by=c('Subject', 'Session'))
# 
# # === Save ==========
# 
# write.csv(df.original.subs, '../../derivatives/oasis3/main_data.csv', na='', quote=F, row.names = F)
  