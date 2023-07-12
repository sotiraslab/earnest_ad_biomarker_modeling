
# === Imports ============

sh <- suppressPackageStartupMessages

sh(library(ggplot2))
sh(library(lme4))
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

# === Calculate PACC ======

subs <- df %>%
  select(Subject)

# 1. MMSE
mmse <- read.csv(PATH.CLINICAL)

mmse <- mmse %>%
  select(OASISID, MMSE, days_to_visit) %>%
  rename(SessionNPS=days_to_visit,
         Subject=OASISID) %>%
  select(Subject, SessionNPS, MMSE)

mmse <- left_join(subs, mmse, by='Subject')

# 2. Neuropsych
nps <- read.csv(PATH.NEUROPSYCH)

nps <- nps %>%
  select(OASISID, days_to_visit, srtfree, MEMUNITS, digsym, ANIMALS, tmb) %>%
  rename(SessionNPS = days_to_visit,
         Subject=OASISID)

nps <- left_join(subs, nps, by='Subject')

# 3. Join all
pacc.df <- df %>%
  select(Subject) %>%
  full_join(mmse, by='Subject') %>%
  full_join(nps, by=c('Subject', 'SessionNPS')) %>%
  arrange(Subject, SessionNPS) %>%
  drop_na(SessionNPS)

# Group assessments which occur at nearby dates
LINK.THR <- 60 # days

pacc.df <- pacc.df %>%
  group_by(Subject) %>%
  mutate(NPSSessionDiff = c(0, diff(SessionNPS))) %>%
  ungroup()

pacc.df$NPSSessionGroup <- cumsum(abs(pacc.df$NPSSessionDiff) > LINK.THR)

select.first <- function(col) {
  return(first(col[! is.na(col)]))
}

pacc.df.group <- pacc.df %>%
  group_by(Subject, NPSSessionGroup) %>%
  summarise(
    SessionPACC = mean(SessionNPS),
    SessionPACCMin = min(SessionNPS),
    SessionPACCMax = max(SessionNPS),
    SessionPACCDiff = SessionPACCMax - SessionPACCMin,
    MMSE = as.numeric(select.first(MMSE)),
    srtfree = as.numeric(select.first(srtfree)),
    MEMUNITS = as.numeric(select.first(MEMUNITS)),
    digsym = as.numeric(select.first(digsym)),
    ANIMALS = as.numeric(select.first(ANIMALS)),
    tmb = as.numeric(select.first(tmb)),
  ) %>%
  ungroup()

# merge into df
df <- left_join(df, pacc.df.group, by="Subject") %>%
  mutate(DiffImagingPACC = MeanImagingDate - SessionPACC) %>%
  group_by(TauID) %>%
  slice_min(abs(DiffImagingPACC), with_ties = F) %>%
  filter(abs(DiffImagingPACC) < THRESHOLD.COGNITIVE.DAYS) %>%
  ungroup() %>%
  arrange(Subject)

# compute
source(PATH.PACC.SCRIPT)

df$PACC <- compute.pacc(df,
                        pacc.columns = c('srtfree', 'MEMUNITS', 'digsym', 'MMSE'),
                        cn.mask = df$Dementia == 'No',
                        higher.better = c(T, T, T, T),
                        min.required = 2)

# === remove NAs =========

all.roi.cols <- colnames(df)[str_detect(colnames(df), '_SUVR$|_VOLUME$')]
na.cols <- c('Age', 'Sex', 'PACC', 'HasE4', 'CDRBinned', all.roi.cols)

df.withna <- df
df <- df %>%
  drop_na(all_of(na.cols))

# === Compute longitudinal change in PACC =======

# this will only be available for some

cn.data <- df %>%
  filter(Dementia == "No")

pacc.long <- df %>%
  select(Subject, SessionPACC, Age, CDRBinned) %>%
  rename(SessionPACC.BL=SessionPACC) %>%
  left_join(pacc.df.group, by="Subject") %>%
  group_by(Subject) %>%
  filter(SessionPACC >= SessionPACC.BL) %>%
  ungroup()

pacc.long$PACC <- compute.pacc(pacc.long,
                               pacc.columns = c('srtfree', 'MEMUNITS', 'digsym', 'MMSE'),
                               cn.data = cn.data,
                               higher.better = c(T, T, T, T),
                               min.required = 2)
pacc.long <- pacc.long %>%
  filter(! is.na(PACC)) %>%
  group_by(Subject) %>%
  filter(n() >= 2) %>%
  mutate(DeltaPACCSession = (SessionPACC - SessionPACC.BL) / 365.25,
         Long.Age = Age + DeltaPACCSession) %>% 
  ungroup()

m <- lmer(PACC ~ DeltaPACCSession + (1+DeltaPACCSession|Subject), data=pacc.long)
pacc.long$PACC.LMER.Predict <- predict(m, pacc.long)

ggplot(pacc.long, aes(x=Long.Age, y=PACC)) +
  geom_point(aes(color=CDRBinned), alpha = .7) + 
  geom_line(aes(y=PACC.LMER.Predict, group=Subject, color=CDRBinned), alpha= .7)

ggsave("oasis_longitudinal_pacc_model.png", width=8, height=6, units='in')

coefs <- coef(m)$Subject %>%
  select(DeltaPACCSession) %>%
  dplyr::rename(DeltaPACC=DeltaPACCSession) %>%
  rownames_to_column(var="Subject")

df <- left_join(df, coefs, by='Subject')

# === Configure subcortical regions =======

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

# === Save ==========

write.csv(df, '../../data/derivatives/oasis_base_table.csv', na='', quote=F, row.names = F)
  