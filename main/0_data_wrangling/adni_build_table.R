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

PATH.PTDEMOG.CSV <- '../../data/rawdata/PTDEMOG.csv'
PATH.ICV <- '../../data/derivatives/adni_icvs.csv'

PATH.PACC.SCRIPT <- '../../scripts/pacc.R'
PATH.EXAMDATE.SCRIPT <- '../../scripts/adni_examdate.R'

PATH.OUTPUT <- '../../data/derivatives/adni_base_table.csv'
PATH.DERIVATIVES = '../../data/derivatives'

source(PATH.PACC.SCRIPT)
source(PATH.EXAMDATE.SCRIPT)

# === Set variables ======

THRESHOLD.IMAGING.DAYS = 365
THRESHOLD.COGNITIVE.DAYS = 365

# === Find amyloid/tau overlap ======

scan.id <- function(RID, EXAMDATE) {
  return (paste(as.character(RID), gsub('-', '', EXAMDATE), sep='-'))
}

# merge on tau as that is less available
# separate merges for amyloid/FBB and concatenation
# just need to check that some individuals didn't get duplicate scans

tau <- ucberkeleyav1451 %>%
  mutate(DateTau = as_datetime(ymd(EXAMDATE)),
         TauID = scan.id(RID, EXAMDATE)) %>%
  select(RID, DateTau, TauID) %>%
  group_by(RID) %>%
  slice_min(DateTau, with_ties = F) %>%
  ungroup()

# Centiloid conversion for ADNI: 
# https://adni.loni.usc.edu/wp-content/themes/freshnews-dev-v2/documents/pet/ADNI%20Centiloids%20Final.pdf

av45 <- ucberkeleyav45 %>%
  mutate(DateAmyloid=as_datetime(ymd(EXAMDATE)),
         AmyloidID = scan.id(RID, EXAMDATE),
         AmyloidTracer = 'AV45',
         AmyloidPositive = as.numeric(SUMMARYSUVR_WHOLECEREBNORM_1.11CUTOFF),
         AmyloidID = scan.id(RID, EXAMDATE),
         Centiloid = as.numeric((196.9*SUMMARYSUVR_WHOLECEREBNORM) - 196.03)) %>%
  select(RID, DateAmyloid, AmyloidID, AmyloidTracer,
         AmyloidPositive, AmyloidID, Centiloid)

# fbb <- ucberkeleyfbb %>%
#   mutate(DateAmyloid=as_datetime(ymd(EXAMDATE)),
#          AmyloidID = scan.id(RID, EXAMDATE),
#          AmyloidTracer = 'FBB',
#          AmyloidPositive = as.numeric(SUMMARYSUVR_WHOLECEREBNORM_1.08CUTOFF),
#          AmyloidID = scan.id(RID, EXAMDATE),
#          Centiloid = as.numeric((159.08*SUMMARYSUVR_WHOLECEREBNORM) - 151.65)) %>%
#   select(RID, DateAmyloid, AmyloidID, AmyloidTracer, AmyloidPositive, AmyloidID, Centiloid)

link.av45 <- left_join(tau, av45, by='RID') %>%
  mutate(DiffTauAmyloid = as.numeric(difftime(DateTau, DateAmyloid, units = 'days'))) %>%
  group_by(TauID) %>%
  slice_min(abs(DiffTauAmyloid), with_ties = F) %>%
  filter(abs(DiffTauAmyloid) < THRESHOLD.IMAGING.DAYS)

# link.fbb <- left_join(tau, fbb, by='RID') %>%
#   mutate(DiffTauAmyloid = as.numeric(difftime(DateTau, DateAmyloid, units = 'days'))) %>%
#   group_by(TauID) %>%
#   slice_min(abs(DiffTauAmyloid), with_ties = F) %>%
#   filter(abs(DiffTauAmyloid) < THRESHOLD.IMAGING.DAYS)

df <- link.av45 %>%
  mutate(MeanImagingDate = as_date(DateTau - (as.difftime(DiffTauAmyloid, units = 'days') / 2)),
         RID = as.numeric(RID)) %>%
  arrange(RID, DateTau)

# === Verification step =======

# verify that the tau scan can be used to uniquely identify each scan pair

if (sum(duplicated(df$TauID)) == 0) {
  print('Tau scan IDs are unique for each row')
} else {
  stop('Tau scans are NOT unique for each row - revisit how identify rows.')
}

# === Add CDR ======

cdr.record <- cdr %>%
  as.data.frame() %>%
  mutate(DateCDR = ifelse(is.na(EXAMDATE),
                          get.examdate.from.registry(cdr),
                          EXAMDATE),
         DateCDR = as_datetime(ymd(DateCDR))) %>%
  select(RID, DateCDR, CDGLOBAL, CDRSB) %>%
  rename(CDRGlobal=CDGLOBAL, CDRSumBoxes=CDRSB) %>%
  drop_na(CDRGlobal)

cdr.df <- left_join(df, cdr.record, by='RID') %>%
  mutate(DiffMeanImagingDateCDR = as.numeric(difftime(MeanImagingDate, DateCDR, units = 'days')))

cdr.df <- group_by(cdr.df, TauID) %>%
  slice_min(order_by=abs(DiffMeanImagingDateCDR), with_ties = F) %>%
  ungroup()

bad <- is.na(cdr.df$CDRGlobal) | (abs(cdr.df$DiffMeanImagingDateCDR) > THRESHOLD.COGNITIVE.DAYS)
cdr.df[bad, c("DateCDR", "CDRGlobal", "CDRSumBoxes", "DiffMeanImagingDateCDR")] <- NA

cdr.df <- cdr.df %>%
  mutate(CDRBinned=cut(CDRGlobal, breaks=c(0, .5, 1, Inf), right=F))
levels(cdr.df$CDRBinned) <- c("0.0", "0.5", "1.0+")

df <- as.data.frame(cdr.df) %>%
  arrange(RID, DateTau)

df$Dementia <- ifelse(df$CDRGlobal >= 0.5 & ! is.na(df$CDRGlobal), 
                      'Yes',
                      'No')
df[is.na(df$CDRGlobal), 'Dementia'] <- 'Unknown'
df$Control <- ifelse(! df$AmyloidPositive & df$Dementia == 'No', 1, 0)

# === add MMSE ======

# add MMSE
mmse.adni <- mmse %>%
  mutate(DateMMSE = ifelse(is.na(EXAMDATE),
                           get.examdate.from.registry(mmse),
                           EXAMDATE),
         DateMMSE = as_datetime(ymd(DateMMSE))) %>%
  dplyr::select(RID, DateMMSE, MMSCORE) %>%
  rename(MMSE=MMSCORE)

mmse.merged <- left_join(df, mmse.adni, by='RID') %>%
  mutate(MMSEDiff = difftime(MeanImagingDate, DateMMSE, units='days')) %>%
  group_by(TauID) %>%
  slice_min(abs(MMSEDiff), with_ties = F) %>%
  ungroup() %>%
  mutate(MMSE = ifelse(! is.na(MMSEDiff) & abs(MMSEDiff) > THRESHOLD.COGNITIVE.DAYS, NA, MMSE))  

df <- mmse.merged

# === Add APOE ======

a1 <- select(apoeres, RID, APGEN1, APGEN2)
a2 <- select(apoego2, RID, APGEN1, APGEN2)
a3 <- select(apoe3, RID, APGEN1, APGEN2)

all.apoe <- do.call(rbind, list(a1, a2, a3))

all.apoe <- all.apoe %>%
  mutate(APOEGenotype=paste(
    pmin(all.apoe$APGEN1, all.apoe$APGEN2),
    pmax(all.apoe$APGEN1, all.apoe$APGEN2),
    sep='/')
  )

df <- left_join(df, all.apoe, by='RID')
df$HasE4 <- ifelse(is.na(df$APOEGenotype), NA, grepl('4', df$APOEGenotype))

# === Add neuropsych ======

# for computation of PACC, need some neuropsych & ADAS cog Q4
# see https://adni.bitbucket.io/reference/pacc.html

# 1. Neuropsych battery
nps <- neurobat %>%
  mutate(DateNeuropsych = ifelse(is.na(EXAMDATE),
                                 get.examdate.from.registry(neurobat),
                                 EXAMDATE),
         DateNeuropsych = as_datetime(ymd(DateNeuropsych))) %>%
  select(RID, DateNeuropsych, LDELTOTAL, DIGITSCOR, TRABSCOR)

df <- left_join(df, nps, by='RID') %>%
  mutate(DiffNPS = as.numeric(abs(difftime(MeanImagingDate, DateNeuropsych, units='days')))) %>%
  group_by(TauID) %>%
  slice_min(DiffNPS, with_ties = F)

bad.nps <- (df$DiffNPS > 365) | (is.na(df$DiffNPS))
df[bad.nps, c('LDELTOTAL', 'DIGITSCOR', 'TRABSCOR')] <- NA

# 2. ADAS
adascog <- adas %>%
  mutate(DateADAS = ifelse(is.na(EXAMDATE),
                           get.examdate.from.registry(adas),
                           EXAMDATE),
         DateADAS = as_datetime(ymd(DateADAS))) %>%
  select(RID, DateADAS, Q4SCORE) %>%
  rename(ADASQ4=Q4SCORE)

df <- left_join(df, adascog, by='RID') %>%
  mutate(DiffADAS = as.numeric(abs(difftime(MeanImagingDate, DateADAS, units='days')))) %>%
  group_by(TauID) %>%
  slice_min(DiffADAS, with_ties = F)

bad.adas <- (df$DiffADAS > 365) | (is.na(df$DiffADAS))
df[bad.adas, c('ADASQ4')] <- NA

# === Add demographics ======

min.ages <- ptdemog %>%
  mutate(DateDemogBL = ifelse(is.na(EXAMDATE),
                              get.examdate.from.registry(ptdemog),
                              as.character(EXAMDATE)),
         DateDemogBL = as_datetime(ymd(DateDemogBL))) %>%
  select(RID, DateDemogBL, AGE, PTGENDER) %>%
  rename(AgeBL=AGE, Sex=PTGENDER) %>%
  mutate(AgeBL = as.numeric(AgeBL)) %>%
  drop_na(AgeBL) %>%
  group_by(RID) %>%
  slice_min(DateDemogBL)

df.age <- left_join(df, min.ages, by='RID')
df.age$TimeSinceBL <- as.numeric(difftime(df.age$MeanImagingDate, df.age$DateDemogBL, units='days')) / 365.25
df.age$Age <- as.numeric(df.age$AgeBL + df.age$TimeSinceBL)

# manually add some missing ages
# these are not in ADNIMERGE::ptdemog but are in the downloaded study tables
missing.age <- df.age[is.na(df.age$Age), ]
replace.rids <- missing.age$RID

demog.csv <- read.csv(PATH.PTDEMOG.CSV) %>%
  select(RID, PTDOBMM, PTDOBYY) %>%
  mutate(DOB=as.POSIXct(paste(PTDOBMM, 15, PTDOBYY, sep='/'), format='%m/%d/%Y')) %>%
  filter(RID %in% replace.rids) %>%
  select(RID, DOB)

df.age <- left_join(df.age, demog.csv, by='RID') %>%
  mutate(Age = ifelse(
    is.na(Age),
    as.numeric(difftime(MeanImagingDate, DOB, units='days') / 365.25),
    Age)
  )

df <- select(df.age, -c(DOB, AgeBL, DateDemogBL, TimeSinceBL))

# same for some missing sex
df$Sex <- as.character(df$Sex)
sex.missing <- df[is.na(df$Sex), ]

demog.csv <- read.csv(PATH.PTDEMOG.CSV) %>%
  select(RID, PTGENDER) %>%
  filter(RID %in% sex.missing$RID) %>%
  mutate(Sex.Missing=recode(PTGENDER, `1`='Male', `2`='Female')) %>%
  select(-PTGENDER)

df.sex <- left_join(df, demog.csv, by='RID') %>%
  mutate(Sex = ifelse(is.na(Sex), Sex.Missing, Sex)) %>%
  select(-Sex.Missing)

df <- df.sex %>%
  mutate(Sex.Integer = ifelse(Sex == 'Male', 1, 0))

# === Add ICV ======

icv <- read.csv(PATH.ICV)

df <- left_join(df, icv, by = 'TauID')

# === Compute PACC =========

# done relative to baseline CN group
# this is using the modified formula recommended by ADNIMERGE R
# https://adni.bitbucket.io/reference/pacc.html

df$PACC <- compute.pacc(df,
                        pacc.columns = c('ADASQ4', 'LDELTOTAL', 'TRABSCOR', 'MMSE'),
                        cn.mask = df$Dementia == 'No',
                        higher.better = c(F, T, F, T),
                        min.required = 2)

# === remove NAs =========

na.cols <- c('Age', 'Sex', 'PACC', 'HasE4', 'CDRBinned')

df.withna <- df
df <- df %>%
  drop_na(all_of(na.cols))

# === variables for selecting ROIs ======

SUBCORTICAL = c('amygdala',
                'caudate',
                'hippocampus',
                'pallidum',
                'putamen',
                'thalamus',
                'ventraldc')
SUBCORTICAL_PAT = paste(SUBCORTICAL, collapse='|')

# === Add ROIs : AV45 ======

rois <- ucberkeleyav45 %>%
  select(WHOLECEREBELLUM_SUVR,
         matches('CTX.*SUVR') & ! contains('UNKNOWN'),
         matches(SUBCORTICAL_PAT, ignore.case = T) & contains('SUVR')) %>%
  mutate(across(contains('SUVR'), function (x) x / WHOLECEREBELLUM_SUVR)) %>%
  select(-WHOLECEREBELLUM_SUVR)

colnames(rois) <- paste('AV45', colnames(rois), sep='_')

roi.names <- data.frame(Region=colnames(rois)) %>%
  mutate(Cortical=ifelse(str_detect(tolower(Region), SUBCORTICAL_PAT), 'subcortical', 'cortical'))
write.csv(roi.names, file.path(PATH.DERIVATIVES, 'av45_regions.csv'), row.names = F)

rois$AmyloidID <- scan.id(ucberkeleyav45$RID, ucberkeleyav45$EXAMDATE)

df <- left_join(df, rois, by = 'AmyloidID')

# === Add ROIs : FTP ======

rois <- ucberkeleyav1451 %>%
  select(INFERIORCEREBELLUM_SUVR,
         matches('CTX.*SUVR') & ! contains('UNKNOWN'),
         matches(SUBCORTICAL_PAT, ignore.case = T) & contains('SUVR')) %>%
  mutate(across(contains('SUVR'), function (x) x / INFERIORCEREBELLUM_SUVR)) %>%
  select(-INFERIORCEREBELLUM_SUVR)

colnames(rois) <- paste('FTP', colnames(rois), sep='_')

roi.names <- data.frame(Region=colnames(rois)) %>%
  mutate(Cortical=ifelse(str_detect(tolower(Region), SUBCORTICAL_PAT), 'subcortical', 'cortical'))
write.csv(roi.names, file.path(PATH.DERIVATIVES, 'ftp_regions.csv'), row.names = F)

rois$TauID <- scan.id(ucberkeleyav1451$RID, ucberkeleyav1451$EXAMDATE)

df <- left_join(df, rois, by = 'TauID')

# === Add ROIs : Volume ======

rois <- ucberkeleyav1451 %>%
  select(matches('CTX.*VOLUME') & ! contains('UNKNOWN'),
         matches(SUBCORTICAL_PAT, ignore.case = T) & contains('VOLUME'))

roi.names <- data.frame(Region=colnames(rois)) %>%
  mutate(Cortical=ifelse(str_detect(tolower(Region), SUBCORTICAL_PAT), 'subcortical', 'cortical'))
write.csv(roi.names, file.path(PATH.DERIVATIVES, 'gm_regions.csv'), row.names = F)

rois$TauID <- scan.id(ucberkeleyav1451$RID, ucberkeleyav1451$EXAMDATE)

df <- left_join(df, rois, by = 'TauID') %>%
  mutate(across(ends_with('_VOLUME'), function (x) (x * 1000 / ICV)))

# === save ========

dir.create(dirname(PATH.OUTPUT), showWarnings = F)
write.csv(df, PATH.OUTPUT, quote = F, na = '', row.names = F)
