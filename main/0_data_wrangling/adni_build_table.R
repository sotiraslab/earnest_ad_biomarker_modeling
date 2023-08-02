# === imports ======

sh <- suppressPackageStartupMessages

sh(library(ADNIMERGE))
sh(library(ggplot2))
sh(library(lme4))
sh(library(lubridate))
sh(library(stringr))
sh(library(this.path))
sh(library(tidyverse))

# === Set working directory ======

setwd(this.dir())

# === Set paths ======

PATH.PTDEMOG.CSV <- '../../data/rawdata/PTDEMOG.csv'

PATH.PACC.SCRIPT <- '../../scripts/pacc.R'
PATH.EXAMDATE.SCRIPT <- '../../scripts/adni_examdate.R'

PATH.OUTPUT <- '../../data/derivatives/adni_base_table.csv'
PATH.DERIVATIVES = '../../data/derivatives'

source(PATH.PACC.SCRIPT)
source(PATH.EXAMDATE.SCRIPT)

dir.create(PATH.DERIVATIVES, showWarnings = F)

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
         TauID = scan.id(RID, EXAMDATE),
         META_TEMPORAL_TAU = META_TEMPORAL_SUVR,
         BRAAK1_TAU = BRAAK1_SUVR,
         BRAAK34_TAU = BRAAK34_SUVR,
         BRAAK56_TAU = BRAAK56_SUVR,
         META_TEMPORAL_VOL = META_TEMPORAL_VOLUME) %>%
  select(RID, DateTau, TauID, META_TEMPORAL_TAU,
         BRAAK1_TAU, BRAAK34_TAU, BRAAK56_TAU, META_TEMPORAL_VOL) %>%
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
         Centiloid = as.numeric((196.9*SUMMARYSUVR_WHOLECEREBNORM) - 196.03),
         AMYLOID_COMPOSITE = SUMMARYSUVR_WHOLECEREBNORM) %>%
  select(RID, DateAmyloid, AmyloidID, AmyloidTracer,
         AmyloidPositive, AmyloidID, Centiloid, AMYLOID_COMPOSITE)

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
  filter(abs(DiffTauAmyloid) < THRESHOLD.IMAGING.DAYS) %>%
  ungroup()

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
  rename(CDR=CDGLOBAL, CDRSumBoxes=CDRSB) %>%
  drop_na(CDR)

cdr.df <- left_join(df, cdr.record, by='RID') %>%
  mutate(DiffMeanImagingDateCDR = as.numeric(difftime(MeanImagingDate, DateCDR, units = 'days')))

cdr.df <- group_by(cdr.df, TauID) %>%
  slice_min(order_by=abs(DiffMeanImagingDateCDR), with_ties = F) %>%
  ungroup()

bad <- is.na(cdr.df$CDR) | (abs(cdr.df$DiffMeanImagingDateCDR) > THRESHOLD.COGNITIVE.DAYS)
cdr.df[bad, c("DateCDR", "CDR", "CDRSumBoxes", "DiffMeanImagingDateCDR")] <- NA

cdr.df <- cdr.df %>%
  mutate(CDRBinned=cut(CDR, breaks=c(0, .5, 1, Inf), right=F))
levels(cdr.df$CDRBinned) <- c("0.0", "0.5", "1.0+")

df <- as.data.frame(cdr.df) %>%
  arrange(RID, DateTau)

df$Dementia <- ifelse(df$CDR >= 0.5 & ! is.na(df$CDR), 
                      'Yes',
                      'No')
df[is.na(df$CDR), 'Dementia'] <- 'Unknown'
df$Control <- ifelse(! df$AmyloidPositive & df$Dementia == 'No', 1, 0)

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
  slice_min(DateDemogBL) %>%
  ungroup

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

icvs <- adnimerge %>%
  select(RID, ICV) %>%
  group_by(RID) %>%
  summarise(ICV=mean(ICV, na.rm=T)) %>%
  ungroup() %>%
  mutate(RID = as.numeric(RID))

df.icv <- left_join(df, icvs, by='RID')

male.icv <- mean(df.icv[df$Sex == 'Male', 'ICV'], na.rm = T)
female.icv <- mean(df.icv[df$Sex == 'Female', 'ICV'], na.rm = T)

df.icv$ICV <- ifelse(df.icv$Sex == 'Male' & is.na(df.icv$ICV), male.icv, df.icv$ICV)
df.icv$ICV <- ifelse(df.icv$Sex == 'Female' & is.na(df.icv$ICV), female.icv, df.icv$ICV)

df <- df.icv

# === Add PACC =========

# for computation of PACC, need some neuropsych & ADAS cog Q4
# see https://adni.bitbucket.io/reference/pacc.html

subs <- df %>%
  select(RID)

# 1. Neuropsych battery
nps <- neurobat %>%
  mutate(DateNeuropsych = ifelse(is.na(EXAMDATE),
                                 get.examdate.from.registry(neurobat),
                                 EXAMDATE),
         DateNeuropsych = as_datetime(ymd(DateNeuropsych)),
         LDELTOTAL = as.numeric(LDELTOTAL)) %>%
  select(RID, DateNeuropsych, LDELTOTAL, DIGITSCOR, TRABSCOR)

nps <- left_join(subs, nps, by='RID')

# 2. ADAS
adascog <- adas %>%
  mutate(DateNeuropsych = ifelse(is.na(EXAMDATE),
                                 get.examdate.from.registry(adas),
                                 EXAMDATE),
         DateNeuropsych = as_datetime(ymd(DateNeuropsych))) %>%
  select(RID, DateNeuropsych, Q4SCORE) %>%
  rename(ADASQ4=Q4SCORE)

adascog <- left_join(subs, adascog, by='RID')

# 3. MMSE 
mmse.adni <- mmse %>%
  mutate(DateNeuropsych = ifelse(is.na(EXAMDATE),
                                 get.examdate.from.registry(mmse),
                                 EXAMDATE),
         DateNeuropsych = as_datetime(ymd(DateNeuropsych))) %>%
  dplyr::select(RID, DateNeuropsych, MMSCORE) %>%
  rename(MMSE=MMSCORE)

mmse.adni <- left_join(subs, mmse.adni, by='RID')

# 4. Join all
pacc.df <- df %>%
  select(RID) %>%
  full_join(nps, by='RID') %>%
  full_join(adascog, by=c('RID', 'DateNeuropsych')) %>%
  full_join(mmse.adni, by=c('RID', 'DateNeuropsych')) %>%
  arrange(RID, DateNeuropsych) %>%
  drop_na(DateNeuropsych)

# Group assessments which occur at nearby dates

LINK.THR <- 60 # days

pacc.df <- pacc.df %>%
  group_by(RID) %>%
  mutate(NPSDateDiff = c(0, as.numeric(difftime(DateNeuropsych[-1],
                                                DateNeuropsych[-n()],
                                                units="days")))) %>%
  ungroup()

pacc.df$NPSDateGroup <- cumsum(abs(pacc.df$NPSDateDiff) > LINK.THR)

select.first <- function(col) {
  return(first(col[! is.na(col)]))
}


# this also gets used later for computing longitudinal change in pacc
pacc.df.group <- pacc.df %>%
  group_by(RID, NPSDateGroup) %>%
  summarise(
    DatePACC = as_datetime(mean(DateNeuropsych)),
    DatePACCMin = min(DateNeuropsych),
    DatePACCMax = max(DateNeuropsych),
    DatePACCDiff = as.numeric(difftime(DatePACCMax, DatePACCMin, units='days')),
    LDELTOTAL = as.numeric(select.first(LDELTOTAL)),
    DIGITSCOR = as.numeric(select.first(DIGITSCOR)),
    TRABSCOR = as.numeric(select.first(TRABSCOR)),
    ADASQ4 = as.numeric(select.first(ADASQ4)),
    MMSE = as.numeric(select.first(MMSE))
  ) %>%
  ungroup()

# merge into df
df <- left_join(df, pacc.df.group, by="RID") %>%
  mutate(DiffImagingPACC = as.numeric(difftime(MeanImagingDate, DatePACC, units = 'days'))) %>%
  group_by(TauID) %>%
  slice_min(abs(DiffImagingPACC), with_ties = F) %>%
  filter(abs(DiffImagingPACC) < THRESHOLD.COGNITIVE.DAYS) %>%
  ungroup() %>%
  arrange(RID)

# Compute PACC
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

# === Compute longitudinal change in PACC =======

# this will only be available for some

cn.data <- df %>%
  filter(Dementia == "No")

pacc.long <- df %>%
  select(RID, DatePACC, Age, CDRBinned) %>%
  rename(DatePACC.BL=DatePACC) %>%
  left_join(pacc.df.group, by="RID") %>%
  group_by(RID) %>%
  filter(DatePACC >= DatePACC.BL) %>%
  ungroup()

pacc.long$PACC <- compute.pacc(pacc.long,
                               pacc.columns = c('ADASQ4', 'LDELTOTAL', 'TRABSCOR', 'MMSE'),
                               cn.data = cn.data,
                               higher.better = c(F, T, F, T),
                               min.required = 2)
pacc.long <- pacc.long %>%
  filter(! is.na(PACC)) %>%
  group_by(RID) %>%
  filter(n() >= 2) %>%
  mutate(DeltaPACCDate = as.numeric(difftime(DatePACC, DatePACC.BL, units='days')) / 365.25,
         Long.Age = Age + DeltaPACCDate) %>% 
  ungroup()

m <- lmer(PACC ~ DeltaPACCDate + (1+DeltaPACCDate|RID), data=pacc.long)
pacc.long$PACC.LMER.Predict <- predict(m, pacc.long)

ggplot(pacc.long, aes(x=Long.Age, y=PACC)) +
  geom_point(aes(color=CDRBinned), alpha = .7) + 
  geom_line(aes(y=PACC.LMER.Predict, group=RID, color=CDRBinned), alpha= .7)

ggsave("adni_longitudinal_pacc_model.png", width=8, height=6, units='in')

coefs <- coef(m)$RID %>%
  select(DeltaPACCDate) %>%
  dplyr::rename(DeltaPACC=DeltaPACCDate) %>%
  rownames_to_column(var="RID") %>%
  mutate(RID=as.numeric(RID))

df <- left_join(df, coefs, by='RID')

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

# === compute hippocampus volume ======

df$HIPPOCAMPUS_VOL = df$LEFT_HIPPOCAMPUS_VOLUME + df$RIGHT_HIPPOCAMPUS_VOLUME

# === CSF markers =======

# looks like the CSF markers are too sparse for this dataset
# and the ADNI1/2/GO vs. ADNI3 assays have values on very different scales
# might need to do some more reading  to figure this out, but leaving be for now

# a <- df
# b <- df
# 
# old.csf <- upennbiomk_master %>%
#   mutate(RID = as.numeric(RID),
#          DateCSF = as_datetime(mdy(DRAWDTE))) %>%
#   select(RID, DateCSF, ABETA, TAU, PTAU)
# 
# a <- left_join(a, old.csf, by='RID') %>%
#   mutate(DiffTauCSF = as.numeric(difftime(DateTau, DateCSF, units='days')) / 365.25) %>%
#   group_by(TauID) %>%
#   slice_min(abs(DiffTauCSF), with_ties = F) %>%
#   filter(abs(DiffTauAmyloid) < THRESHOLD.IMAGING.DAYS) %>%
#   ungroup() %>%
#   select(TauID, DateTau, ABETA, TAU, PTAU)
# 
# bmk10 <- upennbiomk10 %>%
#   select(RID, DRAWDATE, ABETA40, ABETA42, TAU, PTAU) %>%
#   rename(DateCSF=DRAWDATE) %>%
#   mutate(DateCSF=as_datetime(mdy(DateCSF)))
# 
# bmk12 <- upennbiomk12_2020 %>%
#   select(RID, EXAMDATE, AB40, ABETA, TAU, PTAU) %>%
#   rename(DateCSF=EXAMDATE,
#          ABETA40=AB40,
#          ABETA42=ABETA) %>%
#   mutate(DateCSF = as_datetime(ymd(DateCSF)))
# 
# new.csf <- rbind(bmk10, bmk12)
# 
# b <- left_join(b, new.csf, by='RID') %>%
#   mutate(DiffTauCSF = as.numeric(difftime(DateTau, DateCSF, units='days')) / 365.25) %>%
#   group_by(TauID) %>%
#   slice_min(abs(DiffTauCSF), with_ties = F) %>%
#   filter(abs(DiffTauAmyloid) < THRESHOLD.IMAGING.DAYS) %>%
#   ungroup() %>%
#   select(TauID, DateTau, ABETA42, TAU, PTAU)

# === save ========

dir.create(dirname(PATH.OUTPUT), showWarnings = F)
write.csv(df, PATH.OUTPUT, quote = F, na = '', row.names = F)
