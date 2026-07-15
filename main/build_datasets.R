# === imports ======

sh <- suppressPackageStartupMessages

sh(library(ADNIMERGE))
sh(library(car))
sh(library(factoextra))
sh(library(ggplot2))
sh(library(lme4))
sh(library(lubridate))
sh(library(stringr))
sh(library(tableone))
sh(library(this.path))
sh(library(tidyverse))

# === Set working directory ======

setwd(this.dir())

# === Set paths ======

PATH.LOAD.ADNI <- '../rscripts/load_adni_table.R'
PATH.PACC.SCRIPT <- '../rscripts/pacc.R'
PATH.EXAMDATE.SCRIPT <- '../rscripts/adni_examdate.R'

PATH.OUTPUT <- '.'

source(PATH.LOAD.ADNI)
source(PATH.PACC.SCRIPT)
source(PATH.EXAMDATE.SCRIPT)

# === Create output folder ======

outfolder <- file.path('datasets')
dir.create(outfolder, showWarnings = F)

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

# check for subjects with both PVC and non-PVC
tau.nopvc <- load.adni.table('TAU_6MM', 'inputs')
tau.pvc <- load.adni.table('TAUPVC_6MM', 'inputs')

tau.both <- inner_join(tau.nopvc[, c('RID', 'SCANDATE')],
                       tau.pvc[, c('RID', 'SCANDATE')])

tau <- tau.both %>%
  mutate(DateTau = as_datetime(ymd(SCANDATE)),
         TauID = scan.id(RID, SCANDATE)) %>%
  select(RID, DateTau, TauID) %>%
  group_by(RID) %>%
  slice_min(DateTau, with_ties = F) %>%
  ungroup()

amy <- load.adni.table('AMY_6MM', 'inputs')

# centiloid is manually calculated using
# equation in UC Berkeley Amyloid PET docs PDF
# gives a continuous measure rather than integer
av45 <- amy %>%
  filter(TRACER == 'FBP') %>%
  mutate(DateAmyloid=as_datetime(ymd(SCANDATE)),
         AmyloidID = scan.id(RID, SCANDATE),
         AmyloidTracer = 'AV45',
         AmyloidPositive = AMYLOID_STATUS,
         AmyloidComposite = SUMMARY_SUVR,
         Centiloid = (188.22 * SUMMARY_SUVR) - 189.16) %>%
  select(RID, DateAmyloid, AmyloidID, AmyloidTracer,
         AmyloidPositive, AmyloidID, AmyloidComposite, Centiloid)

link.av45 <- left_join(tau, av45, by='RID') %>%
  mutate(DiffTauAmyloid = as.numeric(difftime(DateTau, DateAmyloid, units = 'days'))) %>%
  group_by(TauID) %>%
  slice_min(abs(DiffTauAmyloid), with_ties = F) %>%
  filter(abs(DiffTauAmyloid) < THRESHOLD.IMAGING.DAYS) %>%
  ungroup()

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

cdr.long <- left_join(df, cdr.record, by='RID') %>%
  mutate(DiffMeanImagingDateCDR = as.numeric(difftime(MeanImagingDate, DateCDR, units = 'days')),
         CDRBinned = cut(CDR, breaks=c(0, .5, 1, Inf), right=F))

cdr.df <- group_by(cdr.long, TauID) %>%
  slice_min(order_by=abs(DiffMeanImagingDateCDR), with_ties = F) %>%
  ungroup()

bad <- is.na(cdr.df$CDR) | (abs(cdr.df$DiffMeanImagingDateCDR) > THRESHOLD.COGNITIVE.DAYS)
cdr.df[bad, c("DateCDR", "CDR", "CDRSumBoxes", "DiffMeanImagingDateCDR", 'CDRBinned')] <- NA

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
  select(RID, DateDemogBL, AGE, PTGENDER, PTEDUCAT, PTRACCAT, PTETHCAT) %>%
  rename(AgeBL=AGE, Sex=PTGENDER, Education=PTEDUCAT, Race=PTRACCAT, Ethnicity=PTETHCAT) %>%
  mutate(AgeBL = as.numeric(AgeBL),
         Hispanic = as.numeric(Ethnicity == 'Hispanic or Latino'),
         Race = recode(
           Race,
           'American Indian or Alaskan Native'='Other',
           'Black or African American' = 'Black',
           'More than one race' = 'Other')
         ) %>%
  drop_na(AgeBL) %>%
  group_by(RID) %>%
  slice_min(DateDemogBL) %>%
  ungroup()

df.age <- left_join(df, min.ages, by='RID')
df.age$TimeSinceBL <- as.numeric(difftime(df.age$MeanImagingDate, df.age$DateDemogBL, units='days')) / 365.25
df.age$Age <- as.numeric(df.age$AgeBL + df.age$TimeSinceBL)

df <- df.age

# === Add ICV ======

icvs <- adnimerge %>%
  select(RID, ICV) %>%
  group_by(RID) %>%
  summarise(ICV=mean(ICV, na.rm=T)) %>%
  ungroup() %>%
  mutate(RID = as.numeric(RID))

df.icv <- left_join(df, icvs, by='RID')

# ICVs are imputed for those missing
male.icv <- mean(df.icv[df$Sex == 'Male', 'ICV'], na.rm = T)
female.icv <- mean(df.icv[df$Sex == 'Female', 'ICV'], na.rm = T)

df.icv$ICV <- ifelse(df.icv$Sex == 'Male' & is.na(df.icv$ICV), male.icv, df.icv$ICV)
df.icv$ICV <- ifelse(df.icv$Sex == 'Female' & is.na(df.icv$ICV), female.icv, df.icv$ICV)

df <- df.icv

# === Add MMSE =====

mm <- mmse %>%
  mutate(DateMMSE = as_datetime(ymd(USERDATE))) %>%
  select(RID, DateMMSE, MMSCORE_EDC) %>%
  rename(MMSE=MMSCORE_EDC)

mmse.long <- left_join(df, mm, by='RID') %>%
  mutate(DiffImagingMMSE = as.numeric(difftime(MeanImagingDate, DateMMSE, units='days')))

df.mmse <- mmse.long %>%
  group_by(RID) %>%
  slice_min(order_by=abs(DiffImagingMMSE), with_ties = F) %>%
  ungroup() %>%
  filter(abs(DiffImagingMMSE) < THRESHOLD.COGNITIVE.DAYS)

df <- df.mmse

# === Add PHC =========

adsp <- load.adni.table('ADSP_PHC', 'inputs')

adsp <- adsp %>%
  mutate(DateADSP = as_datetime(ymd(EXAMDATE))) %>%
  select(RID, DateADSP, PHC_MEM, PHC_EXF, PHC_LAN, PHC_VSP)

#save for calculating longitudinal change
df.adsp.long <- left_join(df, adsp, by='RID') %>%
  mutate(DiffImagingADSP = as.numeric(difftime(MeanImagingDate, DateADSP, units = 'days')))

df.adsp <- df.adsp.long %>%
  group_by(RID) %>%
  slice_min(order_by=abs(DiffImagingADSP), with_ties = F) %>%
  ungroup() %>%
  filter(abs(DiffImagingADSP) < THRESHOLD.COGNITIVE.DAYS)

df <- df.adsp

# === Add ADNI-composites =========

psych <- uwnpsychsum %>%
  mutate(DateUWPSYCH = as_datetime(ymd(EXAMDATE))) %>%
  select(RID, DateUWPSYCH, ADNI_MEM, ADNI_EF, ADNI_LAN, ADNI_VS)

#save for calculating longitudinal change
psych.long <- left_join(df, psych, by='RID') %>%
  mutate(DiffImagingUWPSYCH = as.numeric(difftime(MeanImagingDate, DateUWPSYCH, units = 'days')))

df.psych <- psych.long %>%
  group_by(RID) %>%
  slice_min(order_by=abs(DiffImagingUWPSYCH), with_ties = F) %>%
  ungroup() %>%
  filter(abs(DiffImagingUWPSYCH) < THRESHOLD.COGNITIVE.DAYS)

df <- df.psych

# === Add PACC =========

# Neuropsych battery - LDELTOTAL, TRABSCOR
nps <- neurobat %>%
  mutate(DateNeuropsych = ifelse(is.na(EXAMDATE),
                                 get.examdate.from.registry(neurobat),
                                 EXAMDATE),
         DateNeuropsych = as_datetime(ymd(DateNeuropsych))) %>%
  select(RID, DateNeuropsych, LDELTOTAL, TRABSCOR)

# ADAS - Q4
adascog <- adas %>%
  mutate(DateADAS = ifelse(is.na(EXAMDATE),
                           get.examdate.from.registry(adas),
                           EXAMDATE),
         DateADAS = as_datetime(ymd(DateADAS))) %>%
  select(RID, DateADAS, Q4SCORE) %>%
  rename(ADASQ4=Q4SCORE)

# MMSE
this.mmse <- mmse %>%
  mutate(DateMMSE = as_datetime(ymd(USERDATE))) %>%
  select(RID, DateMMSE, MMSCORE_EDC) %>%
  rename(MMSE=MMSCORE_EDC)

# COMBINE
pacc <- left_join(this.mmse, nps, by='RID', relationship = "many-to-many") %>%
  mutate(DiffMMSENPS = as.numeric(difftime(DateMMSE, DateNeuropsych, units='days')))

pacc <- pacc %>%
  group_by(RID, DateMMSE) %>%
  slice_min(order_by=abs(DiffMMSENPS), with_ties = F) %>%
  ungroup() %>%
  filter(abs(DiffMMSENPS) < THRESHOLD.COGNITIVE.DAYS)

pacc <- left_join(pacc, adascog, by='RID', relationship = "many-to-many") %>%
  mutate(DiffMMSEADAS = as.numeric(difftime(DateMMSE, DateADAS, units='days')))

pacc <- pacc %>%
  group_by(RID, DateMMSE) %>%
  slice_min(order_by=abs(DiffMMSEADAS), with_ties = F) %>%
  ungroup() %>%
  filter(abs(DiffMMSEADAS) < THRESHOLD.COGNITIVE.DAYS) %>%
  mutate(DatePACC = DateMMSE) %>%
  select(RID, DatePACC, MMSE, LDELTOTAL, TRABSCOR, ADASQ4)

# Merge into base DF for computing PACC
pacc.long <- df %>%
  select(RID, DateTau, DateMMSE, AmyloidPositive, CDRBinned) 

pacc.long <- left_join(pacc.long, pacc, by = 'RID') %>%
  filter(DatePACC >= DateMMSE) %>%
  group_by(RID) %>%
  mutate(Baseline = row_number( )== 1) %>%
  ungroup()

cn.mask <- (
  pacc.long$Baseline & 
    (pacc.long$AmyloidPositive == 0 & !is.na(pacc.long$AmyloidPositive)) & 
    (pacc.long$CDRBinned == '0.0' & !is.na(pacc.long$CDRBinned))
  )

pacc.long$PACC <- compute.pacc(
  df = pacc.long,
  pacc.columns = c('ADASQ4', 'LDELTOTAL', 'TRABSCOR', 'MMSE'),
  cn.mask = cn.mask,
  higher.better = c(F, T, F, T),
  min.required = 2
)

# Merge into original DF
pacc.merge <- pacc.long %>%
  select(RID, DatePACC, PACC)

pacc.leftjoin <- left_join(df, pacc.merge, by='RID') %>%
  mutate(DiffImagingPACC = as.numeric(difftime(MeanImagingDate, DatePACC, units='days')))

df.pacc <- pacc.leftjoin %>%
  group_by(RID) %>%
  slice_min(order_by=abs(DiffImagingPACC), with_ties = F) %>%
  ungroup() %>%
  filter(abs(DiffImagingPACC) < THRESHOLD.COGNITIVE.DAYS)

df <- df.pacc

# === Helper for computing longitudinal change =======

calc.longitudinal.change <- function(baseline, longitudinal,
                                     variable, date.column,
                                     id.column='RID', age.column='Age',
                                     plot.by='CDRBinned',
                                     fixed.endpoint.years = NULL,
                                     fixed.endpoint.gap = .5,
                                     destination.column = NULL,
                                     do.bootstrap.se = F,
                                     save.model.path = NULL,
                                     save.data.path = NULL) {

  joiner <- longitudinal %>%
    select(!!id.column, !!date.column, !!variable, !!plot.by) %>%
    rename(ID=!!id.column, DATE=!!date.column, VAR=!!variable)
  
  long.data <- baseline %>%
    select(!!id.column, !!date.column, !!age.column, !!plot.by) %>%
    rename(ID=!!id.column, DATE.BL=!!date.column, AGE=!!age.column, PLOTBY=!!plot.by) %>%
    left_join(joiner, by='ID') %>%
    group_by(ID) %>%
    mutate(DELTA = as.numeric(difftime(DATE, DATE.BL, units='days')) / 365.25,
           LONG.AGE = AGE + DELTA) %>%
    filter(DATE >= DATE.BL) %>%
    filter(n() >= 2) %>%
    drop_na(VAR) %>%
    ungroup()
  
  if (! is.null(fixed.endpoint.years)) {
    long.data <- long.data %>%
      mutate(
        InWindow = (DELTA >= (fixed.endpoint.years - fixed.endpoint.gap)) & 
          (DELTA <= (fixed.endpoint.years + fixed.endpoint.gap))
      ) %>%
      group_by(ID) %>%
      filter(any(InWindow)) %>%
      ungroup() %>%
      filter(DELTA <= (fixed.endpoint.years + fixed.endpoint.gap))
  }
  
  # longitudinal modelling
  m <- lmer(VAR ~ DELTA + (1+DELTA|ID), data=long.data)
  long.data$VAR.PREDICT <- predict(m, long.data)
  
  p <- ggplot(long.data, aes(x=LONG.AGE, y=VAR)) +
    geom_point(aes(color=PLOTBY), alpha = .7) +
    geom_line(aes(y=VAR.PREDICT, group=ID, color=PLOTBY), alpha= .7) +
    ggtitle(variable)
  
  print(p)
  
  if (is.null(destination.column)) {
    final.name <- paste('Delta', variable, sep='')
  } else {
    final.name <- destination.column
  }
  
  # Build output
  coefs <- coef(m)$ID %>%
    select(DELTA) %>%
    rownames_to_column(var='ID') %>%
    mutate(ID = as.numeric(ID))
  colnames(coefs) <- c(id.column, final.name)
  
  result <- left_join(baseline, coefs, by=id.column)
  
  # Bootstrap to estimate model uncertainty
  if (do.bootstrap.se) {
    boot.slopes <- bootMer(x = m,
                           FUN = function(m) coef(m)$ID[, 'DELTA'],
                           nsim = 1000, type = 'parametric', verbose = T)
    subject.slope.se <- apply(boot.slopes$t, 2, sd)
    toadd <- data.frame(ID = coefs[[id.column]], BOOT = subject.slope.se)
    colnames(toadd) <- c(id.column, str_c(final.name, 'BootSE'))
    
    result <- left_join(result, toadd)
  }
  
  # Save other outputs
  if (! is.null(save.model.path)) {
    saveRDS(m, file = save.model.path)
  }
  
  if (! is.null(save.data.path)) {
    write.csv(long.data, save.data.path, row.names = F)
  }
  
  return (result)
}

# # ==== TEST ======
# 
# baseline <- df
# longitudinal <- mmse.long
# variable <- 'MMSE'
# date.column <- 'DateMMSE'
# id.column<-'RID'
# age.column<-'Age'
# plot.by<-'CDRBinned'
# fixed.endpoint.years <- NULL
# fixed.endpoint.gap <- .5
# destination.column <- NULL
# do.bootstrap.se <- T
# save.model.path <- NULL
# save.data.path <- NULL
# 
# joiner <- longitudinal %>%
#   select(!!id.column, !!date.column, !!variable, !!plot.by) %>%
#   rename(ID=!!id.column, DATE=!!date.column, VAR=!!variable)
# 
# long.data <- baseline %>%
#   select(!!id.column, !!date.column, !!age.column, !!plot.by) %>%
#   rename(ID=!!id.column, DATE.BL=!!date.column, AGE=!!age.column, PLOTBY=!!plot.by) %>%
#   left_join(joiner, by='ID') %>%
#   group_by(ID) %>%
#   mutate(DELTA = as.numeric(difftime(DATE, DATE.BL, units='days')) / 365.25,
#          LONG.AGE = AGE + DELTA) %>%
#   filter(DATE >= DATE.BL) %>%
#   filter(n() >= 2) %>%
#   drop_na(VAR) %>%
#   ungroup()
# 
# # longitudinal modelling
# m <- lmer(VAR ~ DELTA + (1+DELTA|ID), data=long.data)
# long.data$VAR.PREDICT <- predict(m, long.data)
# 
# if (is.null(destination.column)) {
#   final.name <- paste('Delta', variable, sep='')
# } else {
#   final.name <- destination.column
# }
# 
# # Create output
# coefs <- coef(m)$ID %>%
#   select(DELTA) %>%
#   rownames_to_column(var='ID') %>%
#   mutate(ID = as.numeric(ID))
# colnames(coefs) <- c(id.column, final.name)
# 
# result <- left_join(baseline, coefs, by=id.column)


# === Compute longitudinal changes =====

# MMSE
df <- calc.longitudinal.change(
  baseline = df,
  longitudinal = mmse.long,
  variable = 'MMSE',
  date.column = 'DateMMSE',
  do.bootstrap.se = T,
  save.data.path = file.path(outfolder, 'mmse_lmem_fit_data.csv'),
  save.model.path = file.path(outfolder, 'mmse_lmem_model.rds')
)

# CDRSB
df <- calc.longitudinal.change(
  baseline = df,
  longitudinal = cdr.long,
  variable = 'CDRSumBoxes',
  date.column = 'DateCDR'
)

# PACC
df <- calc.longitudinal.change(
  baseline = df,
  longitudinal = pacc.long,
  variable = 'PACC',
  date.column = 'DatePACC'
)

# ADSP
df <- calc.longitudinal.change(
  baseline = df,
  longitudinal = df.adsp.long,
  variable = 'PHC_MEM',
  date.column = 'DateADSP'
)

df <- calc.longitudinal.change(
  baseline = df,
  longitudinal = df.adsp.long,
  variable = 'PHC_EXF',
  date.column = 'DateADSP'
)

df <- calc.longitudinal.change(
  baseline = df,
  longitudinal = df.adsp.long,
  variable = 'PHC_LAN',
  date.column = 'DateADSP'
)

df <- calc.longitudinal.change(
  baseline = df,
  longitudinal = df.adsp.long,
  variable = 'PHC_VSP',
  date.column = 'DateADSP'
)

# # UW Psych composites
# df <- calc.longitudinal.change(
#   baseline = df,
#   longitudinal = psych.long,
#   variable = 'ADNI_MEM',
#   date.column = 'DateUWPSYCH'
# )
# 
# df <- calc.longitudinal.change(
#   baseline = df,
#   longitudinal = psych.long,
#   variable = 'ADNI_EF',
#   date.column = 'DateUWPSYCH'
# )
# 
# df <- calc.longitudinal.change(
#   baseline = df,
#   longitudinal = psych.long,
#   variable = 'ADNI_LAN',
#   date.column = 'DateUWPSYCH'
# )
# 
# df <- calc.longitudinal.change(
#   baseline = df,
#   longitudinal = psych.long,
#   variable = 'ADNI_VS',
#   date.column = 'DateUWPSYCH'
# )

# === MMSE changes within fixed windows =====

# MMSE - 2 year
df <- calc.longitudinal.change(
  baseline = df,
  longitudinal = mmse.long,
  variable = 'MMSE',
  date.column = 'DateMMSE',
  fixed.endpoint.years = 2,
  fixed.endpoint.gap = 0.5,
  destination.column = 'DeltaMMSE2Year'
)

df <- calc.longitudinal.change(
  baseline = df,
  longitudinal = mmse.long,
  variable = 'MMSE',
  date.column = 'DateMMSE',
  fixed.endpoint.years = 3,
  fixed.endpoint.gap = 0.5,
  destination.column = 'DeltaMMSE3Year'
)

df <- calc.longitudinal.change(
  baseline = df,
  longitudinal = mmse.long,
  variable = 'MMSE',
  date.column = 'DateMMSE',
  fixed.endpoint.years = 4,
  fixed.endpoint.gap = 0.5,
  destination.column = 'DeltaMMSE4Year'
)


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

base.rois <- load.adni.table('AMY', 'inputs')

rois <- base.rois %>%
  select(matches('CTX_(LH|RH)_.*SUVR') & ! contains('UNKNOWN'),
         matches(SUBCORTICAL_PAT, ignore.case = T) & contains('SUVR') & matches('^(LEFT|RIGHT)'))

colnames(rois) <- paste('AV45', colnames(rois), sep='_')

roi.names <- data.frame(Region=colnames(rois)) %>%
  mutate(Cortical=ifelse(str_detect(tolower(Region), SUBCORTICAL_PAT), 'subcortical', 'cortical'))
write.csv(roi.names, file.path(outfolder, 'av45_regions.csv'), row.names = F)

rois$AmyloidID <- scan.id(base.rois$RID, base.rois$SCANDATE)

df <- left_join(df, rois, by = 'AmyloidID')

# === Add ROIs : FTP (non-PVC) ======

base.rois <- load.adni.table('TAU_6MM', 'inputs')

rois <- base.rois %>%
  select(matches('CTX_(LH|RH)_.*SUVR') & ! contains('UNKNOWN'),
         matches(SUBCORTICAL_PAT, ignore.case = T) & contains('SUVR') & matches('^(LEFT|RIGHT)'))

colnames(rois) <- paste('FTP', colnames(rois), sep='_')

roi.names <- data.frame(Region=colnames(rois)) %>%
  mutate(Cortical=ifelse(str_detect(tolower(Region), SUBCORTICAL_PAT), 'subcortical', 'cortical'))
write.csv(roi.names, file.path(outfolder, 'ftp_regions.csv'), row.names = F)

rois$TauID <- scan.id(base.rois$RID, base.rois$SCANDATE)
rois$META_TEMPORAL_TAU <- base.rois$META_TEMPORAL_SUVR

df <- left_join(df, rois, by = 'TauID')

# === Add ROIs : Volume ======

base.rois <- load.adni.table('TAU_6MM', 'inputs')

rois <- base.rois %>%
  select(matches('CTX_(LH|RH)_.*VOLUME') & ! contains('SUVR'),
         matches(SUBCORTICAL_PAT, ignore.case = T) & contains('VOLUME') & matches('^(LEFT|RIGHT)'))

roi.names <- data.frame(Region=colnames(rois)) %>%
  mutate(Cortical=ifelse(str_detect(tolower(Region), SUBCORTICAL_PAT), 'subcortical', 'cortical'))
write.csv(roi.names, file.path(outfolder, 'gm_regions.csv'), row.names = F)

rois$TauID <- scan.id(base.rois$RID, base.rois$SCANDATE)
rois$META_TEMPORAL_VOL <- base.rois$META_TEMPORAL_VOLUME

df <- left_join(df, rois, by = 'TauID') %>%
  mutate(across(ends_with('_VOLUME'), function (x) (x * 1000 / ICV)))

# === Add ROIs : PVC Tau ======

base.rois <- load.adni.table('TAUPVC_6MM', 'inputs')

rois <- base.rois %>%
  select(matches('CTX_(LH|RH)_.*SUVR') & ! contains('UNKNOWN'),
         matches(SUBCORTICAL_PAT, ignore.case = T) & contains('SUVR') & matches('^(LEFT|RIGHT)'))

colnames(rois) <- paste('FTPPVC', colnames(rois), sep='_')

roi.names <- data.frame(Region=colnames(rois)) %>%
  mutate(Cortical=ifelse(str_detect(tolower(Region), SUBCORTICAL_PAT), 'subcortical', 'cortical'))
write.csv(roi.names, file.path(outfolder, 'ftppvc_regions.csv'), row.names = F)

rois$TauID <- scan.id(base.rois$RID, base.rois$SCANDATE)
rois$META_TEMPORAL_TAUPVC <- base.rois$META_TEMPORAL_SUVR

df <- left_join(df, rois, by = 'TauID')

# === compute hippocampus volume ======

df$HIPPOCAMPUS_VOL = df$LEFT_HIPPOCAMPUS_VOLUME + df$RIGHT_HIPPOCAMPUS_VOLUME

# === compute bilateral PET uptakes =======

# needed for some staging models

bilateral.pet.rois <- function(df, tracer) {
  if (! tracer %in% c('tau', 'amyloid', 'pvc')) {
    stop('`tracer` must be "tau" or "amyloid"')
  }

  if (tracer == 'tau') {
    startpat <- 'FTP'
  } else  if (tracer == 'amyloid') {
    startpat <- 'AV45'
  } else if (tracer == 'pvc') {
    startpat <- 'FTPPVC'
  }

  all.cols <- colnames(df)
  pet.lh.cols <- all.cols[str_detect(all.cols, sprintf('%s_CTX_LH_|%s_LEFT_', startpat, startpat))]
  pet.rh.cols <- all.cols[str_detect(all.cols, sprintf('%s_CTX_RH_|%s_RIGHT_', startpat, startpat))]
  pet.lh <- df[, pet.lh.cols]
  pet.rh <- df[, pet.rh.cols]
  vol.lh <- df[, all.cols[str_detect(all.cols, 'CTX_LH_.*_VOLUME|LEFT_.*_VOLUME')]]
  vol.rh <- df[, all.cols[str_detect(all.cols, 'CTX_RH_.*_VOLUME|RIGHT_.*_VOLUME')]]
  tot.vol <- vol.lh + vol.rh
  lh.weight <- vol.lh / tot.vol
  rh.weight <- vol.rh / tot.vol
  print(length(pet.lh.cols))
  print(length(pet.rh.cols))
  weighted.pet <- (lh.weight * pet.lh) + (rh.weight * pet.rh)
  colnames(weighted.pet) <- colnames(pet.lh)
  colnames(weighted.pet) <- str_replace(colnames(weighted.pet), 'LH_|LEFT_', 'TOT_')

  return (weighted.pet)
}

amy.bl <- bilateral.pet.rois(df, 'amyloid')
tau.bl <- bilateral.pet.rois(df, 'tau')
tau.pvc.bl <- bilateral.pet.rois(df, 'pvc')

df <- df %>%
  bind_cols(amy.bl, tau.bl, tau.pvc.bl)

# === add Mattsson composites =======

cols <- colnames(df)
av45.df <- df[, str_detect(cols, 'AV45_') & str_detect(cols, 'TOT_')]
acols <- colnames(av45.df)

MattssonEarly.regions <- c('PRECUNEUS',
                           'POSTERIORCINGULATE',
                           'ISTHMUSCINGULATE',
                           'INSULA',
                           'MEDIALORBITOFRONTAL',
                           'LATERALORBITOFRONTAL')
MattssonEarly.df <- av45.df[, str_detect(acols, paste(MattssonEarly.regions, collapse='|'))]

MattssonIntermediate.regions <- c('BANKSSTS',
                                  'CAUDALMIDDLEFRONTAL',
                                  'CUNEUS',
                                  'FRONTALPOLE',
                                  'FUSIFORM',
                                  'INFERIORPARIETAL',
                                  'INFERIORTEMPORAL',
                                  'LATERALOCCIPITAL',
                                  'MIDDLETEMPORAL',
                                  'PARAHIPPOCAMPAL',
                                  'PARSOPERCULARIS',
                                  'PASORBITALIS',
                                  'PARSTRIANGULARIS',
                                  'PUTAMEN',
                                  'ROSTRALANTERIORCINGULATE',
                                  'ROSTRALMIDDLEFRONTAL',
                                  'SUPERIORFRONTAL',
                                  'SUPERIORPARIETAL',
                                  'SUPERIORTEMPORAL',
                                  'SUPRAMARGINAL')
MattssonIntermediate.df <- av45.df[, str_detect(acols, paste(MattssonIntermediate.regions, collapse='|'))]

MattssonLate.regions <- c('LINGUAL',
                          'PERICALCARINE',
                          'PARACENTRAL',
                          'PRECENTRAL',
                          'POSTCENTRAL')
MattssonLate.df <- av45.df[, str_detect(acols, paste(MattssonLate.regions, collapse='|'))]

df$MattssonEarlySUVR <- rowMeans(MattssonEarly.df)
df$MattssonIntermediateSUVR <- rowMeans(MattssonIntermediate.df)
df$MattssonLateSUVR <- rowMeans(MattssonLate.df)

# === Add Collij 2020 merged regions =======

volume.weighted.mean <- function(pet.data, vol.data, search.columns) {

  pattern <- paste(search.columns, collapse='|')
  pet.cols <- colnames(pet.data)[str_detect(colnames(pet.data), pattern)]
  vol.cols <- colnames(vol.data)[str_detect(colnames(vol.data), pattern)]

  ncols <- length(pet.cols)
  print(sprintf('%s column(s) selected for averaging:', ncols))
  print(pet.cols)

  pet <- pet.data[, pet.cols]
  vol <- vol.data[, vol.cols]

  volumes.norm <- vol / rowSums(vol)
  pet.norm <- pet * volumes.norm
  result <- rowSums(pet.norm)

  return(result)
}

pet.data <- df %>% select(matches('AV45_CTX_(LH|RH)_.*_SUVR'))
vol.data <- df %>% select(matches('CTX_(LH|RH)_.*_VOLUME'))
compare.cols <- data.frame(a=colnames(pet.data), b=colnames(vol.data))

x <- c('ANTERIORCINGULATE')
df$CollijAnteriorCingulate <- volume.weighted.mean(pet.data, vol.data, x)

x <- c('PARSOPERCULARIS', 'PARSTRIANGULARIS', 'PARSORBITALIS')
df$CollijInferiorFrontal <- volume.weighted.mean(pet.data, vol.data, x)

x <- c('CAUDALMIDDLEFRONTAL', 'ROSTRALMIDDLEFRONTAL')
df$CollijMiddleFrontal <- volume.weighted.mean(pet.data, vol.data, x)

# === Braak regions =======

braak1.regs <- c('ENTORHINAL')

braak3.regs <- c('PARAHIPPOCAMPAL',
                 'FUSIFORM',
                 'LINGUAL',
                 'AMYGDALA')

braak4.regs <- c('MIDDLETEMPORAL',
                 'CAUDALANTERIORCINGULATE',
                 'ROSTRALANTERIORCINGULATE',
                 'POSTERIORCINGULATE',
                 'ISTHMUSCINGULATE',
                 'INSULA',
                 'INFERIORTEMPORAL',
                 'TEMPORALPOLE')

braak5.regs <- c('SUPERIORFRONTAL',
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
                 'TRANSVERSETEMPORAL')

braak6.regs <- c('PERICALCARINE',
                 'POSTCENTRAL',
                 '_CUNEUS', # underscore so you don't precuneus also
                 'PRECENTRAL',
                 'PARACENTRAL')

df.tau <- df %>%
  select(matches('FTP_.*_SUVR') & ! contains('TOT_'))
df.taupvc <- df %>%
  select(matches('FTPPVC_.*_SUVR') & ! contains('TOT_'))
df.vol <- df %>%
  select(matches('_VOLUME'))

# tau w/o PVC
df$BRAAK1_TAU <- volume.weighted.mean(df.tau, df.vol, braak1.regs)
df$BRAAK3_TAU <- volume.weighted.mean(df.tau, df.vol, braak3.regs)
df$BRAAK4_TAU <- volume.weighted.mean(df.tau, df.vol, braak4.regs)
df$BRAAK5_TAU <- volume.weighted.mean(df.tau, df.vol, braak5.regs)
df$BRAAK6_TAU <- volume.weighted.mean(df.tau, df.vol, braak6.regs)

df$BRAAK34_TAU <- volume.weighted.mean(df.tau, df.vol, c(braak3.regs, braak4.regs))
df$BRAAK56_TAU <- volume.weighted.mean(df.tau, df.vol, c(braak5.regs, braak6.regs))

# tau w/ PVC
df$BRAAK1_TAUPVC <- volume.weighted.mean(df.taupvc, df.vol, braak1.regs)
df$BRAAK3_TAUPVC <- volume.weighted.mean(df.taupvc, df.vol, braak3.regs)
df$BRAAK4_TAUPVC <- volume.weighted.mean(df.taupvc, df.vol, braak4.regs)
df$BRAAK5_TAUPVC <- volume.weighted.mean(df.taupvc, df.vol, braak5.regs)
df$BRAAK6_TAUPVC <- volume.weighted.mean(df.taupvc, df.vol, braak6.regs)

df$BRAAK34_TAUPVC <- volume.weighted.mean(df.taupvc, df.vol, c(braak3.regs, braak4.regs))
df$BRAAK56_TAUPVC <- volume.weighted.mean(df.taupvc, df.vol, c(braak5.regs, braak6.regs))

# === remove NAs =========

all.cols <- colnames(df)
na.cols <- c('Age', 'Sex', 'DeltaMMSE', 'DeltaPHC_MEM', 'HasE4', 'CDRBinned')
roi.cols <- all.cols[str_detect(all.cols, '_SUVR|_VOLUME')]

df.withna <- df
df <- df %>%
  drop_na(all_of(na.cols), all_of(roi.cols))

# === PCA on PHC ======

# train PCA on the baseline data
# VSP is omitted due to data seeming wonky
pca.data <- df %>%
  select(PHC_MEM, PHC_EXF, PHC_LAN)
pca <- prcomp(pca.data, center = T, scale = T)

# get PC1 loading for baseline data
predicts.bl <- as.data.frame(predict(pca, pca.data))
df$PHC_PC1 <- predicts.bl$PC1

# get PC1 loadings for longitudinal data
predicts.long <- as.data.frame(predict(pca, df.adsp.long))
df.adsp.long$PHC_PC1 <- predicts.long$PC1

# model longitudinal change in PC1
# ADSP
df <- calc.longitudinal.change(
  baseline = df,
  longitudinal = df.adsp.long,
  variable = 'PHC_PC1',
  date.column = 'DateADSP'
)

# === ML-friendly variables =======

df$SexBinary <- ifelse(df$Sex == 'Male', 1, 0)
df$HasE4Binary <- ifelse(df$HasE4, 1, 0)

# === Add longitudinal followup metrics ======

# calc.longitudinal.change <- function(baseline, longitudinal,
#                                      variable, date.column,
#                                      id.column='RID', age.column='Age',
#                                      plot.by='CDRBinned') {}


# Reusing some code from calc.longitudinal change
# to look at the distributions of longitudinal visits
joiner <- mmse.long %>%
  select(RID, DateMMSE, MMSE) 

longmerge <- df %>%
  select(RID, DateMMSE, Age) %>%
  rename(DateMMSEBL=DateMMSE) %>%
  left_join(joiner, by='RID') %>%
  group_by(RID) %>%
  mutate(DELTA = as.numeric(difftime(DateMMSE, DateMMSEBL, units='days')) / 365.25,
         LONG.AGE = Age + DELTA) %>%
  filter(DateMMSE >= DateMMSEBL) %>%
  filter(n() >= 2) %>%
  drop_na(MMSE) %>%
  ungroup()

longstats <- longmerge %>%
  group_by(RID) %>%
  summarise(NVisits = n(),
            FollowupYears = max(DELTA)) %>%
  ungroup()

df <- df %>%
  left_join(longstats, by = c('RID'))

# === CSF markers =======

# looks like the CSF markers are sparser for this dataset
# and the ADNI1/2/GO vs. ADNI3 assays have values on very different scales

# --->  this code merges the older values (1/2/GO)
# old.csf <- upennbiomk_master %>%
#   mutate(RID = as.numeric(RID),
#          DateCSF = as_datetime(mdy(DRAWDTE))) %>%
#   select(RID, DateCSF, ABETA, TAU, PTAU)
#
# df <- left_join(df, old.csf, by='RID') %>%
#   mutate(DiffTauCSF = as.numeric(difftime(DateTau, DateCSF, units='days')) / 365.25) %>%
#   group_by(TauID) %>%
#   slice_min(abs(DiffTauCSF), with_ties = F) %>%
#   filter(abs(DiffTauAmyloid) < THRESHOLD.IMAGING.DAYS) %>%
#   ungroup() %>%
#   select(TauID, DateTau, ABETA, TAU, PTAU)
# --->

bmk10 <- upennbiomk10 %>%
  select(RID, DRAWDATE, ABETA40, ABETA42, TAU, PTAU) %>%
  rename(DateCSF=DRAWDATE,
         CSF_ABETA40=ABETA40,
         CSF_ABETA42=ABETA42,
         CSF_TAU=TAU,
         CSF_PTAU=PTAU) %>%
  mutate(DateCSF=as_datetime(mdy(DateCSF)))

bmk12 <- upennbiomk12_2020 %>%
  select(RID, EXAMDATE, AB40, ABETA, TAU, PTAU) %>%
  rename(DateCSF=EXAMDATE,
         ABETA40=AB40,
         CSF_ABETA40=AB40,
         CSF_ABETA42=ABETA,
         CSF_TAU=TAU,
         CSF_PTAU=PTAU) %>%
  mutate(DateCSF = as_datetime(ymd(DateCSF)))

new.csf <- rbind(bmk10, bmk12)

df.csf <- left_join(df, new.csf, by='RID') %>%
  mutate(DiffTauCSF = as.numeric(difftime(DateTau, DateCSF, units='days')) / 365.25) %>%
  group_by(TauID) %>%
  slice_min(abs(DiffTauCSF), with_ties = F) %>%
  filter(abs(DiffTauAmyloid) < THRESHOLD.IMAGING.DAYS) %>%
  ungroup() %>%
  mutate(CSF_AB42OVER40=CSF_ABETA42/CSF_ABETA40) %>%
  filter(! is.na(CSF_PTAU))

# === Plasma markers =======

plasma <- load.adni.table('UPENN_PLASMA', 'inputs')

plasma <- plasma %>%
  mutate(
    DatePlasma = as_datetime(ymd(EXAMDATE)),
    PLASMA_AB40 = AB40_F,
    PLASMA_AB42 = AB42_F,
    PLASMA_AB42OVER40 = AB42_AB40_F,
    PLASMA_PTAU217 = pT217_F,
    PLASMA_PTAU217OVERAB42 = pT217_AB42_F,
    PLASMA_NFL = NfL_Q
    ) %>%
  select(RID, DatePlasma, 
         starts_with('PLASMA'))

df.plasma <- left_join(df, plasma, by = 'RID') %>%
  mutate(DiffTauPlamsa = as.numeric(difftime(DateTau, DatePlasma, units='days')) / 365.25) %>%
  group_by(TauID) %>%
  slice_min(abs(DiffTauPlamsa), with_ties = F) %>%
  filter(abs(DiffTauPlamsa) < THRESHOLD.IMAGING.DAYS) %>%
  ungroup() %>%
  filter(! is.na(PLASMA_PTAU217), PLASMA_NFL != -4)

# === save ========

write.csv(df, file.path(outfolder, 'maindata.csv'), quote = F, na = '', row.names = F)
write.csv(df.csf, file.path(outfolder, 'maindata_csf.csv'), quote = F, na = '', row.names = F)
write.csv(df.plasma, file.path(outfolder, 'maindata_plasma.csv'), quote = F, na = '', row.names = F)

# === Table 1 =======

df$Race <- factor(df$Race, levels = c('White', 'Black', 'Asian', 'Other', 'Unknown'))
df$Hispanic <- as.character(df$Hispanic)

vars <- c(
  'Age', 'Sex', 'Race', 'Hispanic', 'Education',
  'HasE4', 'Centiloid', 'META_TEMPORAL_TAU', 'HIPPOCAMPUS_VOL',
  'MMSE', 'DeltaMMSE',
  'CDRSumBoxes', 'DeltaCDRSumBoxes',
  'PACC', 'DeltaPACC',
  'NVisits', 'FollowupYears'
  )

tbl1 <- CreateTableOne(vars=vars,
                       strata='CDRBinned',
                       data=df)
print(tbl1, showAllLevels=T)

# === Table 1: CSF =======

vars.csf <- c(vars,  'CSF_AB42OVER40', 'CSF_PTAU', 'CSF_TAU')

tbl1.csf <- CreateTableOne(vars=vars.csf,
                       strata='CDRBinned',
                       data=df.csf)
print(tbl1.csf, showAllLevels=T)


# === Table 1: Plasma =======

vars.plasma <- c(vars,  'PLASMA_AB42OVER40', 'PLASMA_PTAU217', 'PLASMA_NFL')

tbl1.plasma <- CreateTableOne(vars=vars.plasma,
                           strata='CDRBinned',
                           data=df.plasma)
print(tbl1.plasma, showAllLevels=T)

# === Look at multicollinearity ========

continuous <- df %>%
  select(Centiloid, META_TEMPORAL_TAU, BRAAK1_TAU, BRAAK34_TAU, BRAAK56_TAU,
         HIPPOCAMPUS_VOL, META_TEMPORAL_VOL)

cormat <- cor(continuous)

m <- lm(DeltaMMSE ~ Centiloid + META_TEMPORAL_TAU + HIPPOCAMPUS_VOL, data=df)

