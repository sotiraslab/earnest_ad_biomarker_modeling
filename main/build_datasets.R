# === imports ======

sh <- suppressPackageStartupMessages

sh(library(ADNIMERGE))
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

# === Add cognitive measures =========

adsp <- load.adni.table('ADSP_PHC', 'inputs')

adsp <- adsp %>%
  mutate(DateADSP = as_datetime(ymd(EXAMDATE))) %>%
  select(RID, DateADSP, PHC_MEM, PHC_EXF, PHC_LAN, PHC_VSP)
adsp$PHC_GLOBAL <- rowMeans(adsp[, c('PHC_MEM', 'PHC_EXF', 'PHC_LAN', 'PHC_VSP')])

#save for calculating longitudinal change
df.adsp.long <- left_join(df, adsp, by='RID') %>%
  mutate(DiffImagingADSP = as.numeric(difftime(MeanImagingDate, DateADSP, units = 'days')))

df.adsp <- df.adsp.long %>%
  group_by(RID) %>%
  slice_min(order_by=abs(DiffImagingADSP), with_ties = F) %>%
  ungroup() %>%
  filter(abs(DiffImagingADSP) < THRESHOLD.COGNITIVE.DAYS)

df <- df.adsp

# === Compute longitudinal change in PHC =======

joiner <- df.adsp.long %>%
  select(RID, DateADSP, PHC_GLOBAL, PHC_MEM, PHC_EXF, PHC_LAN, PHC_VSP) %>%
  filter(! is.na(PHC_GLOBAL))

long.data <- df %>%
  select(RID, DateADSP, Age, CDRBinned) %>%
  rename(DateADSP.BL=DateADSP) %>%
  left_join(joiner, by='RID') %>%
  group_by(RID) %>%
  mutate(DeltaADSPDate = as.numeric(difftime(DateADSP, DateADSP.BL, units='days')) / 365.25,
         Long.Age = Age + DeltaADSPDate) %>% 
  filter(DateADSP >= DateADSP.BL) %>%
  filter(n() >= 2) %>%
  ungroup()

# longitudinal modelling
m <- lmer(PHC_GLOBAL ~ DeltaADSPDate + (1+DeltaADSPDate|RID), data=long.data)
long.data$PHC_GLOBAL.LMER.Predict <- predict(m, long.data)

ggplot(long.data, aes(x=Long.Age, y=PHC_GLOBAL)) +
  geom_point(aes(color=CDRBinned), alpha = .7) + 
  geom_line(aes(y=PHC_GLOBAL.LMER.Predict, group=RID, color=CDRBinned), alpha= .7)

coefs <- coef(m)$RID %>%
  select(DeltaADSPDate) %>%
  dplyr::rename(DeltaADSP=DeltaADSPDate) %>%
  rownames_to_column(var="RID") %>%
  mutate(RID=as.numeric(RID))

df <- left_join(df, coefs, by='RID')

# scale up deltaADSP for training purposes
# found to be getting wonky SVM results without this step
# factor of 10 puts it on similar scale as cross-sectional value
df <- df %>%
  mutate(DeltaADSP = DeltaADSP * 10)

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
na.cols <- c('Age', 'Sex', 'PHC_GLOBAL', 'HasE4', 'CDRBinned')
roi.cols <- all.cols[str_detect(all.cols, '_SUVR|_VOLUME')]

df.withna <- df
df <- df %>%
  drop_na(all_of(na.cols), all_of(roi.cols))

# === ML-friendly variables =======

df$SexBinary <- ifelse(df$Sex == 'Male', 1, 0)
df$HasE4Binary <- ifelse(df$HasE4, 1, 0)

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

# === save ========

df.long <- df %>%
  filter(! is.na(DeltaADSP))

write.csv(df, file.path(outfolder, 'maindata.csv'), quote = F, na = '', row.names = F)
write.csv(df.long, file.path(outfolder, 'maindata_long.csv'), quote = F, na = '', row.names = F)
write.csv(df.csf, file.path(outfolder, 'maindata_csf.csv'), quote = F, na = '', row.names = F)

# === Table 1 =======

vars <- c('Age', 'Sex', 'HasE4', 'Centiloid', 'PHC_GLOBAL')

tbl1 <- CreateTableOne(vars=vars,
                       strata='CDRBinned',
                       data=df)
print(tbl1, showAllLevels=T)

# === Table 1: longitudinal =======

tbl1.long <- CreateTableOne(vars=vars,
                       strata='CDRBinned',
                       data=df.long)
print(tbl1.long showAllLevels=T)

# === Table 1: CSF =======

vars <- c('Age', 'Sex', 'HasE4', 'Centiloid', 'PHC_GLOBAL')

tbl1.csf <- CreateTableOne(vars=vars,
                       strata='CDRBinned',
                       data=df.csf)
print(tbl1.csf, showAllLevels=T)
