# This script proceeds stepwise and instructs the 
# user to conduct additional steps in order to generate the
# final ICV table.  This table is provided in the GitHub for
# this project to avoid image downloading & preprocessing steps.

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

# input
PATH.DF <- '../../derivatives/adni_base_table.csv'

# intermediate files
PATH.PTIDS <- 'PTIDS.txt'
PATH.DOWNLOAD.IMAGES <- 'DOWNLOAD_IDS.txt'
PATH.MRISEARCH <- 'MRISEARCH.csv'
PATH.MASKSIZES <- 'MASKSIZES.csv'
PATH.DOWNLOADLIST <- 'DOWNLOADLIST.txt'

# output
PATH.ICV.TABLE <- '../../derivatives/adni_icvs.csv'


# === Step 1: Produce search for MRI ======

# get a list of MRIs for each subject
# must run a search on IDA for this

# ICV is deselected here to prevent problems
# regenerating the column
df <- read.csv(PATH.DF) %>%
  select(-ICV)

inv <- inventory %>%
  select(RID, PTID) %>%
  filter(!duplicated(PTID))

df <- left_join(df, inv, by='RID')

write.table(df$PTID, PATH.PTIDS, sep=',', row.names = F, quote = F, eol=',', col.names = F)

# === check for MRI search results ======

if (! file.exists(PATH.MRISEARCH)) {
  msg <- sprintf('"%s" not found.  Use "%s" to get a list of MRIs for all subjects.', PATH.MRISEARCH, PATH.PTIDS)
  stop(msg)
}

# === Step 2: Filter MRI results ======

# first look based on excluding unwanted scans
exclude_search <- read.csv(PATH.MRISEARCH) %>%
  filter(! str_detect(tolower(Description), 'localizer'),
         ! str_detect(tolower(Description), 'calibration'),
         ! str_detect(tolower(Description), 'fgre'),
         ! str_detect(tolower(Description), 'scount'),
         ! str_detect(tolower(Description), 'loc'),
         ! str_detect(tolower(Description), 'scout'),
         ! str_detect(tolower(Description), 'take off auto send'),
         ! str_detect(tolower(Description), 'field mapping'))

# then look based on wanted scans
inclusions <- 'mprage|t1|mt1|spgr|mpr|mp-rage|mp rage'
include_search <- read.csv(PATH.MRISEARCH) %>%
  filter(str_detect(tolower(Description), inclusions))

# this can be used to check what is the difference
differences = exclude_search[! exclude_search$Image.ID %in% include_search$Image.ID, ]

# include search seems to be close enough
# some additional exclusionary filtering
search <- include_search %>%
  filter(! str_detect(tolower(Description), '_nd'),
         ! str_detect(tolower(Description), 'scout'),
         ! str_detect(tolower(Description), 'mask')) %>%
  rename(PTID=Subject.ID,
         ImageID=Image.ID,
         DateMRI=Study.Date,
         DescriptionMRI=Description,
         TypeMRI=Type) %>%
  select(PTID, DateMRI, DescriptionMRI, TypeMRI, ImageID) %>%
  mutate(DateMRI=mdy(DateMRI))

# === Step 3: Find closest scan for each visit ========

# the selection of scans is done by picking the one with the 
# highest ImageID, which should roughly correspond to the most recently
# acquired and most processed

# note that in ADNI1, the processed image steps are NOT necessarily in increasing
# order.  So for earlier ADNI data, this process would not be recommended.
# But seems to be okay for mosly ADNI2/3 data

df.merged <- left_join(df, search, by='PTID', relationship="one-to-many") %>%
  mutate(DiffToMRI = as.numeric(difftime(MeanImagingDate, DateMRI, units='days'))) %>%
  group_by(TauID) %>%
  slice_min(abs(DiffToMRI)) %>%
  slice_max(ImageID) %>%
  ungroup()

missing <- df[!df$TauID %in% df.merged$TauID, ]
print(sprintf('No ICV images found for %s', paste(missing$TauID, collapse=',')))

# === Step 4: Save images for download =======

# use the image list to search on ICV and download
write.table(df.merged$ImageID, PATH.DOWNLOAD.IMAGES, sep=',', row.names=F,
            col.names = F, eol=',', quote=F, na='')

msg <- sprintf('A list of images to download has been generated (%s)',
               PATH.DOWNLOAD.IMAGES)

msg <- sprintf("Download the images in '%s' in NII format and generate ICV measurements for them using DLICV.", PATH.DOWNLOAD.IMAGES)
print(msg)

# === check for ICV results ======

if (! file.exists(PATH.MASKSIZES)) {
  msg <- sprintf('"%s" not found.  Process the images and then use mean_intensity_tool.py to generate mask sizes for DLICV outputs.', PATH.MASKSIZES)
  stop(msg)
}

if (! file.exists(PATH.DOWNLOADLIST)) {
  msg <- sprintf('"%s" not found.  In image download folder, run `find . -name "ADNI*.nii.gz" > %s` and copy the output to this folder', PATH.DOWNLOADLIST)
  stop(msg)
}

# === Step 5: Link ICVS =======

# add imageid column back to df
tmp <- select(df.merged, TauID, DateMRI, DescriptionMRI, ImageID, DiffToMRI) 
df.final <- left_join(df, tmp, by='TauID')

icv <- read.csv(PATH.MASKSIZES) %>%
  mutate(start = str_locate(Path, '\\d{3}_S_\\d{4}')[, 1],
         end = str_locate(Path, 'S\\d+')[, 2],
         Linker = substr(Path, start, end)) %>%
  select(-start, -end)

download.images <- read.table(PATH.DOWNLOADLIST, sep = ' ', header = F, col.names = c('Path'))  %>%
  mutate(start = str_locate(Path, '\\d{3}_S_\\d{4}')[, 1],
         end = str_locate(Path, 'S\\d+')[, 2],
         Linker = substr(Path, start, end),
         ImageID = str_extract(Path, 'I\\d+')) %>%
  select(-start, -end, -Path)

icv <- left_join(icv, download.images, by='Linker') %>%
  mutate(ImageID = ifelse(is.na(ImageID),
                          str_extract(Path, 'I\\d+'),
                          ImageID),
         ICV = NumberVoxels * VoxelVolume)

# not sure why there are more than one image for
# some downloads
link.icv <- select(icv, ImageID, ICV) %>%
  filter(! duplicated(ImageID))

df.final <- df.final %>%
  mutate(ImageID = ifelse(! is.na(ImageID),
                          paste('I', ImageID, sep=''),
                          ImageID)) %>%
  left_join(link.icv, by='ImageID')

# === Step 6: Impute missing ICVS =======

male.icv <- mean(df.final[!is.na(df.final$ICV) & df.final$Gender == 'Male', 'ICV'])
female.icv <- mean(df.final[!is.na(df.final$ICV) & df.final$Gender == 'Female', 'ICV'])

df.final <- df.final %>%
  mutate(ICV = ifelse(is.na(ICV) & Gender == 'Male',
                      male.icv,
                      ICV),
         ICV = ifelse(is.na(ICV) & Gender == 'Female',
                      female.icv,
                      ICV))

print('Everyone has an ICV value?')
print(sum(is.na(df$ICV)) == 0)

# === Step 7: Save =======

df.final <- df.final %>%
  select(TauID, ICV)

write.csv(df.final, PATH.ICV.TABLE, na = '', quote = F, row.names = F)