
# --- Imports -----
sh <- suppressPackageStartupMessages

sh(library(stringr))
sh(library(this.path))
sh(library(tidyverse))

# --- Working directory -----
setwd(this.dir())

# --- Load data -----

df <- read.csv('datasets/maindata.csv')
mmse <- read.csv('datasets/mmse_lmem_fit_data.csv')

mmse <- mmse %>%
  filter(ID %in% df$RID) %>%
  group_by(ID) %>%
  filter(n() >=2) %>%
  mutate(AGEBL = first(LONG.AGE)) %>%
  ungroup() %>%
  arrange(AGEBL, ID) %>%
  mutate(CogStatus = ifelse(PLOTBY == '0.0', 'CU', 'CI'),
         Age = LONG.AGE,
         ID.FACTOR = fct_inorder(as.character(ID))) %>%
  group_by(ID.FACTOR) %>%
  mutate(Subject = cur_group_id()) %>%
  ungroup()

ggplot(mmse, aes(x=Age, y=Subject, group=Subject,
                 color=CogStatus, fill=CogStatus)) +
  geom_line() +
  geom_point()
  