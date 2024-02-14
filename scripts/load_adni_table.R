
sh <- suppressPackageStartupMessages

sh(library(this.path))

DEFAULT.DIRECTORY <- normalizePath(file.path(this.dir(), '..', 'inputs'))

load.adni.table <- function(table, directory = NULL) {
  if (is.null(directory)) directory <- DEFAULT.DIRECTORY
  table <- tolower(table)
  files <- list.files(directory)
  matches <- files[grepl(table, files, ignore.case = T, perl = T)]
  if (length(matches) == 0) {
    stop(sprintf('Cannot find table matching "%s" in "%s"', table, directory))
  } else if (length(matches) > 1) {
    stop(sprintf('More than one table matching "%s" in "%s"', table, directory))
  } 
  path <- file.path(directory, matches[1])
  df <- read.csv(path)
  return (df)
}
