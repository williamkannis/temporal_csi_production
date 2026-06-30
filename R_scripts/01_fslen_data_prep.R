#-------------------------------------------------------------------------------
#
#  Instantaneous growth rate production estimation 
#
#-------------------------------------------------------------------------------

# AUTHOR: William K. Annis

# CREATED: June 25, 2026

# DESCRIPTION: 


# Housekeeping  ----------------------------------------------------------------
rm(list = ls())

# Load in packages
library(dplyr)
library(stringr)

# directories
data_dir <- file.path(
  "~/Documents/Work/Everglades post-doc/Data analysis/Data cleaning",
  "cleaned_data"
)

readRDS_newest <- function(base.name,dir) {
  files <- list.files(dir,base.name)
  max_date <-max(as.Date(str_extract(files, "\\d{4}-\\d{2}-\\d{2}")))
  file.name <- paste0(
    dir,
    "/",
    base.name,
    "_",
    max_date,
    ".rds")
  print(paste0("Date = ", max_date))
  readRDS(file.name)
}

# Data

len_df <- readRDS_newest(data_dir,"fslen_cleaned")
len_df <- readRDS(file.path(data_dir,"fslen_cleaned_2026-06-30.rds"))

# Filter data  -----------------------------------------------------------------

# For this manuscript we will only focus on the 6 most abundant species:
sp <- c(
  "FUNCHR",
  "HETFOR",
  "GAMHOL",
  "JORFLO",
  "LUCGOO",
  "POELAT",
  "NOFISH"  # inlcude sites with no fish as true zeros
  )
len_filtered <- len_df %>% 
  filter(species %in% sp)

# How many fish are missing lengths per species
len_filtered %>% 
  group_by(species) %>% 
  summarize(
    tot = n(),
    n_missing= length(length[is.na(length)]),
    prop_missing = 100*n_missing/tot
    )
# less than 1 percent of the fish per each species is missing lengths
# these will be imputed below


# Length imputation  -----------------------------------------------------------


