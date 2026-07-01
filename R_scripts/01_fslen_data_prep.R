#-------------------------------------------------------------------------------
#
#  Fish length data preperation
#
#-------------------------------------------------------------------------------

# AUTHOR: William K. Annis

# CREATED: June 29, 2026

# DESCRIPTION: Prepares fish length data collected for 30 years at 24 site in 
# the Florida Everglades for secondary production estimation. Data is filtered
# to only include the study's focal species and missing lengths are imputed 
# (> 1% of filtered data)


# Housekeeping  ----------------------------------------------------------------
rm(list = ls())

# Load in packages
library(dplyr)
library(stringr)
devtools::load_all("~/Documents/work/R packages/secProd")

# directories
data_dir <- file.path(
  "~/Documents/Work/Everglades post-doc/Data analysis/Data cleaning",
  "cleaned_data")
input_dir <- "input_data"

# Data
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
set.seed(999)

# First impute missing lengths for site and sampling period combinations with
# sufficient individuals of a species
len_tier1 <- group_impute(
  data = len_filtered,
  column = length,
  by = c(species,site,cum),
  limit = 40,
  impute.tag = "tier1") 

# Next impute missing lengths using entire sampling periods for site/periods  
# with a insufficient number of individuals
len_tier2 <- group_impute(
  data = len_tier1 %>% filter(is.na(impute)),
  column = length,
  by = c(species,cum),
  limit = 40,
  impute.tag = "tier2")

# Finally, impute missing lengths using data for the entire species if the 
# sampling period has insufficient sample sizes
len_tier3 <- group_impute(
  data = len_tier2 %>% filter(is.na(impute)),
  column = length,
  by = c(species),
  limit = 40,
  impute.tag = "tier3")

# Consolidate the imputation tiers  to one data frame
len_impute <- bind_rows(
  len_tier1 %>% filter(!is.na(impute)),
  len_tier2 %>% filter(!is.na(impute)),
  len_tier3 %>% filter(!is.na(impute))
) %>% 
  mutate(old_length = NA)

# Replace missing data with imputed data
len_final <- len_filtered %>% 
  anti_join(
    len_impute,
    by = join_by(
      cum,
      region,
      site,
      plot,
      throw,
      species,
      length == old_length
    )
  ) %>% 
  bind_rows(len_impute) %>% 
  select(-old_length)

# check that number of rows is the same
nrow(len_final) == nrow(len_filtered)

# Imputation summary
len_final %>% 
  filter(!is.na(impute)) %>% 
  group_by(species,impute) %>% 
  summarize(n = n())


# Export  ----------------------------------------------------------------------
file_name <- paste0("fslen_imputed_",Sys.Date(),".rds")
saveRDS(len_final,file.path(input_dir,file_name))
