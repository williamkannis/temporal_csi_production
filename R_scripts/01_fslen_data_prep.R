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

# directories
data_dir <- file.path(
  "~/Documents/Work/Everglades post-doc/Data analysis/Data cleaning",
  "cleaned_data")
input_dir <- "input_data"

# Data
len_df <- readRDS(file.path(data_dir,"fslen_cleaned_2026-11-03.rds"))


# Data preparation  ------------------------------------------------------------

# For this manuscript we will only focus on the 6 most abundant species:
sp <- c(
  "FUNCHR",
  "HETFOR",
  "GAMHOL",
  "JORFLO",
  "LUCGOO",
  "POELAT",
  "NOFISH"  # include sites with no fish as true zeros
)
len_filtered <- len_df %>% 
  filter(
    species %in% sp,
    wateryear > 1996  
    )


# Length summary  --------------------------------------------------------------

# Proportion of data made up by big 6
len_df %>% 
  filter(wateryear > 1996) %>% 
  mutate(
    focal = case_when(
      species %in% sp ~ T,
      T ~ F
    )
  ) %>% 
  summarize(
    tot = n(),
    n_focal= sum(focal),
    prop_focal = 100*n_focal/tot
  )


len_df %>% 
  filter(
    wateryear > 1996,
    !species %in% c(sp,"NOFISH")
    ) %>% 
  distinct(site,cum) %>% 
  nrow()

# Percent of nonfocal species at sites
prop_focal <-len_df %>% 
  filter(wateryear > 1996) %>% 
  mutate(
    type = case_when(
      species %in% sp[-7] ~ "focal",
      species == "NOFISH" ~ "nofish",
      T ~ "rare"
    )
  ) %>% 
  group_by(site,cum) %>% 
  summarise(
    n_focal = sum(type=="focal"),
    n_rare = sum(type == "rare"),
    n = n(),
    prop_rare = n_rare/n
  )

prop_focal %>% 
  filter(prop_rare == 1) %>% 
  nrow()
# No rare species made up an entire sites density. All missing sites after
#  filtering are truly missing, rather than just missing the focal species

median(prop_focal$prop_rare)
1-sum(prop_focal$prop_rare <.5)/nrow(prop_focal)
# The median sample has  ~7% of its density made up by rare species
# and only ~3% of samples have rare specis making up half its density

# Missing lengths
len_filtered %>% 
  summarize(
    tot = n(),
    n_missing= length(length[is.na(length)]),
    prop_missing = 100*n_missing/tot
  )

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

# Find median values for each species and season
miss_meds <- len_filtered %>% 
  filter(species != "NOFISH") %>% 
  group_by(species,cum) %>% 
  summarise(
    length_median = median(length,na.rm = T)
  )

# Replace missing values with seasonal medians
len_impute <- len_filtered %>% 
  left_join(miss_meds) %>% 
  mutate(
    impute = case_when(
      is.na(length) ~ T,
      T ~ F
    ),
    length = case_when(
      is.na(length) ~ length_median,
      T ~ length
    )
  ) %>% 
  select(-length_median)

# Did all NAs get imputed?
len_impute %>% 
  filter(is.na(length))
summary(len_impute)

# check that number of rows is the same
nrow(len_impute) == nrow(len_filtered)


# Export  ----------------------------------------------------------------------
file_name <- paste0("fslen_imputed_",Sys.Date(),".rds")
saveRDS(len_impute,file.path(input_dir,file_name))

