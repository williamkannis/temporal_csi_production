#-------------------------------------------------------------------------------
#
#  Hydroperiod Classification  
#
#-------------------------------------------------------------------------------

# AUTHOR: William Annis

# CREATED: Feb 2, 2026

# DESCRIPTION: This script inputs the hydrological site data from the cleaned
# data folder and uses these data to classify sampling events (site+years) into
# 3 classifications based on annual inundated days. Data are available at period
# level for plots and throws, and these need to be aggregated to the year and
# site level. 


# Housekeeping  ----------------------------------------------------------------
rm(list = ls())

# Load in packages
library(dplyr)

# Directories
data_dir <- paste0(
  "~/Documents/Work/Everglades post-doc/",
  "Data analysis/Data cleaning/cleaned_data"
)
input_dir <- "input_data"


# Data
phy_df <- readRDS(file.path(data_dir,"phys_cleaned_2026-02-25.rds"))


# Create hydroperiod  ----------------------------------------------------------
# Create Site level hydroperiod classification for each period. For growth curve
# input data, we will use period 4 as the data were collected in October. We
# will create hydroperiods for each period for the use with the main length data
# set

hyd_df <- phy_df %>% 
  
  # Summarize hydroperiods from throw to site level for each period
  group_by(cum,year,period,wateryear,waterperiod,region,site,) %>% 
  summarize(hydro_annual = mean(wet_sum_365day,na.rm=T)) %>% 
  ungroup() %>% 
  
  # Create hydroperiod groupings    
  mutate(hydroperiod = case_when(
      hydro_annual > 360 ~ "long",
      hydro_annual >=320 & hydro_annual <= 360 ~ "intermediate",
      hydro_annual < 320 ~ "short",
      T~NA
    ),
    hydroperiod = factor(hydroperiod,levels=c("short","intermediate","long"))) 

# Are the number of site/year/periods the same
phy_df %>% 
  distinct(cum,region,site) %>% 
  nrow() == 
  hyd_df %>% 
  distinct(cum,region,site) %>% 
  nrow()

# Export data  -----------------------------------------------------------------

# Final data check
summary(hyd_df)
# NAs are from unsampled sites, should not be an issue

# Export
saveRDS(
  hyd_df,
  file.path(input_dir,paste0("hydr_class_annual_",Sys.Date(),".rds")
    )
  )

