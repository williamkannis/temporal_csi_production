#-------------------------------------------------------------------------------
#
#  Instantaneous growth rate production estimation 
#
#-------------------------------------------------------------------------------

# AUTHOR: William Annis

# CREATED: June 25, 2026

# DESCRIPTION: 


# Housekeeping  ----------------------------------------------------------------
rm(list = ls())

# Load in packages
library(dplyr)

# Packages under development (switch to github for publication)
devtools::load_all("~/Documents/work/R packages/growthstack")
devtools::load_all("~/Documents/work/R packages/secProd")

# Directories
data_dir <- paste0(
  "~/Documents/Work/Everglades post-doc/",
  "Data analysis/Data cleaning/cleaned_data"
)
input_dir <- "input_data"

# Data
len_df <- readRDS(file.path(data_dir,"fslen_cleaned_2026-02-25.rds"))
hyd_df <- readRDS(file.path(input_dir,"hydr_class_annual_2026-06-25.rds"))
wt_df <- read.csv(file.path(grow_dir,"length_weight_parameters.csv"))


# Data preparation  ------------------------------------------------------------

# Create data.frame containing sampling event, date, and sampling area (i.e.,
# number of traps)
samp_df <- len_df %>% 
  left_join(hyd_df) %>% 
  mutate(group_id = as.numeric(hydroperiod)) %>% 
  group_by(site,cum,group_id) %>% 
  summarise(
    date = mean(date,na.rm=T),
    area = n_distinct(plot,throw)
    ) 

sampleInterval(samp_df)

# Estimate biomass for each fish
bio_df <- lend_df %>% 
  left_join(wt_df) %>% 
  mutate(
    wet_wt = 10^(a + b * log10(length*c)),
    dry_wt = wet_wt*.19
    ) %>% 
  select(region,site,wateryear,cum, species,length,dry_wt)
  

growth_intervalInput(samp_df,seq(0:60))

# Species-specific production estimates  ---------------------------------------
sp <- unique(bio_df$species)

growth_list <- lapply(sp, function(sp){
  
  # Filter all data for one species
  samp <- samp_df %>% filter(species == sp)
  bio <- bio_df %>% filter(species == sp)
  stack <- stack_df %>% filter(species == sp)
  sp_dir
  wt <- wt_df %>% filter(species == sp)
  
  # Create species specific age classes
  age_classes
  
  # Species specific growth estimation
  growth <- stack_predict(
    stack.df = staci,
    mod.dir = sp_dir,
    sim = 100,
    summarize = F,
    sum.fun = median,
    type = "predict",
    group.id = "cat",
    pred.input = age_classes,
    create.input = t,
    pred.group = samp$hydroperiod,
    pred.interval = samp$interval,
    stack = F,
    input.var = "length",
    output.var = "interval_growth",
    wt.df = wt,
    dry.wt = .19
  )
  
  # Species specific production
  igr_prod()
  
})


# Export  ----------------------------------------------------------------------




