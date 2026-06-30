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

# Packages under development (switch to github for publication)
devtools::load_all("~/Documents/work/R packages/growthstack")
devtools::load_all("~/Documents/work/R packages/secProd")

# Production directories
input_dir <- "input_data"
data_dir <- file.path(
  "~/Documents/Work/Everglades post-doc/Data analysis/Data cleaning",
  "cleaned_data"
)

# Growth model directories
grow_dir <- "~/Documents/Work/Everglades post-doc/Data analysis/growth curves"
stack_dir <- file.path(grow_dir,"loo_outputs_cat")
mod_dir <- file.path(grow_dir,"stan_outputs/model_out")

# Data
len_df <- readRDS(file.path(data_dir,"fslen_cleaned_2026-06-30.rds"))
hyd_df <- readRDS(file.path(input_dir,"hydr_class_annual_2026-06-25.rds"))
wt_df <- read.csv(file.path(input_dir,"length_weight_parameters.csv"))
stack_list <- readRDS(file.path(stack_dir,"stack_wt_out_2026-06-22.rds"))


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

# Estiamte number of days between each sampling event
interval_df <- sample_interval(samp_df) %>% 
  filter(!is.na(interval)) %>% 
  left_join(samp_df)


# Estimate biomass for each fish and attach sampling info
bio_df <- len_df %>% 
  left_join(hyd_df) %>% 
  left_join(samp_df %>% select(-date)) %>% 
  left_join(wt_df) %>% 
  left_join(interval_df) %>% 
  mutate(
    wet_wt = 10^(a + b * log10(length*c)),
    dry_wt = wet_wt*.19
    ) %>% 
  select(region,site,wateryear,cum,interval,group_id,area, species,length,dry_wt)
  


# Species-specific production estimates  ---------------------------------------
sp <- names(stack_list)[[4]]
cohort <- 30

growth_list <- lapply(sp, function(sp){
  
  # Filter all data for one species
  bio <- bio_df %>% filter(species == sp)
  stack <- stack_list[[sp]]
  sp_dir <- file.path(mod_dir,sp)
  wt <- wt_df %>% filter(species == sp)
  
  # Create species specific length input data
  length_vec <- min(bio$length,na.rm = T):max(bio$length,na.rm = T)
  
  # Species specific growth and age at length estimation
  growth_post <- stack_predict(
    stack.df = stack,
    mod.dir = sp_dir,
    sim = 100,
    summarize = F,
    sum.fun = "median",
    type = "prediction",
    group.id = "cat",
    pred.input = length_vec,
    create.input = T,
    pred.group = interval_df$group_id,
    pred.interval = interval_df$interval,
    stack = T,
    input.var = "length",
    output.var = c("interval_growth","age"),
    wt.df = wt,
    dry.wt = .19
  )
  
  # Species specific production
  prod <- production_igr(
    growth = as.data.frame(growth_post[,,1]),
    biomass = bio,
    size.class = cohort
    # n.sim = 10,
    # return.raw=F,
    # parallel=F,
    # mc.cores=NULL
  )
  
})


# Export  ----------------------------------------------------------------------




