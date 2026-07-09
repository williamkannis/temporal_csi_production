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
export_dir <- "prod_data"

# Growth model directories
grow_dir <- "~/Documents/Work/Everglades post-doc/Data analysis/growth curves"
stack_dir <- file.path(grow_dir,"loo_outputs_cat")
stackJ_dir <- file.path(grow_dir,"loo_outputs")
mod_dir <- file.path(grow_dir,"stan_outputs/model_out")

# Data
len_df <- readRDS(file.path(input_dir,"fslen_imputed_2026-07-09.rds"))
hyd_df <- readRDS(file.path(input_dir,"hydr_class_annual_2026-07-09.rds"))
wt_df <- read.csv(file.path(input_dir,"length_weight_parameters.csv"))
stack_list <- readRDS(file.path(stack_dir,"stack_wt_out_2026-06-22.rds"))
stackJ <- readRDS(file.path(stackJ_dir,"stack_wt_out_2026-06-22.rds"))["JORFLO"]


# Data preparation  ------------------------------------------------------------

# Create data.frame containing sampling event, date, sampling area (i.e.,
# number of traps), and interval between sampling periods
samp_df <- len_df %>% 
  left_join(hyd_df) %>% 
  mutate(group_id = as.numeric(hydroperiod)) %>% 
  group_by(site,cum,group_id) %>% 
  summarise(
    date = mean(date,na.rm=T),
    area = n_distinct(plot,throw),
    .groups = "drop"
    ) %>% 
  sample_interval()

# Check if any missing interval has sequential sampling event
samp_df %>% 
  filter(is.na(interval)) %>% 
  distinct(cum,site) %>% 
  mutate(cum = cum+1) %>% 
  inner_join(samp_df)
# No issues


### REQUIRE USERS TO HAVE CONSEQUATIVE SAMPLE COUNTER, MAKE A FUNCTION
# TO GENREATE THIS . MAYBE USE DATES AND USE SAMPLE INTERVAL FUNCTION TO DENOTE
# MISSING INTERVALS IF DESIRED

# Estimate number of days between each sampling event
interval_df <- samp_df %>%
  filter(!is.na(interval)) %>% 
  distinct(group_id,interval)



# Estimate biomass for each fish and attach sampling info
bio_df <- len_df %>% 
  left_join(wt_df) %>% 
  mutate(
    wet_wt = 10^(a + b * log10(length*c)),
    wt = wet_wt*.19
    ) %>% 
  select(
    region,
    site,
    wateryear,
    cum,
    species,
    length,
    wt
    )


# Species-specific production estimates  ---------------------------------------

# JORFLO has no group speficic growth esitmates and stacking wts are contained 
# in seperate file. Adds these to main stacking weigth list
stack_list <- c(stack_list,stackJ)

# Production input settings
sp <- names(stack_list)
cohort <- 30
growth_iter <- 10
prod_iter <- 10

prod_list <- lapply(sp, function(s){
  
  # Filter all data for one species
  samp <- samp_df
  bio <- bio_df %>% filter(species == s)
  stack <- stack_list[[s]]
  sp_dir <- file.path(mod_dir,s)
  wt <- wt_df %>% filter(species == s)
  group_type = "cat"
  pred_group <- interval_df$group_id
  pred_interval <- interval_df$interval
  
  # Create species specific length input data
  length_vec <- min(bio$length,na.rm = T):max(bio$length,na.rm = T)
  
  # JORFLO does not have group specific growth rates, use population
  if(s == "JORFLO") {
    group_type <- "mu"
    pred_group <- NULL
    samp <- samp %>% mutate(group_id = 1)
  }
  
  # Species specific growth and age at length estimation
  growth_post <- stack_predict(
    stack.df = stack,
    mod.dir = sp_dir,
    sim = growth_iter,
    summarize = F,
    sum.fun = "median",
    type = "prediction",
    group.id = group_type,
    pred.input = length_vec,
    create.input = T,
    pred.group = pred_group,
    pred.interval = pred_interval,
    stack = T,
    input.var = "length",
    output.var = c("interval_growth","age"),
    wt.df = wt,
    dry.wt = .19,
    parallel = T,
    mc.cores = 10
  )
  
  # Species specific production
  production(
    method = "igr",
    sample = samp,
    biomass = bio,
    growth = growth_post, 
    class.type = "size",
    class.size = 1,
    bio.boot = T, 
    growth.boot = T,
    iter = prod_iter,
    return.raw = F,
    parallel = T,
    mc.cores = 10
  ) %>% 
    mutate(species = s)
})
prod_df <- bind_rows(prod_list)


# Export  ----------------------------------------------------------------------
prod_file <- paste0("fsprod_igr_",Sys.Date(),".rds")
saveRDS(prod_df,file.path(export_dir,prod_file))



