#-------------------------------------------------------------------------------
#
#  Instantaneous growth rate production estimation 
#
#-------------------------------------------------------------------------------

# AUTHOR: William K. Annis

# CREATED: June 25, 2026

# DESCRIPTION: 

## !!  RERUN PRODUCTION WITH NEW LEADED GROUP IDS  !! ##


# Housekeeping  ----------------------------------------------------------------
rm(list = ls())

# Load in packages
library(dplyr)
library(tidyr)

# Packages under development (switch to github for publication)
devtools::load_all("~/Documents/work/R packages/growthstack")
devtools::load_all("~/Documents/work/R packages/secProd")

# Production directories
input_dir <- "input_data"
export_dir <- "prod_data"

# Growth model directories
grow_dir <- 
  "~/Documents/Work/Everglades post-doc/Data analysis/growth curves/outputs"
stack_dir <- file.path(grow_dir,"loo_outputs")
mod_dir <- file.path(grow_dir,"stan_outputs")

# Data
len_df <- readRDS(file.path(input_dir,"fslen_imputed_2026-07-09.rds"))
hyd_df <- readRDS(file.path(input_dir,"hydr_class_annual_2026-07-09.rds"))
wt_df <- read.csv(file.path(input_dir,"length_weight_parameters.csv"))
stack_list <- readRDS(file.path(stack_dir,"_cat-stack_wt_out.rds"))
stackJ <- readRDS(file.path(stack_dir,"stack_wt_out.rds"))["JORFLO"]


# Data preparation  ------------------------------------------------------------

# Create data.frame containing sampling event, date, sampling area (i.e.,
# number of traps), and interval between sampling periods
samp_df <- len_df %>% 

  # Estimate total sampling area and average sampling date
  group_by(site,cum) %>% 
  summarise(
    date = mean(date,na.rm=T),
    area = n_distinct(plot,throw),
    .groups = "drop"
    ) %>% 
  
  # Create group id based on hydroperiod at end of sampling interval
  left_join(
    hyd_df,
    by = join_by(site, cum)
    ) %>% 
  group_by(site) %>% 
  mutate(group_id = lead(as.numeric(hydroperiod),order_by = cum)) %>% 
  ungroup() %>% 
  
  # Estimate sampling interval
  select(site,cum,date,area,group_id) %>% 
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
  left_join(
    wt_df,
    by = join_by(species)
    ) %>% 
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
growth_iter <- 10000
prod_iter <- 10000

prod_list2 <- lapply(sp, function(s){
  
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


# Species-specific Interval density and mean size  -----------------------------

den_size_df <- bio_df %>% 
  filter(species != "NOFISH") %>% 
  
  # Create new data for all species (total response)
  mutate(species = "all") %>% 
  
  # add back in species specific data
  bind_rows(bio_df %>% filter(species != "NOFISH")) %>% 
  
  # Introduce missing species using NA length and wt. No site that contains that
  # species has an NA for length or WT
  right_join(
    samp_df %>% crossing(species = c(sp,"all")),
    by = join_by(species,site, cum)
  ) %>% 
  mutate(
    length = case_when(
      is.na(length)~ 0,
      T ~ length
    ),
    wt = case_when(
      is.na(wt)~ 0,
      T ~ wt
    )
  ) %>% 
  
  # Estimate density and sum of size at each sampling event
  group_by(site,cum,species,area,interval) %>% 
  summarise(
    n = sum(length>1),
    length_sum = sum(length),
    wt_sum = sum(wt),
    .groups = "drop"
  ) %>% 
  
  # Estimate average density and size across sampling interval
  group_by(site,species) %>% 
  mutate(
    density = n/area,
    across(
      c(n, density, wt_sum, length_sum),
      ~ case_when(
        cum + 1 != lead(cum, order_by = cum) ~ NA,
        TRUE ~ lead(.x, order_by = cum)
      ),
      .names = "lead_{.col}"
    ),
    interval_n = n+lead_n,
    interval_density = (density+lead_density)/2,
    interval_mean_length = case_when(
      interval_n != 0 ~ (length_sum+lead_length_sum)/interval_n,
      interval_n == 0 ~ 0
    ),
    interval_mean_wt = case_when(
      interval_n != 0 ~ (wt_sum+lead_wt_sum)/interval_n,
      interval_n == 0 ~ 0
    )
  ) %>% 
  
  # Prepare output
  ungroup() %>% 
  select(
    site,
    cum,
    species,
    interval_density,
    interval_mean_length,
    interval_mean_wt
  )
  

# Export  ----------------------------------------------------------------------

# sample data
samp_file <- paste0("fs_sample_info_",Sys.Date(),".rds")
saveRDS(samp_df,file.path(export_dir,samp_file))

# Production
prod_file <- paste0("fsprod_igr_",Sys.Date(),".rds")
saveRDS(prod_df,file.path(export_dir,prod_file))

# Interval density and size
den_size_file <- 'fs_densize.rds'
saveRDS(den_size_df,file.path(export_dir,den_size_file))


