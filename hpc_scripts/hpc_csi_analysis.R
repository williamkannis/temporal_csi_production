#-------------------------------------------------------------------------------
#
#   Run stan models using HPC cluster
#
#-------------------------------------------------------------------------------

# AUTHOR: William K. Annis

# CREATED: July 13, 2026

# DESCRIPTION: 


# Housekeeping  ----------------------------------------------------------------
rm(list = ls()) 

message("Start model run script")

# Load in packages
library(dplyr)
library(purrr)
library(abind)
library(cmdrstan)

# directories
mod_dir <- "stan_scripts"
input_dir <- "stan_inputs"
out_dir <- "stan_outputs"

# Data
samp_df <- 
  readRDS(file.path(input_dir,""))
phy_reg_year <- 
  readRDS(file.path(input_dir,"phys_regionyear_predictors_2026-08-17.rds"))

message("Data files loaded in")

# Species and response selection  ----------------------------------------------

# Shell argument for selecting species and response of choice
arg <- 1
arg <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
all_sp <- unique(samp_df$species)
all_rs <- c("production_mean","ptob","biomass")
cb <- expand.grid(all_sp,all_rs)

# Select species and response of choice
sp <- cb[arg,1]
rs <- cb[arg,2]


# Model input lists  -----------------------------------------------------------

## Data filtering  ##

# Create column y for response of choice
samp_df$y <- samp_df[[rs]]

# Extract data for selected species
site_df <- samp_df %>% 
  filter(species == sp) %>% 
  
  # Create site and year id for random effects
  arrange(region,site,cum) %>% 
  mutate(site_id = cur_group_id(),.by = c(region,site)) %>% 
  mutate(year_id = cur_group_id(),.by = c(wateryear)) %>% 
  mutate(reg_id = cur_group_id(),.by = c(region)) %>% 
  
  # merge in site level predictors
  left_join(phy_site,by = join_by(wateryear,region,site,cum))

## Region-site bridge
reg_bridge <- site_df %>% 
  distinct(reg_id,site_id) %>% 
  pull(reg_id)


## Year-level data prep  ##

# Add year id to year level predictors
regyear_df <- site_df %>% 
  distinct(
    region,
    wateryear,
    reg_id,
    year_id
    ) %>% 
  arrange(reg_id,year_id) %>% 
  left_join(
    phy_reg_year,
    by = join_by(region,wateryear)
    )


## Predictor prep  ##

# Select 1st-level predictors (create column of 1 for intercept)
x_df <- site_df %>% 
  mutate(int = 1) %>% 
  select(
    int,
    depth,
    dsldd_int,
    plt_cov_int,
    peri_vol_int
  ) %>% 
  mutate(
    across(!int,~as.numeric(scale(.x)))   # scale and center data
    )

# Select 2nd-level predictors (create column of 1 for intercept)
z_list <- lapply(1:max(reg_bridge), function(r){
  z_df <- regyear_df %>% 
    filter(reg_id == r) %>% 
    mutate(int = 1) %>% 
    select(
      int,
      wet_sum_365day,
      pisc_index
    ) %>% 
    mutate(
      across(!int,~as.numeric(scale(.x)))   # scale and center data
      )
})
z_bind <- abind(z_list,along = 3)
z_data <- aperm(z_bind,c(3,1,2))


## Bundle data for stan  ##

# Stan input list
stan_data <- list(
  M = 60,
  N = nrow(site_df),
  `T` = n_distinct(site_df$year_id),
  S = n_distinct(site_df$site_id),
  R = n_distinct(site_df$reg_id),
  K = ncol(x_df),
  L = ncol(z_data[1,,]),
  y= site_df$y,
  yr = site_df$year_id,
  st = site_df$site_id,
  rg = reg_bridge,
  x = x_df,
  z = z_data
)

# Predictor bridges (scaled to raw used for plotting)
xbridge <-site_df 
zbridge <- regyear_df

message("Stan data prepared")

# Run models  ------------------------------------------------------------------

message("Start stan model run")
stan_out <- stan(
    file = file.path(
      mod_dir,
      "tweedie_mvn_second_level_regyear_effects_opt.stan"
      ),
    data = stan_data,
    iter = 3000,
    warmup = 1000,
    chains =4,
    cores = 4
  )

message("stan model run complete")

# Export model  ----------------------------------------------------------------

# Package model outputs and bridges
out <- list(
  stan_out = stan_out,
  x_bridge = x_bridge,
  z_bridge = z_bridge
)

# Name file and export
out_name <- paste(sp,rs,"out_list.rds",sep = "_")
saveRDS(out,file.path(out_dir,out_name))

