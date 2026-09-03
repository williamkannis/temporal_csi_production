#-------------------------------------------------------------------------------
#
#   Prepare data for HPC analysis scripts
#
#-------------------------------------------------------------------------------

# AUTHOR: William K. Annis

# CREATED: July 13, 2026

# DESCRIPTION: 


# Housekeeping  ----------------------------------------------------------------
rm(list = ls()) 

# Load in packages
library(dplyr)
library(purrr)
library(abind)
library(jsonlite)

# directories
mod_dir <- "stan_scripts"
input_dir <- "prod_data"
export_dir <- "hpc/stan_input"

# Data
prod_df <- 
  readRDS(file.path(input_dir,"fsprod_igr_2026-07-13.rds"))
phy_site <- 
  readRDS(file.path(input_dir,"phys_site_predictors_2026-08-17.rds"))
phy_reg_year <- 
  readRDS(file.path(input_dir,"phys_regionyear_predictors_2026-08-17.rds"))


# Prep data  -------------------------------------------------------------------

# Remove NA production estimates
prod_for <- prod_df %>% filter(!is.na(production_mean))

# Create a composite production measure using all species
prod_all <- prod_for %>% 
  group_by(site,cum,date,area,interval) %>% 
  summarise(
    across(contains("sample_den"),sum),
    across(contains("biomass"),sum),
    across(contains("production"),sum)
  ) %>% 
  mutate(species = "all") %>% 
  bind_rows(prod_for) %>% 
  mutate(
    ptob = production_mean/interval_biomass_mean,
    ptob = case_when(
      is.nan(ptob) ~ 0,
      T ~ ptob
    )
  )

# Add sample info to production data
samp_df <- phy_site %>% 
  distinct(wateryear,region,site,cum) %>% 
  filter(wateryear != 2024) %>%   ## TEMPORARY ASK NATE FOR NEWEST SHARK RIVER DATA (2025)
  filter(wateryear !=1995) %>%  ## TEMP ASK JOEL FOR LAG DATA FOR THIS YEAR
  left_join(prod_all, by = join_by(site,cum)) %>% 
  filter(!is.na(production_mean)) %>% 
  
  # transform response varibales to improv convergence
  mutate(
    production_mean = production_mean*1000,
    biomass_mean = biomass_mean*1000,
    ptob = ptob*1000 
  )

# Number of species response combinations
cb <- expand.grid(
  unique(samp_df$species),
  c("production_mean","biomass_mean","ptob")
)
n_cb <- nrow(cb)


# Model input lists  -----------------------------------------------------------

# Create a stan input data list for each species
input_list <- lapply(1:n_cb, function(i){
  
  # Extract species and response name
  s <- cb[i,1]
  r <- cb[i,2]
  
  # Create column y for response of choice
  samp_df$y <- samp_df[[r]]

  # Extract data for selected specices
  site_df <- samp_df %>% 
    filter(species == s) %>% 
    
    # Create site and year id for random effects
    arrange(region,site,cum) %>% 
    mutate(site_id = cur_group_id(),.by = c(region,site)) %>% 
    mutate(year_id = cur_group_id(),.by = c(wateryear)) %>% 
    # mutate(regyear_id = cur_group_id(),.by = c(region,wateryear)) %>% 
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
    # distinct(region,wateryear,regyear_id) %>% 
    distinct(region,wateryear,reg_id,year_id) %>% 
    # arrange(regyear_id) %>% 
    arrange(reg_id,year_id) %>% 
    left_join(phy_reg_year,by = join_by(region,wateryear))
  
  
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
    
    # scale and center data
    mutate(across(!int,~as.numeric(scale(.x))))
  
  # Select 2nd level predictors (create column of 1 for intercept) and create
  # array for region specific
  # z_df <- regyear_df %>% 
  #   mutate(int = 1) %>% 
  #   select(
  #     int,
  #     wet_sum_365day
  #   ) %>% 
  # 
  # # scale and center data
  # mutate(across(!int,~as.numeric(scale(.x))))
  
  z_list <- lapply(1:max(reg_bridge), function(r){
    z_df <- regyear_df %>% 
      filter(reg_id == r) %>% 
      mutate(int = 1) %>% 
      select(
        int,
        wet_sum_365day,
        # wet_sum_365day_delta,
        pisc_index
        # pisc_index_lag
        # PC1,
        # PC2
      ) %>% 
      
      # scale and center data
      mutate(across(!int,~as.numeric(scale(.x))))
  })
  z_bind <- abind(z_list,along = 3)
  z_data <- aperm(z_bind,c(3,1,2))
  
  
  ## Stan list  ##
  # stan <-list(
  #   N = nrow(site_df),
  #   L = n_distinct(site_df$year_id),
  #   O = n_distinct(site_df$site_id),
  #   P = n_distinct(site_df$regyear_id),
  #   K = ncol(x_df),
  #   J = ncol(z_df),
  #   y= log(1000*site_df$production_mean+1),
  #   # y=site_df$production_mean*1000,
  #   pp = site_df$regyear_id,
  #   ll = site_df$year_id,
  #   oo = site_df$site_id,
  #   x = x_df,
  #   z = z_df
  # )
  
  stan_data <- list(
    # M = 50,
    M = 60,
    N = nrow(site_df),
    `T` = n_distinct(site_df$year_id),
    S = n_distinct(site_df$site_id),
    R = n_distinct(site_df$reg_id),
    K = ncol(x_df),
    L = ncol(z_data[1,,]),
    y = site_df$y,
    yr = site_df$year_id,
    st = site_df$site_id,
    rg = reg_bridge,
    x = x_df,
    z = z_data
  )
  
  list(stan_data = stan_data, xbridge = site_df, zbridge = regyear_df)
  
}
)


# Format input lists  ----------------------------------------------------------

names(input_list) <- apply(cb, 1, paste, collapse = "_")
input_list_t <- transpose(input_list)

# For anlaysis on HPC
stan_list <- input_list_t$stan_data

# For plot labels
xbridge_list <-input_list_t$xbridge
zbridge_list <-input_list_t$zbridge


# Export  ----------------------------------------------------------------------

# For analysis on HPC
lapply(1:n_cb,function(i){
  file_name <- paste0(
    names(stan_list)[i],
    "_input_data.json"
    )
  write_json(
    stan_list[[i]],
    file.path(export_dir,file_name)
  )
}
)

