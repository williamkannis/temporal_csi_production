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
library(ssh)

# directories
mod_dir <- "stan_scripts"
input_dir <- "prod_data"
export_dir <- "hpc/stan_input"

# Data
prod_df <- 
  readRDS(file.path(input_dir,"fsprod_igr_2026-07-13.rds"))
phy_site <- 
  readRDS(file.path(input_dir,"phys_site_predictors_2026-08-17.rds"))


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
nrow(
  expand.grid(
    unique(samp_df$species),
    c("production_mean","biomass_mean","ptob")
    )
  )
# Use this value in batch array shell script

# HPC configuration  -----------------------------------------------------------

library(ssh)
# Connect to cluster
session <- ssh_connect("wka25@hpc-login.rcc.fsu.edu")


# Hpc uploads  -----------------------------------------------------------------
# Upload package scripts (ADD ALL SCRIPTS)
scp_upload(
  session = session,
  files = "hpc_scripts/hpc_install_packages.R",
  to = "/gpfs/home/wka25/hpc_install_packages.R"
)

# Upload data
scp_upload(
  session = session,
  files = "hpc_scripts/hpc_install_packages.R",
  to = "/gpfs/home/wka25/hpc_install_packages.R"
)

# Install packages  ------------------------------------------------------------
ssh_exec_wait(
  session,
  command = paste(
    "srun",
    "--cpus-per-task=4",
    "--mem=8G",
    "--time=00:20:00",
    "bash -lc",
    shQuote(
      paste0(
        "module load gnu/13 && module load R/4.4.0 && module load webproxy &&",
        " Rscript /gpfs/home/wka25/hpc_install_packages.R"
      )
    )
  )
)




# Run analyses on cluster  -----------------------------------------------------

# Compile models
ssh_exec_wait(
  session,
  command = paste(
    "srun",
    "--cpus-per-task=1",
    "--mem=8G",
    "--time=00:02:00",
    "bash -lc",
    shQuote(
      paste0(
        "module load gnu/13 && module load R/4.4.0 && Rscript /gpfs/home/wka25",
        "/hpc_compile_model.R"
      )
    )
  )
)

# Model runs
ssh_exec_wait(
  session,
  command = paste(
    "sbatch",
    "--array=1-21",
    "--cpus-per-task=4",
    "--mem=8G",
    "--time=2:00:00",
    "--job-name=production",
    "--output=/gpfs/home/wka25/my_project/results/slurm-%A_%a.out",
    "--wrap",
    shQuote(
      paste0(
        "module load gnu/13 && module load R/4.4.0 && Rscript /gpfs/home/",
        "wka25/my_project/R/hpc_csi_analysis.R"
      )
    )
  )
)

# Retrieve results  ------------------------------------------------------------

# Download model outputs


# export model outputs
saveRDS()


# Clean up cluster  ------------------------------------------------------------

# Delete package installation script
ssh_exec_wait(
  session,
  command = "rm /gpfs/home/wka25/hpc_install_packages.R"
)

