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

# Select specied dataset  ------------------------------------------------------

# Shell argument for selecting species and response of choice
arg <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

data <- list.files()[arg]



# Run models  ------------------------------------------------------------------

message("Start stan model run")
stan_out <- mod$
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

