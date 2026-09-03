#-------------------------------------------------------------------------------
#   Run stan models using HPC cluster
#-------------------------------------------------------------------------------

message("Start model run script")

# For use with SBATCH on Slurm Scheduler. 
# DO NO RUN OUTSIDE OF SLURM.

# Housekeeping  ----------------------------------------------------------------

# Load in packages
library(cmdstanr)

# directories
mod_dir <- "hpc/stan_scripts"
input_dir <- "hpc/data"
out_dir <- "hpc/stan_outputs"

# Select specified dataset  ------------------------------------------------------

# Shell argument for selecting species and response of choice
arg <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

file_name <- list.files(
  input_dir,
  pattern = "\\.json$",
  full.names = TRUE
  )[arg]


# Run models  ------------------------------------------------------------------

# Load in model
message("Compiling code")
mod <- cmdstan_model()

# Run model
message("Start stan model run")
out <- mod$sample(
  data = file_name,
  iter_sampling = 2000, 
  iter_warmup = 1000,
  chains = 4,
  parallel_chains = 4 
)
message("stan model run complete")


# Export model  ----------------------------------------------------------------

# Name file and export
out_name <-  file_name |>
  gsub("hpc/data",out_dir,x=_) |>
  gsub("_input_data.json","_stan_out.rds",x=_)
  
saveRDS(out,out_name)

