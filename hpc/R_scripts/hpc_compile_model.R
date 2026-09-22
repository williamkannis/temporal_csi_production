#-------------------------------------------------------------------------------
#   Compile Stan models using HPC cluster
#-------------------------------------------------------------------------------

message("============================================================")
message("Starting hpc_compile_model.R")
message("Start time: ", Sys.time())
message("Slurm job ID: ", Sys.getenv("SLURM_JOB_ID"))
message("============================================================")

# For use with SBATCH on Slurm Scheduler. 
# DO NO RUN OUTSIDE OF SLURM.

# Housekeeping  ----------------------------------------------------------------

# Load packages
message("[", Sys.time(), "] Loading packages...")
library(cmdstanr)
message("[", Sys.time(), "] Packages loaded")

# Model directory
mod_dir <- "hpc/stan_scripts"


# Compile code  ----------------------------------------------------------------

mod_name <- "tweedie_mvn_regyear_effects.stan"
message("[", Sys.time(), "] Compiling/loading Stan model: ",mod_name,"...")
mod <- cmdstan_model(file.path(mod_dir,mod_name))
message("[", Sys.time(), "] Stan model ready")

mod_name <- "hurdle_mvn_regyear_effects.stan"
message("[", Sys.time(), "] Compiling/loading Stan model: ",mod_name,"...")
mod <- cmdstan_model(file.path(mod_dir,mod_name))
message("[", Sys.time(), "] Stan model ready")

# End script  ------------------------------------------------------------------

message("============================================================")
message("Finished hpc_csi_analysis.R run")
message("End time: ", Sys.time())
message("============================================================")

