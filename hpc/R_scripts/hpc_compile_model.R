#-------------------------------------------------------------------------------
#   Compile Stan models using HPC cluster
#-------------------------------------------------------------------------------

# For use with SBATCH on Slurm Scheduler. 
# DO NO RUN OUTSIDE OF SLURM.

# Housekeeping  ----------------------------------------------------------------

# Load packages
library(cmdstanr)

# Model directory
mod_dir <- "hpc/stan_scripts"

# Compile code  ----------------------------------------------------------------

message("Start code compiling")
mod <- cmdstan_model()