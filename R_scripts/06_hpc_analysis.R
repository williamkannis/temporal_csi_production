#-------------------------------------------------------------------------------
#
#   Prepare data for HPC analysis scripts
#
#-------------------------------------------------------------------------------

# AUTHOR: William K. Annis

# CREATED: Sept 3, 2026

# DESCRIPTION: 


# Housekeeping  ----------------------------------------------------------------
rm(list = ls())

# Load in packages
library(ssh)


# HPC configuration  -----------------------------------------------------------

# Connect to cluster
session <- ssh_connect("wka25@hpc-login.rcc.fsu.edu")


# Hpc uploads  -----------------------------------------------------------------

scp_upload(
  session = session,
  files = "hpc",
  to = "/gpfs/home/wka25/"
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
        " Rscript /gpfs/home/wka25/hpc/R_scripts/hpc_install_packages.R"
      )
    )
  )
)


# Run analyses on cluster  -----------------------------------------------------

# Set M value for all jobs
M <- 60

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
        "/hpc/R_scripts/hpc_compile_model.R"
      )
    )
  )
)

# Model runs
ssh_exec_wait(
  session,
  command = paste(
    "sbatch",
    "--array=1-14",
    "--cpus-per-task=4",
    "--mem=8G",
    "--time=3:00:00",
    "--job-name=production",
    "--output=/gpfs/home/wka25/my_project/results/slurm-%A_%a.out",
    "--wrap",
    shQuote(
      paste0(
        "module load gnu/13 && module load R/4.4.0 && Rscript /gpfs/home/",
        "wka25/hpc/R_scripts/hpc_csi_analysis.R ",
        M
      )
    )
  )
)


# Clean up cluster  ------------------------------------------------------------

# Download model outputs
scp_download(
  session = session,
  files = "/gpfs/home/wka25/hpc/stan_outputs",
  to = "."
)

# Remove all scripts and data. 
# WARNING: ENSURE THAT RESULTS HAVE DOWNLOADED FIRST
ssh_exec_wait(
  session,
  command = "rm -rf /gpfs/home/wka25/hpc"
)

