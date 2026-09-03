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

# directories
mod_dir <- "hpc_scripts"
export_dir <- "hpc/stan_output"


# HPC configuration  -----------------------------------------------------------

# Connect to cluster
session <- ssh_connect("wka25@hpc-login.rcc.fsu.edu")


# Hpc uploads  -----------------------------------------------------------------

# # Upload package scripts (ADD ALL SCRIPTS)
# scp_upload(
#   session = session,
#   files = "hpc_scripts/hpc_install_packages.R",
#   to = "/gpfs/home/wka25/hpc_install_packages.R"
# )
# 
# # Upload data
# scp_upload(
#   session = session,
#   files = "hpc_scripts/hpc_install_packages.R",
#   to = "/gpfs/home/wka25/hpc_install_packages.R"
# )

# Upload data and scripts
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


# Clean up cluster  ------------------------------------------------------------

# Download model outputs
ssh_exec_wait()

# Remove all scripts and data. 
# WARNING: ENSURE THAT RESULTS HAVE DOWNLOADED FIRST
ssh_exec_wait(
  session,
  command = "rm /gpfs/home/wka25/hpc"
)

