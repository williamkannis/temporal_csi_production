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


# Prepare for analysis  --------------------------------------------------------

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




# Run analyses on cluster  -----------------------------------------------------

session <- ssh_connect("wka25@hpc-login.rcc.fsu.edu")

# Set M value for all jobs
M <- 100

# PtoB Model runs
ssh_exec_wait(
  session,
  command = paste(
    "sbatch",
    "--array=1-7",
    "--cpus-per-task=4",
    "--mem=16G",
    "--time=0:30:00",
    "--job-name=ptob",
    "--output=/gpfs/home/wka25/my_project/results/slurm-%x_%A_%a.out",
    "--mail-type=ALL",
    "--wrap",
    shQuote(
      paste0(
        "module load gnu/13 && module load R/4.4.0 && Rscript /gpfs/home/",
        "wka25/hpc/R_scripts/hpc_csi_analysis.R ptob"
      )
    )
  )
)


# Size Model runs
ssh_exec_wait(
  session,
  command = paste(
    "sbatch",
    "--array=1-7",
    "--cpus-per-task=4",
    "--mem=16G",
    "--time=0:30:00",
    "--job-name=interval_mean_wt",
    "--output=/gpfs/home/wka25/my_project/results/slurm-%x_%A_%a.out",
    "--mail-type=ALL",
    "--wrap",
    shQuote(
      paste0(
        "module load gnu/13 && module load R/4.4.0 && Rscript /gpfs/home/",
        "wka25/hpc/R_scripts/hpc_csi_analysis.R interval_mean_wt"
      )
    )
  )
)


# Density model runs
ssh_exec_wait(
  session,
  command = paste(
    "sbatch",
    "--array=1-7",
    "--cpus-per-task=4",
    "--mem=16G",
    "--time=24:00:00",
    "--job-name=interval_density",
    "--output=/gpfs/home/wka25/my_project/results/slurm-%x_%A_%a.out",
    "--mail-type=ALL",
    "--wrap",
    shQuote(
      paste0(
        "module load gnu/13 && module load R/4.4.0 && Rscript /gpfs/home/",
        "wka25/hpc/R_scripts/hpc_csi_analysis.R interval_density ",
        M
      )
    )
  )
)



# Biomass Model runs
ssh_exec_wait(
  session,
  command = paste(
    "sbatch",
    "--array=1-7",
    "--cpus-per-task=4",
    "--mem=16G",
    "--time=24:00:00",
    "--job-name=interval_biomass_mean",
    "--output=/gpfs/home/wka25/my_project/results/slurm-%x_%A_%a.out",
    "--mail-type=ALL",
    "--wrap",
    shQuote(
      paste0(
        "module load gnu/13 && module load R/4.4.0 && Rscript /gpfs/home/",
        "wka25/hpc/R_scripts/hpc_csi_analysis.R interval_biomass_mean ",
        M
      )
    )
  )
)

# Production Model runs
ssh_exec_wait(
  session,
  command = paste(
    "sbatch",
    "--array=2",
    "--cpus-per-task=4",
    "--mem=16G",
    "--time=24:00:00",
    "--job-name=production_mean",
    "--output=/gpfs/home/wka25/my_project/results/slurm-%x_%A_%a.out",
    "--mail-type=ALL",
    "--wrap",
    shQuote(
      paste0(
        "module load gnu/13 && module load R/4.4.0 && Rscript /gpfs/home/",
        "wka25/hpc/R_scripts/hpc_csi_analysis.R production_mean ",
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

