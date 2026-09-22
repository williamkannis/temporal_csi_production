#-------------------------------------------------------------------------------
#   Run stan models using HPC cluster
#-------------------------------------------------------------------------------

message("============================================================")
message("Starting hpc_cis_analysis.R")
message("M = ",commandArgs(trailingOnly = TRUE))
message("Start time: ", Sys.time())
message("Slurm job ID: ", Sys.getenv("SLURM_JOB_ID"))
message("Slurm array job ID: ", Sys.getenv("SLURM_ARRAY_JOB_ID"))
message("Slurm array task ID: ", Sys.getenv("SLURM_ARRAY_TASK_ID"))
message("============================================================")

# For use with SBATCH on Slurm Scheduler. 
# DO NO RUN OUTSIDE OF SLURM.

# Housekeeping  ----------------------------------------------------------------
rm(list = ls())

# Load in packages
message("[", Sys.time(), "] Loading packages...")
library(cmdstanr)
library(posterior)
message("[", Sys.time(), "] Packages loaded")

# directories
mod_dir <- "hpc/stan_scripts"
input_dir <- "hpc/data"
out_dir <- "hpc/stan_outputs"


# Select specified dataset  ------------------------------------------------------

# Shell job id for selecting species and response of choice
task <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# Load in command arguments
message("[", Sys.time(), "] Load in command arguments... ")
arg <- commandArgs(trailingOnly = TRUE)
if (length(arg) < 1) {
  stop("Missing required arguments")
}
message("[", Sys.time(), "] Loaded in following arguments: ",arg)

# Select response type
message("[", Sys.time(), "] Selecting model response for arg ...")
mod_arg <- arg[1]
if(!mod_arg %in% c("sample_den", "ptob","biomass","production")) {
  stop("mod argument must be 'sample_den' 'ptob', 'biomass', or 'production'")
}
message("[", Sys.time(), "] Model arguement: ",mod_arg," extracted")


message("[", Sys.time(), "] Selecting input file for array task ", task,
        " and mod arg: ",mod_arg)
file_name <- list.files(
  input_dir,
  # pattern = "\\.rds$",
  pattern = paste0(mod_arg,".*\\.rds$"),
  full.names = TRUE
)[task]
message("[", Sys.time(), "] Input file: ", file_name)


# Load in data
message("[", Sys.time(), "] Loading data...")
data <- readRDS(file_name)
message("[", Sys.time(), "] Data loaded")

# Select model type

message(
  "[", Sys.time(), "] Selecting model type based on response arg: ", mod_arg
  )
if(mod_arg %in% c("sample_den","biomass","production")) {
  mod_name <- "tweedie_mvn_regyear_effects.stan" 
  }

if(mod_arg == "ptob"){
  mod_name <- "hurdle_mvn_regyear_effects.stan"
}
message("[", Sys.time(), "] ",mod_name," selected")

if(grepl("tweedie",mod_name)) {
  # Set selected value for m
  message("[", Sys.time(), "] Setting M value based on shell argument...")
  m_arg <- arg[2]
  if (is.na(m_arg)) {
    stop("Missing required argument M")
  }
  # if (length(m_arg) > 1) {
  #   stop("Multiple values supplied for argument M")
  # }
  M <- as.numeric(m_arg)
  if (is.na(M)) {
    stop("Argument 'value' must be numeric")
  }
  data$M <- M
  message("[", Sys.time(), "] Set M in data to ", M)
}


# Run models  ------------------------------------------------------------------

# Load in model
# mod_name <- "tweedie_mvn_regyear_effects.stan"
message("[", Sys.time(), "] Compiling/loading Stan model...")
mod <- cmdstan_model(file.path(mod_dir,mod_name))
message("[", Sys.time(), "] Stan model ready")

# Run model
n_iter <- 2000
n_warm <- 1000
chains <- 4
message("------------------------------------------------------------")
message(paste0("[", Sys.time(), "] Starting Stan sampling for ", mod_name))
message(paste0("Chains: ", chains))
message(paste0("Warmup iterations: ", n_warm))
message(paste0("Sampling iterations: ", n_iter))
message("------------------------------------------------------------")

start_time <- Sys.time()

out <- mod$sample(
  data = data,
  iter_sampling = n_iter, 
  iter_warmup = n_warm,
  chains = chains,
  parallel_chains = 4,
  init = 0.5
)

end_time <- Sys.time()

message("------------------------------------------------------------")
message("[", end_time, "] Stan sampling complete")
message("Sampling time: ", round(
  difftime(end_time, start_time, units = "mins"), 2
), " minutes")
message("------------------------------------------------------------")


# Diagonstics  -----------------------------------------------------------------

# Check for convergence or sampling issues
message("[", Sys.time(), "] Checking convergence and sampling diagnostics...")

# convergence
conv <- out$summary(NULL, c("rhat", "ess_bulk"))
high_rhat <- sum(conv$rhat > 1.1, na.rm = TRUE)
low_ess   <- sum(conv$ess_bulk < 400, na.rm = TRUE)

# sampling
diag <- out$diagnostic_summary()
div <- sum(diag$num_divergent)
tree <- sum(diag$num_max_treedepth)

# combined issues
issues <- div + high_rhat + low_ess

if (issues == 0) {
  message("============================================================")
  message("No issues dectected:")
  message("  Rhat > 1.1: ", high_rhat)
  message( "Bulk ESS < 400: ", low_ess)
  message("  Divergent transitions: ", div)
  message(" Max tree depth exceeded: ", tree)
  message("============================================================")
} else {
  message("============================================================")
  message("WARNING")
  message("Issues detected:")
  message("  Rhat > 1.1: ", high_rhat)
  message("  Divergent transitions: ", div)
  message("  Max tree depth exceeded: ", tree)
  message("============================================================")
}

# Export  ----------------------------------------------------------------------

# Name file and export

if(grepl("tweedie",mod_name)) {
  out_name <-  file_name |>
    gsub("hpc/data",out_dir,x=_) |>
    gsub(
      "_input_data",
      paste0("_stan_out_M",M),
      x=_)
} else {
  out_name <-  file_name |>
    gsub("hpc/data",out_dir,x=_) |>
    gsub(
      "_input_data",
      "_stan_out",
      x=_)
}

out$save_object(out_name)
message("[", Sys.time(), "] Model successfully saved:")
message("  ", out_name)


# End script  ------------------------------------------------------------------

message("============================================================")
message("Finished hpc_csi_analysis.R run")
message("End time: ", Sys.time())
message("============================================================")

