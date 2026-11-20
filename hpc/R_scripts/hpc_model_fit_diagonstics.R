message("============================================================")
message("Start hpc_model_fit_diagonstics.R")
message("Start time: ", Sys.time())
message("Slurm job ID: ", Sys.getenv("SLURM_JOB_ID"))
message("Slurm array job ID: ", Sys.getenv("SLURM_ARRAY_JOB_ID"))
message("Slurm array task ID: ", Sys.getenv("SLURM_ARRAY_TASK_ID"))
message("============================================================")

# For use with SBATCH on Slurm Scheduler. 
# DO NO RUN OUTSIDE OF SLURM.


# House keeping  ===============================================================
rm(list = ls())

### Load Packages and custum functions  ###
message("[", Sys.time(), "] Loading packages...")

library(cmdstanr)

source("hpc/R_scripts/hpc_diag_functions.R")

message("[", Sys.time(), "] Packages loaded")

### Directories  ###
out_dir <- "hpc/stan_outputs"
out_dir <- "stan_outputs"  # TEMP DELETE
data_dir <- "hpc/data"
diag_dir <- "hpc/diagnostics"


# Load in data  ================================================================

###  Extract task id to select model file  ### 
message("[", Sys.time(), "] Extracting task id...")

arg <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

if (length(arg) < 1) {
  stop("Missing required arguments")
}

message("[", Sys.time(), "] Task id ", ar, " loaded")

###  Select model file  ### 
message("[", Sys.time(), "] Selecting model file...")

out_file <- list.files(out_dir)[arg]

message("[", Sys.time(), "] ", out_file, " selected")

###  Load in model output  ### 
message("[", Sys.time(), "] Loading in model output...")

out <- readRDS(file.path(out_dir,out_file))

message("[", Sys.time(), "] Model output loaded")

### Find corresponding data list for model output  ###
message("[", Sys.time(), "] Searching for model input data file...")

mods <- gsub("_stan_out_M100.rds|_stan_out.rds","",out_file)
all_data_files <- list.files(data_dir)
data_file <- all_data_files[grepl(mods,all_data_files)]

if (length(data_file) > 1) {
  stop("Multiple data files contain model name: ", mods )
}

message("[", Sys.time(), "] ", data_file, " selected")


### Load in input data  ###
message("[", Sys.time(), "] Loading input data")

data_list <- readRDS(file.path(data_dir,data_file))

message("[", Sys.time(), "] ", data_file, " loaded")


### Extract model type  ###
message("[", Sys.time(), "] Extracting model type...")

if (grepl("_M",out_file)) {
  mod.type <- "tweedie"
}else{
  mod.type <- "hurdle"
}

message("[", Sys.time(), "] Extracted model type: ", mod.type)


# Check for convergence and sampling issues  ===================================

### convergence  ###
message("[", Sys.time(), "] Checking for convergence issues..")

conv_df <- out$summary()

if(any(!is.na(conv_df$rhat) & conv_df$rhat > 1.01)) {
  stop("Some rhats greater than 1.01")
}

if(any(!is.na(conv_df$ess_bulk) & conv_df$ess_bulk < 400)) {
  stop("Some ess_bulk less than 400")
}

if(any(!is.na(conv_df$ess_tail) & conv_df$ess_tail < 400)) {
  stop("Some ess_tail less than 400")
}

message("[", Sys.time(), "] No convergence issues found")


### Sampling issues  ###
message("[", Sys.time(), "] Checking for sampling issures...")

diag <- out$diagnostic_summary()

if(sum(diag$num_divergent) > 0) {
  
  # stop(sum(diag$num_divergent)," samples ended in divergence")
  warning(sum(diag$num_divergent)," samples ended in divergence")
  write.csv(NULL,file.path(diag_dir,paste0("DIVERGENT_",mods)))
  
} else {
  message("[", Sys.time(), "] No sampling issues found")
}


# Tweedie model fit diagnostics  =======================================================

if (mod.type == "tweedie") {
  
  message("[", Sys.time(), "] Starting tweedie diagnostics...")
  
  ### Extract M value  ###
  message("[", Sys.time(), "] Extracting M value ...")
  
  M <- as.numeric(sub(".*_M([0-9.]+)\\..*", "\\1", out_file))
  
  if (is.na(M)) {
    stop("make sure M is included in file name")
  }
  
  message("[", Sys.time(), "]  M of ", M, " extracted")
  
  ### M value diagnostics  ###
  message("[", Sys.time(), "] Starting M value diagnostics......")
  
  lambda_draws <- posterior::as_draws_matrix(out$draws("lambda"))
  max_lambda <- quantile(lambda_draws, 0.999)
  M_required <- qpois(
    1 - 1e-10,
    lambda = max_lambda
  )
  
  print(data.frame(
    tolerance = c(1e-6, 1e-8, 1e-10, 1e-12),
    M = c(
      qpois(1 - 1e-6, lambda_max),
      qpois(1 - 1e-8, lambda_max),
      qpois(1 - 1e-10, lambda_max),
      qpois(1 - 1e-12, lambda_max)
    )
  )
  )
  
  if (M_required > M) {
    warning(
      "Provided M ", 
      M, 
      " is lower than recommened M of ", 
      M_required, 
      " at the level of 1e-10"
      )
    write.csv(NULL,file.path(diag_dir,paste0("LOW_M_",mods)))
  } else {
    
  }
  
  message("[", Sys.time(), "] No issues with M")
  
  ### Run model fit function  ###
  message("[", Sys.time(), "] Running model fit function...")
  
  diag <- tweedie_diagnostics(
    fit = out,
    y = data_list$all_production_mean$y,
    st = data_list$all_production_mean$st,
    yr = data_list$all_production_mean$yr,
    rg = data_list$all_production_mean$rg
  )
  
  message("[", Sys.time(), "] Model fit function ran successfully")
  
}




# Hurdle model fit diagnostics  ================================================
if (mod.type == "hurdle") {
  message("[", Sys.time(), "] Starting Hurdle diagnostics...")
  
  ### Run model fit function  ###
  message("[", Sys.time(), "] Running model fit function...")
  
  diag <- hurdle_diagnostics(
    fit = out,
    y = data_list$all_production_mean$y,
    st = data_list$all_production_mean$st,
    yr = data_list$all_production_mean$yr,
    rg = data_list$all_production_mean$rg
  )
  
  message("[", Sys.time(), "] Model fit function ran successfully")
}


# Export model fit diagnostics  ================================================

message("[", Sys.time(), "] Exporting model fit diagnostics...")

diag_file <- paste0("diag_",mods,".rds")
saveRDS(diag, file.path(diag_dir, diag_file))

message("[", Sys.time(), "] ", diag_file, " exported")


# End of script  ===============================================================

message("============================================================")
message("Finished hpc_model_fit_diagonstics.R run")
message("End time: ", Sys.time())
message("============================================================")

