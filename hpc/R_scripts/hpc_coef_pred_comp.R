message("============================================================")
message(".R")
message("Start time: ", Sys.time())
message("Slurm job ID: ", Sys.getenv("SLURM_JOB_ID"))
message("Slurm array job ID: ", Sys.getenv("SLURM_ARRAY_JOB_ID"))
message("Slurm array task ID: ", Sys.getenv("SLURM_ARRAY_TASK_ID"))
message("============================================================")

# For use with SBATCH on Slurm Scheduler. 
# DO NO RUN OUTSIDE OF SLURM.

# House keeping  ===============================================================

# Load in packages
library(dplyr)

# Directories
input_dir <- "hpc/intermediate"
export_dir <- "hpc/" 


# Load in data  ================================================================

name_prefix <- c(
  "coef_gamma",
  "coef_theta",
  "predictions_theta",
  "predictions_response"
)

### Extract file names for each file type  ###
file_list <- lapply(name_prefix, function(x) {
  
  message("[", Sys.time(), "] Loading in file names for ", x, " ...")
  
  files <- list.files(input_dir,x)
  
  message("[", Sys.time(), "] ",length(files), " files names for ", x, " found")
  
  files

  })
names(file_list) <- name_prefix


### Check that all file types have the same number of files  ###
if(length(unique(lengths(a))) !=1) {
  warning("Some file types have unequal number of species/response files")
}


### Load in data  ###
data_list <- lapply(file_list, function(x) {
  
  lapply(x, function(y) {
    
    message("[", Sys.time(), "] Loading in data for ", y, " ...")
    
    df <- readRDS(file.path(input_dir,y))
    
    message("[", Sys.time(), "] ", y, " loaded")
    
    df
    }) %>% 
    bind_rows()

})



# Exporting  ===================================================================

### Export each file type  ###
message("[", Sys.time(), "] Starting exporting...")

lapply(name_prefix, function(x){
  
  df <- data_list[[x]]
  file_name <- file.path(export_dir,(paste0(x,".rds")))
  saveRDS(df,file_name)
  
  message("[", Sys.time(), "] ", x, " exported")
  
})

message("[", Sys.time(), "] Exports finished")


message("============================================================")
message("Finished .R run")
message("End time: ", Sys.time())
message("============================================================")


