# ------------------------------------------------------------------------------
# HPC package installation
# ------------------------------------------------------------------------------

# For use with SBATCH on Slurm Scheduler. 
# DO NO RUN OUTSIDE OF SLURM.

message("Starting cluster-side dependency check...")


# Establish package directories  -----------------------------------------------

# Number of cores specified in slurm scheduler
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))

# Use the user's configured R library
user_lib <- Sys.getenv("R_LIBS_USER")

if (user_lib == "") {
  stop("R_LIBS_USER is not set.")
}

# Create user library if necessary
dir.create(
  user_lib,
  recursive = TRUE,
  showWarnings = FALSE
)

# Put user library first
.libPaths(c(user_lib, .libPaths()))

message("R version: ", R.version.string)
message("User library: ", user_lib)
message("Library paths:")
print(.libPaths())


# Check if packages are installed  ---------------------------------------------

req_pck <- c("cmdstanr")
installed <- installed.packages()[, "Package"]
missing_pkgs <- req_pck[!(req_pck %in% installed)]


# Install missing packages  ----------------------------------------------------

if (length(missing_pkgs) > 0) {
  message("Installing missing packages: ", paste(missing_pkgs, collapse = ", "))
  if("cmdstanr" %in% missing_pkgs) {
    install.packages(
      "cmdstanr", 
      lib = user_lib,
      repos = c('https://stan-dev.r-universe.dev', getOption("repos")),
      Ncpus = n_cores
      )
  }
} else {
  message("All dependencies are already installed on the cluster!")
}

