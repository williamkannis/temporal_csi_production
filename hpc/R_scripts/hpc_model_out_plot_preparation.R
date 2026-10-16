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

# Packages
message("[", Sys.time(), "] Loading packages...")

library(dplyr)
library(stringr)
# library(tidytable) use :: notation instead to avoid conflicts
library(tidyr)

message("[", Sys.time(), "] Packages loaded")


# Directories
out_dir <- "stan_outputs"
data_dir <- "hpc/data"
export_dir <- "hpc/intermediate"



# Data loading  ================================================================

message("============================================================")
message("Start time: ", Sys.time())
message("Starting data loading section")
message("============================================================")


### Load in command arguments  ###
message("[", Sys.time(), "] Load in command arguments... ")

arg <- as.numeric(commandArgs(trailingOnly = TRUE))
# arg <- 2  ### REMOVE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#!#!#!#!#!#!#!#!!#!#!#!#!!#
#!!!!!#!!!!!!!!#!!!!!!!!#!!!!!!!!!!!!#!#!!!#!

if (length(arg) < 1) {
  stop("Missing required arguments")
}

message("[", Sys.time(), "] Loaded in following arguments: ",arg)


### Select in model file based on slurm argument  ###
message(
  "[", Sys.time(), "] Selecting model file based on command arguments... "
  )

out_file <- list.files(out_dir)[arg]

message("[", Sys.time(), "] ", out_file, " selected")


### Load in model outputs  ###
message("[", Sys.time(), "] Loading model outputs for ", out_file)

out <-readRDS(file.path(out_dir,out_file))

message("[", Sys.time(), "] ", out_file, " loaded")


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


### Create model indexing names  ###
message("[", Sys.time(), "] Creating model indexing names...")

sp <- stringr::str_split_fixed(mods, "_", n = 2)[1]
response <- stringr::str_split_fixed(mods, "_", n = 2)[2]

message(
  "[", Sys.time(), 
  "] Indices for species: ", sp, " and response: ", response, " created"
  )


### Model prediction settings
message("[", Sys.time(), "] Setting model prediction settings...")

pred_len <- 100
n_iter <- 1000

message(
  "[", Sys.time(),
  "] Model prediction will be created along a range of ",
  pred_len,
  " predictor values using ",
  n_iter,
  " posterior draws"
  )



# Input data preparation  ======================================================

message("============================================================")
message("Start time: ", Sys.time())
message("Starting Input data preparation")
message("============================================================")


# Z data  ---------------------------------------------------------------------

# Standardized and raw second level predictor data

### Change into a single data.frame with columns for region and wateryear  ###
message("[", Sys.time(), "] Preparing z input data...")

z_data <- lapply(1:dim(data_list$z)[[1]], function(x) {
  rg <- dimnames(data_list$z)[[1]][x]
  as.data.frame(data_list$z[x,,]) %>% 
    tibble::rownames_to_column("wateryear") %>% 
    mutate(
      region = rg,
      wateryear = as.numeric(wateryear)
    )
  
}) %>% 
  bind_rows()
z_data_raw <- data_list$z_raw

message("[", Sys.time(), "] Z input data prepared")


# X data  ----------------------------------------------------------------------

# Standardized and raw second level predictor data

### Merge raw data into standardized data  ###
message("[", Sys.time(), "] Preparing x input data...")

x_data <- data_list$x
x_data_raw <- data_list$x_raw

message("[", Sys.time(), "] x input data prepared")



# Effect coefficient preparation  ==============================================

message("============================================================")
message("Start time: ", Sys.time())
message("Effect coefficient preparation")
message("============================================================")


# Gamma coefficient summary  ---------------------------------------------------

# Mean and 95% credible intervals for global first level effects, and cross-
# scale interaction coefficients.

### Extract mean and 95% CIs  ###
message("[", Sys.time(), "] Summarizing gamma coefficients...")

gamma_coef <- out$summary(
  "gamma",
  mean,
  quantile, 
  .args = list(probs = c(0.025, 0.975))
) %>% 
  rename(
    lwr = `2.5%`,
    upr = `97.5%`
  ) %>% 
  mutate(
    overlap0 = case_when(
      lwr*upr >0 ~ F,
      T ~ T
    ),
    
    # Create model indexing columns
    species = sp,
    response = response,
  ) %>% 
  
  # Bridge coef names with x and z variable names
  tidytable::separate_wider_regex(
    cols = variable,
    patterns = c(
      ".*",          
      "\\[",         
      z_idx = "\\d+", 
      ",",           
      x_idx = "\\d+", 
      "\\]"          
    ),
    cols_remove = FALSE
  ) %>% 
  mutate(
    z_idx = as.numeric(z_idx),
    x_idx = as.numeric(x_idx)
  ) %>% 
  
  # Z var bridge
  left_join(
    data.frame(
      z_idx = seq_len(ncol(data_list$z[1,,])),
      z_var = colnames(data_list$z[1,,])
    ), 
    by = "z_idx"
  ) %>% 
  
  # X var bridge
  left_join(
    data.frame(
      x_idx = seq_len(ncol(data_list$x)),
      x_var = colnames(data_list$x)
    ), 
    by = "x_idx"
  ) %>% 
  select(-x_idx,-z_idx)

message("[", Sys.time(), "] Gamma coef summary created")


### Create gamma bridge linking gamma parameters with predictor names  ###
  message("[", Sys.time(), "] Creating gamma-varaible name bridge...")
gamma_bridge <- gamma_coef %>% 
  distinct(variable,z_var,x_var)
message("[", Sys.time(), "] Gamma-varaible name bridge created")


# Theta coefficient summary  ---------------------------------------------------
  
# Region/year specific first level effect coefficients

### Extract mean and 95% CIs  ###
message("[", Sys.time(), "] Summarizing Theta coefficients...")

theta_coef <- out$summary(
  "beta",  # CHANGE TO THETA IN MODEL FOR FINAL RUN
  mean,
  quantile, 
  .args = list(probs = c(0.025, 0.975))
) %>% 
  rename(
    upr = `97.5%`,
    lwr = `2.5%`
  ) %>% 
  mutate(
    overlap0 = case_when(
      lwr*upr >0 ~ F,
      T ~ T
    ),
    
    # Create model indexing columns
    species = sp,
    response = response,
  ) %>% 
  
  # Bridge coef names with x variable names, regions, and wateryear
  tidytable::separate_wider_regex(
    cols = variable,
    patterns = c(
      ".*",          
      "\\[",         
      r_idx = "\\d+", 
      ",",           
      x_idx = "\\d+", 
      ",",           
      y_idx = "\\d+", 
      "\\]"          
    ),
    cols_remove = FALSE
  ) %>% 
  mutate(
    x_idx = as.numeric(x_idx),
    r_idx = as.numeric(r_idx),
    y_idx = as.numeric(y_idx)
  ) %>% 
  
  # X var bridge
  left_join(
    data.frame(
      x_idx = seq_len(ncol(data_list$x)),
      x_var = colnames(data_list$x)
    ), 
    by = "x_idx"
    ) %>% 
  
  # Region name bridge
  left_join(
    r_name <- data.frame(
      r_idx = 1:dim(data_list$z)[[1]],
      region = dimnames(data_list$z)[[1]]
    ), 
    by = "r_idx"
    ) %>% 
  
  # Wateryear name bridge
  left_join(
    data.frame(
      y_idx = 1:dim(data_list$z)[[2]],
      wateryear = dimnames(data_list$z)[[2]]
    ), 
    by = "y_idx"
    ) %>% 
  select(-x_idx,-r_idx,-y_idx)

message("[", Sys.time(), "] Theta coef summary created")

# # Create beta bridge linking gamma parameters with predictor names
# message("[", Sys.time(), "] Creating Theta-varaible name bridge...")
# theta_bridge <- theta_coef %>% 
#   distinct(variable,x_var,region,wateryear)
# message("[", Sys.time(), "] Theta-variable name bridge created")



# Model predictions  ===========================================================

message("============================================================")
message("Start time: ", Sys.time())
message("Effect coefficient preparation")
message("============================================================")


# Predicted thetas (full distribution)  ----------------------------------------

# Use second level gamma parameters to predict 1st level slopes across predictors

### Extract gamma model draws  ###
message("[", Sys.time(), "] Extracting all gamma model draws...")

gamma_draws <- out$draws("gamma",format = "draws_matrix")
n_draws <- nrow(gamma_draws)

message("[", Sys.time(), "] ", n_draws, " gamma model draws extracted")

  
### Find all combinations of variables  ###
message("[", Sys.time(), "] Extracting all combinations of x and z drivers...")

var_combo_theta <- gamma_coef %>% 
  distinct(z_var,x_var) %>% 
  filter(z_var != "int")

message("[", Sys.time(), "] ", nrow(var_combo_theta), " combinations extracted")  


### Create n iterations of theta predictions using gamma coef  ###
message(
  "[", Sys.time(),
  "] Predicting theta values across ", 
  pred_len, 
  " values of each z predictor"
  )
tpred_start <- Sys.time()

theta_pred_raw <- lapply(seq_len(n_iter), function(i) {
  iter <- sample(seq_len(n_draws),1,replace = T)
  
  lapply(1:nrow(var_combo_theta), function(v) {
    vars <- var_combo_theta[v,]
    predictor <- vars$z_var
    
    # Create z predictor inputs
    z_mat <- matrix(
      0,
      nrow = pred_len+1,
      ncol = ncol(z_data %>% select(-region,-wateryear))
    )
    colnames(z_mat) <- colnames(z_data %>% select(-region,-wateryear))
    z <- z_data[,predictor]
    z_mat[,predictor] <- c(
      seq(min(z),max(z),length.out = pred_len),
      0  # global effect
    )
    z_mat[,"int"] <- 1
    
    # Prepare raw data for plotting
    z_raw <- z_data_raw[[predictor]]
    z_raw_range <- c(
      seq(min(z_raw),max(z_raw),length.out = pred_len),
      mean(z_raw)  # global effect
    )
    
    
    # Extract a draw of gamma, and format
    gamma <- as.data.frame(gamma_draws[iter,]) %>% 
      pivot_longer( 
        cols = everything(), 
        names_to = "variable", 
        values_to = "coef"
      ) %>% 
      left_join(
        gamma_coef %>% distinct(variable,x_var,z_var),
        by = join_by(variable)
      ) %>% 
      filter(x_var == vars$x_var) %>% 
      pull(coef)
    
    # Create prediction
    pred <- z_mat %*% gamma
    
    # Save prediction and other info in data.frame
    data.frame(
      iter = i,
      response = response,
      species = sp,
      x_var = vars$x_var,
      z_var = predictor,
      z_range = z_mat[,predictor],
      z_raw_range = z_raw_range,
      pred = pred
      ) %>%
      mutate(
        z_cat = case_when(
          z_range == min(z_range) ~ "min",
          z_range == 0 ~ "mean",  #global effect
          z_range == max(z_range) ~ "max"
        )
    )
  }
  ) %>% 
    bind_rows()
}
)

tpred_end <- Sys.time()
message(
  "[", Sys.time(), 
  "] Theta predictions created in ", 
  round(difftime(tpred_end, tpred_start, units = "mins"), 2),
  " minutes"
)


# Predicted thetas (summarized)  -----------------------------------------------

# Estimate mean and 95% Credible intervals for each theta prediction
message("[", Sys.time(), "] Summarizing theta predictions...")

theta_pred <- theta_pred_raw %>% 
  bind_rows() %>% 
  select(-iter) %>% 
  group_by(across(-pred)) %>% 
  filter(is.na(z_cat) | z_cat != "mean") %>%   # remove the global effect
  summarize(
    pred_md = quantile(pred, probs = 0.5, na.rm = TRUE),
    pred_lo = quantile(pred, probs = 0.025, na.rm = TRUE),
    pred_up = quantile(pred, probs = 0.975, na.rm = TRUE),
    .groups = "drop"
    )

message("[", Sys.time(), "] Theta predictions summarized")

# Predicted responses (raw; min,med,max CSI)  ----------------------------------

# Use predicted thetas to predict responses across 1st level predictors. 
# Predictions are summarized with mean and 95% credible intervals

# X and Z combos used for response predictions
message(
  "[", Sys.time(),
  "] Extracting all combinations of x and z drivers for response predictions..."
  )

var_combo_response <- var_combo_theta %>% 
  filter(x_var !="int") %>% 
  tidyr::crossing(z_cat = c("min","mean","max")) %>% 
  as.data.frame()

message(
  "[", Sys.time(),
  "] Extracted ", 
  nrow(var_combo_response), 
  " combinations"
  )


# Predict responses for each x driver, z interaction, and z 
# category (i.e., min, max, mean) 
message(
  "[", Sys.time(),
  "] Predicting response values across ", 
  pred_len, 
  " values of each x predictor"
)
rpred_start <- Sys.time()

response_pred_raw <- lapply(1:nrow(var_combo_response), function(p) {
  
  # Extract indicated variables
  x_v <- var_combo_response[p,2]
  z_v <- var_combo_response[p,1]
  z_r <- var_combo_response[p,3]
  
  # Prepare input predictor data
  x_mat <- matrix(
    0,
    nrow = pred_len,
    ncol = ncol(x_data)
  )
  colnames(x_mat) <- colnames(x_data)
  x_vec <- x_data[,x_v]
  x_mat[,x_v] <- seq(min(x_vec),max(x_vec),length.out = pred_len)
  x_mat[,"int"] <- 1
  
  # Prepare raw data for plotting
  x_raw <- x_data_raw[[x_v]]
  x_raw_range <- seq(min(x_raw),max(x_raw),length.out = pred_len)
  
  # Create predictions for every iteration of gamma predictions
  pred_list <- lapply(seq_along(theta_pred_raw), function(i){
    
    # Extract theta coef
    df <- theta_pred_raw[[i]]
    theta <- df %>% 
      filter(
        z_var == z_v,
        z_cat == z_r
      ) %>% 
      pull(pred)
    
    # # Create prediction on normal scale
    # exp(x_mat %*% theta)
    
    # Create prediction on log scale
    pred <- x_mat %*% theta
    
    # Prepare output into dataframe
    data.frame(
      iter = unique(df$iter),
      species = sp,
      response = response,
      x_var = x_v,
      z_var = z_v,
      z_cat = z_r,
      x_range=x_mat[,x_v],
      x_raw_range = x_raw_range,
      pred = pred
    )
  }) 
})%>% 
  bind_rows()

rpred_end <- Sys.time()
message(
  "[", Sys.time(), 
  "] Response predictions created in ", 
  round(difftime(rpred_end, rpred_start, units = "mins"), 2),
  " minutes"
  )


# Predicted responses (summary; min,med,max CSI)  ------------------------------

# Estimate mean and 95% Credible intervals for each response prediction
message("[", Sys.time(), "] Summarizing response predictions...")

response_pred <- response_pred_raw %>% 
  select(-iter) %>% 
  group_by(across(-pred)) %>% 
  summarize(
    pred_md = quantile(pred, probs = 0.5, na.rm = TRUE),
    pred_lo = quantile(pred, probs = 0.025, na.rm = TRUE),
    pred_up = quantile(pred, probs = 0.975, na.rm = TRUE),
    .groups = "drop"
  )
 
message("[", Sys.time(), "] Response predictions summarized") 


# Proportion change predictions  -----------------------------------------------

# Estimate the proportion change of each response across the entire range of 
# each driver based on predicted responses 
message("[", Sys.time(), "] Estimating proportion change and difference...")

prop_change <- response_pred_raw %>% 
  select(-x_raw_range) %>% 
  
  # Assign tags for min and max x values
  group_by(across(-c(x_range,pred))) %>% 
  mutate(
    x_cat = case_when(
      x_range == min(x_range) ~ "min",
      x_range == max(x_range) ~ "max",
    )
  ) %>% 
  ungroup() %>% 
  
  # Add in annual drivers, formatted to match seasonal drivers
  bind_rows(
    theta_pred_raw %>% 
      bind_rows() %>% 
      filter(x_var == "int") %>% 
      mutate(
        x_var = z_var,
        x_cat = z_cat,
        z_cat = "mean" # TEMP
      ) 
  ) %>% 
  
  # Estimate the log proportion change and differnce between highest and lowest 
  # value of each driver
  filter(x_cat %in% c("max", "min")) %>% 
  select(iter,species,response,z_var,z_cat,x_var,x_cat,pred) %>% 
  pivot_wider(names_from = x_cat,values_from = pred) %>% 
  mutate(
    prop_change = log(exp(max)/exp(min)),
    diff = exp(max)-exp(min)
    ) %>% 
  
  # Summarize with 95% credible intervals
  select(-iter,-min,-max) %>% 
  group_by(across(-c(prop_change,diff))) %>% 
  summarize(
    across(
      .cols = c(prop_change, diff),
      .fns = list(
        md = ~ quantile(.x, probs = 0.500, na.rm = TRUE),
        lo = ~ quantile(.x, probs = 0.025, na.rm = TRUE),
        up = ~ quantile(.x, probs = 0.975, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) %>% 
  
  # Indicate if credible intervals include zero
  mutate(
    overlap0 = case_when(
      change_lo*change_up >0 ~ F,
      T ~ T
    )
    )

message("[", Sys.time(), "] Proportion change and difference estimated") 
  


# Exporting  ===================================================================

message("============================================================")
message("Start time: ", Sys.time())
message("All extractions and predictions finished")
message("Start exporting section")
message("============================================================")

### Create file names  ###
message("[", Sys.time(), "] Creating master file name scheme... ")

master_name <- paste0("_",sp,"_",response,".rds")

message("[", Sys.time(), "] Master file name scheme set as: ", master_name)


### Export each file  ###
message("[", Sys.time(), "] Starting export...")

out_data <- list(
  gamma_coef,
  theta_coef,
  theta_pred,
  response_pred,
  prop_change
)

name_prefix <- c(
  "coef_gamma",
  "coef_theta",
  "predictions_theta",
  "predictions_response",
  "predictions_prop_change"
)

lapply(seq_along(out_data), function(x) {

  file_name <- file.path(export_dir,paste0(name_prefix[x],master_name))
  saveRDS(out_data[[x]],file_name)
  
  message("[", Sys.time(), "] ", file_name, "exported")
  }
  )
message("[", Sys.time(), "] Exports finished")


message("============================================================")
message("Finished .R run")
message("End time: ", Sys.time())
message("============================================================")

