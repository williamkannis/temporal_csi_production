
# House keeping  ===============================================================
rm(list = ls())   

# Packages
library(abind)
library(dplyr)
library(tidyr)
library(ggplot2)

# Directories
out_dir <- "stan_outputs"
data_dir <- "hpc/data"
plot_dir <- "figures"
prod_dir <- "prod_data"


# Load in model out puts
# out_files <- list.files(out_dir,"M100.rds$")
out_files <- list.files(out_dir,".rds$")
out_list <- lapply(out_files, function(x) readRDS(file.path(out_dir,x)))
names(out_list) <- gsub("_stan_out_M100.rds|_stan_out.rds","",out_files)

# Load in data
data_files <- list.files(data_dir)
data_list <- lapply(data_files, function(x) readRDS(file.path(data_dir,x)))
names(data_list) <- gsub("_input_data.rds","",data_files)

# Load in bridges
x_bridge_list <- readRDS(file.path(prod_dir,"csi_model_x_bridge.rds"))
z_bridge_list <- readRDS(file.path(prod_dir,"csi_model_z_bridge.rds"))

# Model names
mods <- names(out_list)


# Data preparation  ============================================================


# Create gamma bridge  ---------------------------------------------------------

# Link gamma parameters with predictor names
gamma_bridge <- lapply(mods, function(m){

  # Extract var names
  z_name <- data.frame(
    z_idx = seq_len(ncol(data_list[[m]]$z[1,,])),
    z_var = colnames(data_list[[m]]$z[1,,])
  )
  x_name <- data.frame(
    x_idx = seq_len(ncol(data_list[[m]]$x)),
    x_var = colnames(data_list[[m]]$x)
  )
  
  
  # Extract model coefficients
  out <- out_list[[m]]
  
  out$summary("gamma",mean) %>% 
  
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
    left_join(z_name, by = join_by(z_idx)) %>% 
    left_join(x_name, by = join_by(x_idx)) %>% 
    select(variable,z_var,x_var)
}
)
names(gamma_bridge) <- names(out_list)

# Create beta bridge  ---------------------------------------------------------

# link region:year slopes with sample info
beta_bridge <- lapply(mods, function(m){
  
  # Extract var names
  x_name <- data.frame(
    x_idx = seq_len(ncol(data_list[[m]]$x)),
    x_var = colnames(data_list[[m]]$x)
  )
  r_name <- data.frame(
    r_idx = 1:3,
    region = c("SRS","TSL","WCA")
  )
  y_name <- data.frame(
    y_idx = 1:29,
    wateryear = 1995:2023
  )
  
  
  # Extract model coefficients
  out <- out_list[[m]]
  
  out$summary("beta",mean) %>% 
    
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
    left_join(x_name, by = join_by(x_idx)) %>% 
    left_join(r_name, by = join_by(r_idx)) %>% 
    left_join(y_name, by = join_by(y_idx)) %>% 
    select(variable,region,wateryear,x_var)
}
)
names(beta_bridge) <- names(out_list)


# Beta list   -----------------------------------------------------------------

# Extract region year slopes
beta_list <- lapply(mods, function(m){
  
  # Extract species and response
  sp <- stringr::str_split_fixed(m, "_", n = 2)[1]
  response <- stringr::str_split_fixed(m, "_", n = 2)[2]
  
  # Extract model coefficients
  out <- out_list[[m]]
  out$summary(
    "beta",
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
      )
    ) %>% 
    left_join(
      beta_bridge[[m]],
      by = join_by(variable)
      )
  
}
)
names(beta_list) <- names(out_list)


# Z data list  -----------------------------------------------------------------

# Extract second level predictor values
z_list <- lapply(mods, function(m){
  z_data <- data_list[[m]]$z
  ### TEMP. FIX THIS DATA PREP SCRIPT
  dimnames(z_data)[[1]] <- c("SRS","TSL","WCA")
  dimnames(z_data)[[2]] <- 1995:2023
  
  # Change into a single data.frame with columns for region and wateryear
  lapply(1:dim(z_data)[[1]], function(x) {
    rg <- dimnames(z_data)[[1]][x]
    as.data.frame(z_data[x,,]) %>% 
      tibble::rownames_to_column("wateryear") %>% 
      mutate(
        region = rg,
        wateryear = as.numeric(wateryear)
        )
    
  }) %>% 
    bind_rows()
}
)
names(z_list) <- names(out_list)


# Predicted slopes  ------------------------------------------------------------
pred_len <- 100
n_iter <- 1000

# Use second level gamma parameters to predict 1st level slopes across predictors
predicted_slopes <- lapply(mods, function(m){
  
  # Extract species and response
  sp <- stringr::str_split_fixed(m, "_", n = 2)[1]
  response <- stringr::str_split_fixed(m, "_", n = 2)[2]
  
  # Extract model coefficients and data
  out <- out_list[[m]]
  gamma_draws <- out$draws("gamma",format = "draws_matrix")
  n_draws <- nrow(gamma_draws)
  z_data <- z_list[[m]]
  z_bridge <- z_bridge_list[[m]]
  
  # Find all combinations of variables
  var_combo <- gamma_bridge[[m]] %>% 
    distinct(z_var,x_var) %>% 
    filter(z_var != "int")
  
  lapply(1:nrow(var_combo), function(v) {
    vars <- var_combo[v,]
    predictor <- vars$z_var
    
    # Create z predictor inputs
    z_mat <- matrix(
      0,
      nrow = pred_len,
      ncol = ncol(z_data %>% select(-region,-wateryear))
    )
    colnames(z_mat) <- colnames(z_data %>% select(-region,-wateryear))
    z <- z_data[,predictor]
    z_mat[,predictor] <- seq(min(z),max(z),length.out = pred_len)
    z_mat[,"int"] <- 1
    
    # Prepare raw data for plotting
    z_raw <- z_bridge[[predictor]]
    z_raw_range <- seq(min(z_raw),max(z_raw),length.out = pred_len)
    
    # Create n iterations of beta predictions using gamma params
    pred_vec <- NULL
    for(i in seq_len(n_iter)) {
      
      # Extract a draw of gamma, and format
      iter <- sample(seq_len(n_draws),1,replace = T)
      gamma <- as.data.frame(gamma_draws[iter,]) %>% 
        pivot_longer( 
          cols = everything(), 
          names_to = "variable", 
          values_to = "coef"
        ) %>% 
        left_join(
          gamma_bridge[[m]],
          by = join_by(variable)
        ) %>% 
        filter(x_var == vars$x_var) %>% 
        pull(coef)
      
      # Create prediction
      pred <- z_mat %*% gamma
      pred_vec <- cbind(pred_vec,pred)
    }
    
    # Summarize prediction in data.frame
    pred_df <- data.frame(
      response = response,
      species = sp,
      x_var = vars$x_var,
      z_var = predictor,
      z_range=z_mat[,predictor],
      z_raw_range = z_raw_range,
      pred_md = apply(pred_vec, 1, quantile, probs=0.5),
      pred_up = apply(pred_vec, 1, quantile, probs=0.975),
      pred_lo = apply(pred_vec, 1, quantile, probs=0.025)
    )
  }
  ) 
}) %>% 
  bind_rows()


# Coefficients  ----------------------------------------------------------------

# summarize model coefficients across all species
coef_df <- lapply(mods, function(m){
  
  # Extract species and response
  sp <- stringr::str_split_fixed(m, "_", n = 2)[1]
  response <- stringr::str_split_fixed(m, "_", n = 2)[2]
  
  
  # Extract var names
  z_name <- data.frame(
    z_idx = seq_len(ncol(data_list[[m]]$z[1,,])),
    z_var = colnames(data_list[[m]]$z[1,,])
  )
  x_name <- data.frame(
    x_idx = seq_len(ncol(data_list[[m]]$x)),
    x_var = colnames(data_list[[m]]$x)
  )
  
  
  # Extract model coefficients
  out <- out_list[[m]]
  gamma_out <- out$summary(
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
      species = sp,
      response = response,
      overlap0 = case_when(
        lwr*upr >0 ~ F,
        T ~ T
      )
    ) %>% 
    
    tidytable::separate_wider_regex(
      cols = variable,
      patterns = c(
        ".*",          
        "\\[",         
        z_idx = "\\d+", 
        ",",           
        x_idx = "\\d+", 
        "\\]"          
      )
    ) %>% 
    mutate(
      z_idx = as.numeric(z_idx),
      x_idx = as.numeric(x_idx)
    ) %>% 
    left_join(z_name) %>% 
    left_join(x_name)
  
}
) %>% bind_rows() %>% 
  mutate(
    coef = factor(
      x_var,
      levels= rev(c("int","depth","dsldd_int","plt_cov_int","peri_vol_int"))
    ),
    species = forcats::fct_rev(species)
  )
row.names(coef_df) <- NULL


# Plotting  ====================================================================


# CSI Plots (v4)  --------------------------------------------------------------

# USE THIS ONE. Includes annual effects

# Find only significant CSI's
sig_df <- coef_df %>% 
  filter(
    z_var != "int",
    x_var != "int",
    !overlap0
  ) %>% 
  mutate(sig = T) %>% 
  distinct(
    species,
    response,
    x_var,
    z_var,
    sig
  )

# Population level coefficients, use when no sig CSI
pop_coef <- coef_df %>% 
  filter(
    z_var == "int",
    x_var != "int"
  ) %>% 
  rename(
    overlap0_all = overlap0
  ) %>% 
  select(
    species,
    response,
    x_var,
    z_var,
    mean,
    lwr,
    upr,
    overlap0_all
  ) %>% 
  anti_join(sig_df %>% distinct(species,response,x_var))

# Extract and format annual drivers
yr_vars <-coef_df %>% 
  filter(
    z_var != "int",
    x_var == "int"
  ) %>% 
  mutate(
    x_var = z_var,
    z_var = "int"
  ) %>% 
  rename(
    overlap0_all = overlap0
  ) %>% 
  select(
    species,
    response,
    x_var,
    z_var,
    mean,
    lwr,
    upr,
    overlap0_all
  )

# Extract season slopes at min and max value of annual predictors
minmax_slopes <-predicted_slopes %>% 
  group_by(
    response,
    species,
    x_var,
    z_var
  ) %>% 
  mutate(
    z_range = case_when(
      z_range == min(z_range) ~ "min",
      z_range == max(z_range) ~ "max"
    ),
    overlap0 = case_when(
      pred_up*pred_lo > 0 ~ F,
      T ~T
    )
  ) %>% 
  ungroup() %>% 
  filter(
    !is.na(z_range),
    x_var != "int"
  ) %>% 
  select(-pred_up,-pred_lo) %>% 
  pivot_wider(
    names_from = z_range,
    values_from = c(pred_md,overlap0)
  ) %>% 
  mutate(
    overlap0_all =case_when(
      overlap0_min & overlap0_max ~ T,
      T~F
    )) %>% 
  
  # Remove significant CSI's and ...
  right_join(
    sig_df,
    by = join_by(response,species, x_var, z_var)
  ) %>% 
  # mutate(
  #   pred_md_min = case_when(sig ~ pred_md_min),
  #   pred_md_max = case_when(sig ~ pred_md_max),
  #   overlap0_min = case_when(sig ~ overlap0_min),
  #   overlap0_max = case_when(sig ~ overlap0_max),
  #   overlap0_all = case_when(sig ~ overlap0_all)
  # ) %>% 
  
  # replace with global coef
  bind_rows(pop_coef) %>% 
  
  # Add in annual coef
  bind_rows(yr_vars)

# Plotting parameters
sp_colors <- 
  c("all" = "black",
    "FUNCHR" = "#b5a331",
    "GAMHOL" = "#339d38",
    "HETFOR" = "#c26a77",
    "JORFLO" = "#8c6d3f",
    "LUCGOO" = "#2f2585",
    "POELAT" = "#2b695c"
  )
pd <- position_dodge(width = .8)
responses <- c("sample_den","biomass_mean","production_mean","ptob")

# Create plot for each response
lapply(responses, function(r) {
  
  # Prepare data for plots
  plot_df <- minmax_slopes %>% 
    filter(response == r) %>% 
    mutate(
      x_var = factor(
        x_var,
        levels= rev(
          c(
            "depth",
            "dsldd_int",
            "plt_cov_int",
            "peri_vol_int",
            "wet_sum_365day",
            "pisc_index"
          )
        )
      ),
      species = forcats::fct_rev(species),
      z_var = factor(
        z_var,
        levels = rev(c("int","wet_sum_365day", "pisc_index"))
      ),
      dodge_group = interaction(
        species,
        z_var,
        sep = "_",
        lex.order = TRUE
      )
    )
  
  # Create plot 
  plot <- ggplot(
    data=plot_df, 
    aes(
      y=x_var,
      x=pred_md_max,
      color = species,
      alpha = overlap0_all,
      group = dodge_group,
      linewidth = species,
      size = species
    )
  ) +
    geom_vline(xintercept = 0, color = "red",linewidth =2) + 
    
    # global 95CI
    geom_errorbarh(
      aes(
        xmin = pred_md_min, 
        xmax = pred_md_max,
        linetype = z_var
      ),
      height =0,
      position = pd,
      # linewidth = 2,
    )+
    
    # CSI coef range
    geom_errorbarh(
      aes(
        xmin = upr, 
        xmax = lwr,
      ),
      height =0,
      position = pd,
      # linewidth = 2,
    )+
    
    # Global coef mean
    geom_point(
      aes(
        y=x_var,
        x=mean,
        color = species
      ),
      position = pd,
      # size = 5
    )+ 
    
    # Effect at max annual predictors
    geom_point(
      aes(
        y=x_var,
        x=pred_md_max,
        color = species,
        alpha = overlap0_max
      ),
      position = pd,
      # size = 5,
      shape = 17)+ 
    
    # Effect at min annual predictors
    geom_point(
      aes(
        y=x_var,
        x=pred_md_min,
        color = species,
        alpha = overlap0_min
      ),
      position = pd,
      # size = 5,
      shape = 15
    )+ 
    
    # Set custom plotting parameters
    scale_alpha_manual(
      values = c(`TRUE` = 0.3, `FALSE` = 1),
      guide = "none"
    ) +
    scale_linetype_manual(
      values = c("41","11","solid")
    )+
    scale_linewidth_manual(
      values = rev(c(3,rep(2,6)))
    )+
    scale_size_manual(
      values = rev(c(8,rep(5,6)))
    )+
    scale_color_manual(values =sp_colors)+
    
    # Theme and display settings
    xlab("")+
    theme(legend.position="none")+
    theme(
      axis.title.y = element_blank(),
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      plot.border  = element_blank(),
      axis.line.x = element_line(color = "black", linewidth = 2),
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.background  = element_rect(fill = "transparent", color = NA),
      axis.text.x = element_text(size = 24)
    )#+
  #labs(title = r);print(plot)
  
  # Export plot
  plot_name <- paste0("coef_csi_plot_",r,"_v4.png")
  plot_height <- 13*1.5
  ggsave(
    file.path(
      plot_dir,
      "csi_coef",
      plot_name
    ),
    plot = plot,
    bg = "transparent",
    width = 5,
    height = plot_height,
    dpi = 300
  )
})


# Gamma plots  -----------------------------------------------------------------  
gamma_plots <- lapply(mods,function(m){ 
  
  # Load in parameters and data
  z_data <- z_list[[m]]
  beta <- beta_list[[m]]
  z_bridge <- z_bridge_list[[m]]
  sp <- stringr::str_split_fixed(m, "_", n = 2)[1]
  res <- stringr::str_split_fixed(m, "_", n = 2)[2]
  
  var_combo <- predicted_slopes %>% 
    distinct(z_var,x_var) %>% 
    filter(z_var != "int")
  
  # Create plot for each predictor combo
  plot_list <- lapply(1:nrow(var_combo), function(v) {
    
    
    ### Filter predicted slopes based on desired plot  ###
    vars <- var_combo[v,]
    predictor <- vars$z_var
    pred_df <- predicted_slopes %>% 
      filter(
        species == sp,
        response == res,
        x_var == vars$x_var,
        z_var == vars$z_var
      )
    
    ### Create data frame with actual slope values  ###
    
    # Extract predictor data for region and wateryear
    # z_df <- z_data
    z_df <- z_bridge
    z_df[,"z"] <- z_df[[predictor]]
    
    # Merge predictor data into slope data
    slope_df <- beta %>% 
      filter(x_var == vars$x_var) %>% 
      left_join(
        z_df,
        by = join_by(region, wateryear)
      ) 
    
    #### Create plots  ###
    g_plot <- ggplot(
      data=pred_df,
      aes(
        x = z_raw_range,
        y = pred_md
        )
    )+
      geom_abline(
        slope = 0,
        intercept = 0, 
        color = "black",
        linewidth = 1.5,
        linetype = "dotted"
        )+
      geom_smooth(color="black")+
      geom_ribbon(
        aes(ymin = pred_lo,ymax =pred_up),
        alpha=0.2,
        fill = "grey"
      )+
      geom_point(
        data=slope_df,
        aes(
          x=z,
          y=mean, 
          alpha = overlap0,
          color = mean > 0
          # color = mean
          ),
        size = 2.5
      )+
      geom_errorbar(
        data=slope_df,
        aes(
          x=z,
          ymin =lwr,
          ymax = upr,
          alpha = overlap0,
          color = mean > 0
          # color = mean
          ),
        linewidth = 1,
        inherit.aes = FALSE
      )+
      scale_alpha_manual(
        values = c(`TRUE` = 0.3, `FALSE` = 1),
        guide = "none"
      ) +
      # scale_color_gradient(low = "red", high = "blue")+
      # scale_color_manual(values = c("black","grey"))+
      # scale_color_gradient2(
      #   low = "blue",       # Color for negative values
      #   # mid = "white",      # Color for zero
      #   high = "red",       # Color for positive values
      #   midpoint = 0        # Forces white to sit exactly at 0
      # )+
      scale_color_manual(
        values = c("TRUE" = "blue", "FALSE" = "red")
        ) +
      theme_classic()+
      theme(
        axis.text.x = element_text(size = 18),  
        axis.text.y = element_text(size = 18),
        legend.position = "none",
        panel.border =  element_rect(color = "black", fill = NA, size = 1)
      )+
      # ylab(paste0("Slope of annual ",p,"-production relationship"))+
      # xlab("Wetted days")
      ylab("")+
      xlab("")
    # print(g_plot)
    g_plot
    
    
  })
  names(plot_list) <- paste(var_combo$z_var,var_combo$x_var,sep = "_")
  plot_list
  
}
)
names(gamma_plots) <- mods

species = "LUCGOO";response = "production_mean";z="pisc_index"; x ="plt_cov_int"

plotter <- function(species,response,z,x) {
  outer_name <- paste(species,response,sep="_")
  inner_name <- paste(z,x,sep="_")
  gamma_plots[[outer_name]][[inner_name]]
}
plotter("POELAT","ptob","pisc_index","peri_vol_int")
plotter("JORFLO","production_mean","pisc_index","peri_vol_int")
plotter("FUNCHR","production_mean","pisc_index","depth")
plotter("GAMHOL","ptob","wet_sum_365day","plt_cov_int")
plotter("LUCGOO","ptob","wet_sum_365day","dsldd_int")
plotter("LUCGOO","biomass_mean","wet_sum_365day","dsldd_int")

# These can be in manuscript
plotter("all","production_mean","pisc_index","plt_cov_int")
plotter("all","production_mean","wet_sum_365day","depth")
plotter("all","biomass_mean","pisc_index","plt_cov_int")
plotter("all","biomass_mean","wet_sum_365day","depth")
plotter("all","sample_den","pisc_index","plt_cov_int")
plotter("all","sample_den","wet_sum_365day","depth")
plotter("all","ptob","wet_sum_365day","plt_cov_int")


# First level plot with CSI  ---------------------------------------------------

## THINK IF YOU WANT ON LOG SCALE

pred_len <- 100
n_iter <- 1000

# Create response predictions using slopes at min and max CSI values
pred_list <- lapply(mods, function(m){
  
  # Extract species and response
  sp <- stringr::str_split_fixed(m, "_", n = 2)[1]
  response <- stringr::str_split_fixed(m, "_", n = 2)[2]
  
  # Extract model coefficients and data
  out <- out_list[[m]]
  gamma_draws <- out$draws("gamma",format = "draws_matrix")
  n_draws <- nrow(gamma_draws)
  x_data <- data_list[[m]]$x
  x_bridge <- x_bridge_list[[m]]
  z_data <- z_list[[m]]
  
  ### 2nd level predictions  ###
  
  # All x and z combos
  var_combo <- gamma_bridge[[m]] %>% 
    distinct(z_var,x_var) %>% 
    filter(z_var != "int")
  
  # Extract i iterations of each slope at min and max CSI values
  minmax_list <- lapply(seq_len(n_iter), function(i) {
    iter <- sample(seq_len(n_draws),1,replace = T)
    
    # Create data.frame with beta coef predictions for every x/z  combo
    lapply(1:nrow(var_combo), function(v) {
      
      
      vars <- var_combo[v,]
      predictor <- vars$z_var
      
      # Prepare z predictiorn input data
      z_mat <- matrix(
        0,
        nrow = 3,
        ncol = ncol(z_data %>% select(-region,-wateryear))
      )
      colnames(z_mat) <- colnames(z_data %>% select(-region,-wateryear))
      z <- z_data[,predictor]
      z_mat[,predictor] <- c(min(z),0,max(z))
      z_mat[,"int"] <- 1
      
      # Extract a draw of gamma paramters
      gamma <- as.data.frame(gamma_draws[iter,]) %>% 
        pivot_longer( 
          cols = everything(), 
          names_to = "variable", 
          values_to = "coef"
        ) %>% 
        left_join(
          gamma_bridge[[m]],
          by = join_by(variable)
        ) %>% 
        filter(x_var == vars$x_var) %>% 
        pull(coef)
      
      # Create beta coef preidction
      coef <- z_mat %*% gamma
      data.frame(
        x_var = vars$x_var,
        z_var = vars$z_var,
        z_range = c("min","mean","max"),
        coef =coef
      )
    }
    ) %>% 
      bind_rows()
  }
  )
  
  ### 1st level predictions  ###
  
  # X and Z combos used for plots
  plot_combo <- var_combo %>% 
    filter(x_var !="int") %>% 
    distinct(z_var, x_var) %>% 
    tidyr::crossing(z_range = c("min","mean","max")) %>% 
    as.data.frame()
  
  # Create data.frame with predictions for every variable combo
  lapply(1:nrow(plot_combo), function(p) {
    
    # Extract indicated variables
    x_v <- plot_combo[p,2]
    z_v <- plot_combo[p,1]
    z_r <- plot_combo[p,3]
    
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
    x_raw <- x_bridge[[x_v]]
    x_raw_range <- seq(min(x_raw),max(x_raw),length.out = pred_len)
    
    # Create predictions for every iteration of gamma predictions
    pred_list <- lapply(minmax_list, function(df){
      
      # Extract beta coef
      beta <- df %>% 
        filter(
          z_var == z_v,
          z_range == z_r
        ) %>% 
        pull(coef)
      
      # # Create prediction on normal scale
      # exp(x_mat %*% beta)
      
      # Create prediction on log scale
      x_mat %*% beta
    })
    
    # Prepare output into dataframe with mean and 95CI
    pred_mat <- do.call(cbind,pred_list)
    data.frame(
      species = sp,
      response = response,
      x_var = x_v,
      z_var = z_v,
      z_range = z_r,
      x_range=x_mat[,x_v],
      x_raw_range = x_raw_range,
      pred_md = apply(pred_mat, 1, quantile, probs=0.5),
      pred_up = apply(pred_mat, 1, quantile, probs=0.975),
      pred_lo = apply(pred_mat, 1, quantile, probs=0.025)
    )
  }) %>% 
    bind_rows()
  
}
)
names(pred_list) <- mods


pred_plots <- lapply(mods, function(m){
  response <- stringr::str_split_fixed(m, "_", n = 2)[2]
  pred_df <- pred_list[[m]] %>% 
    
    # Remove non-significance CSIs
    inner_join(
      sig_df,
      by = join_by(species, response, x_var, z_var)
    ) %>% 
    filter(z_range != "mean")
  
  # Change one set of mean parameters (doesn't matter which, both are equal) 
  # to intercept. These will be used in absence of sig CSI
  for_df <- pred_list[[m]] %>% 
    filter(
      z_var == "pisc_index",
      z_range == "mean"
    ) %>% 
    mutate(z_var = "int") %>% 
    anti_join(
      sig_df,
      by = join_by(species, response, x_var)
    ) %>% 
    bind_rows(pred_df)
  
  # Set maximum y axis values
  y_max <- max(for_df$pred_md, for_df$pred_lo, for_df$pred_up)*1.01
  # y_max <- ceiling(y_max/5)*5
  y_min <- min(for_df$pred_md, for_df$pred_lo, for_df$pred_up)*1.01
  # y_min <- floor(y_min/5)*5
  
  # Create plots for each combination of variables
  x_vars <- unique(for_df$x_var)
  
  plot_list <- lapply(x_vars, function(x_v){
    
    # Filter data
    plot_df <- for_df %>% 
      filter(
        x_var == x_v
      ) %>% 
      
      # Create grouping variable based on z variable and range
      mutate(
        group_idx = paste(z_range,z_var,sep="_")
      )
    
    # Set colors based on z_range and var
     z_colors <- c(
      "mean_int" = "black",
      "min_wet_sum_365day" = "#FFD700",
      "max_wet_sum_365day" = "#1B4F72",
      "max_pisc_index"    = "#CC79A7",
      "min_pisc_index"   = "#009E73"
      )

    # Plotting
    plot <- ggplot(
      plot_df,
      aes(
        x = x_raw_range,
        y = pred_md,
        colour = group_idx,
        fill = group_idx,
        group = group_idx
      ),
    ) +
      
      # Predictions and credible intervals
      geom_smooth()+
      geom_ribbon(
        aes(
          ymin = pred_lo,
          ymax = pred_up
        ),
        alpha=0.2,
        colour = NA
      )+
      
      # Axes and other formatting
      ylim(y_min,y_max) +
      ylab("")+
      xlab("")+
      scale_color_manual(values = z_colors) +
      scale_fill_manual(values = z_colors) +
      theme_classic()+
      theme(
        legend.position = "none",
        panel.border =  element_rect(color = "black", fill = NA, size = 1),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = margin(5, 18, 0, 0)
      )
    
    # Reduce the number of breaks for dsd and peri
    if(x_v %in% c("dsldd_int", "peri_vol_int")){
      plot <- plot + scale_x_continuous(n.breaks = 3)
    }

    # Remove labels for plots that will be inside grid
    if(x_v == "depth"){
      plot <- plot + theme(axis.text.y = element_text(size = 18))
    } else {
      plot <- plot + theme(axis.text.y  = element_text(size = 0))
    }
    if(response == "ptob"){
      plot <- plot + theme(axis.text.x = element_text(size = 18))
    } else{
      plot <- plot + theme(axis.text.x  = element_text(size = 0))
    } 
    plot
  })
  names(plot_list) <- x_vars
  plot_list
})
names(pred_plots) <- mods

# Make into grid for each species
species <- unique(stringr::str_split_fixed(mods, "_", n = 2)[,1])
lapply(species, function(s){
  
  # ordering for grid
  resp_order <- c(
    "sample_den",
    "biomass_mean",
    "production_mean",
    "ptob"
  )
  pred_order <- c(
    "depth",
    "dsldd_int",
    "plt_cov_int",
    "peri_vol_int"
  )
  sp_resp_order <- sapply(
    resp_order, 
    function(x) paste(s,x,sep = "_")
    )
  
  # Extract plots and reorder
  plot_list <- unlist(pred_plots[sp_resp_order])
  order <- apply(expand.grid(sp_resp_order, pred_order), 1, paste, collapse = ".")
  plots_ordered <-plot_list[order]
  
  # Plot
  plot <- cowplot::plot_grid(
    plotlist = plots_ordered,
    ncol = 4,
    nrow = 4,
    align = "hv",
    axis = "none",
    byrow = FALSE
  )
  
  # Export
  plot_name <- paste0("response_pred_",s,".png")
  
  ggsave(
    filename = file.path(
      plot_dir,
      "response_prediction",
      plot_name
      ),
    plot = plot,
    width = 12,
    height = 8,
    dpi = 600
  )
  
})


# Hurdle Coefficient plots (make appendix)  ------------------------------------
h_mods <- mods[grepl("ptob",mods)]
# summarie model coefifcents across all species
hurdle_coef_df <- lapply(h_mods, function(m){
  
  # Extract species and response
  sp <- stringr::str_split_fixed(m, "_", n = 2)[1]
  response <- stringr::str_split_fixed(m, "_", n = 2)[2]
  
  
  # Extract var names
  # z_name <- data.frame(
  #   z_idx = seq_len(ncol(data_list[[m]]$z[1,,])),
  #   z_var = colnames(data_list[[m]]$z[1,,])
  # )
  x_name <- data.frame(
    x_idx = seq_len(ncol(data_list[[m]]$x)),
    x_var = colnames(data_list[[m]]$x)
  )
  
  
  # Extract model coefficients
  out <- out_list[[m]]
  beta_out <- out$summary(
    "beta_hurdle",
    mean,
    quantile, 
    .args = list(probs = c(0.025, 0.975))
  ) %>% 
    rename(
      lwr = `2.5%`,
      upr = `97.5%`
    ) %>% 
    mutate(
      species = sp,
      response = response,
      overlap0 = case_when(
        lwr*upr >0 ~ F,
        T ~ T
      )
    ) %>% 
    
    tidytable::separate_wider_regex(
      cols = variable,
      patterns = c(
        ".*",          
        "\\[",         
        x_idx = "\\d+", 
        "\\]"          
      )
    ) %>% 
    mutate(
      x_idx = as.numeric(x_idx)
    ) %>% 
    left_join(x_name, by = join_by(x_idx))
  
}
) %>% bind_rows() %>% 
  mutate(
    coef = factor(
      x_var,
      levels= rev(c("int","depth","dsldd_int","plt_cov_int","peri_vol_int"))
    ),
    species = forcats::fct_rev(species)
  )
row.names(hurdle_coef_df) <- NULL

# Plot
sp_colors <- 
  rev(c("black","#b5a331","#339d38","#c26a77","#8c6d3f","#2f2585","#2b695c"))
pd <- position_dodge(width = 0.5)

h_plot_df <-hurdle_coef_df %>%
  filter(coef != "int")

h_plot <- ggplot(
  data=h_plot_df, 
  aes(
    y=coef,
    x=mean,
    color = species,
    alpha = overlap0)
) +
  geom_vline(xintercept = 0, color = "red",linewidth =2) + 
  geom_errorbarh(
    aes(xmin = lwr, xmax = upr),
    height =0,
    position = pd,
    linewidth = 2,
  )+
  geom_point(position = pd,size = 5)+ 
  scale_alpha_manual(
    values = c(`TRUE` = 0.3, `FALSE` = 1),
    guide = "none"
  ) +
  scale_color_manual(values =sp_colors)+
  xlab("")+
  theme(legend.position="none")+
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    plot.border  = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 2),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background  = element_rect(fill = "transparent", color = NA),
    axis.text.x = element_text(size = 24)
  ) #+
# labs(title = paste("hurdle_ptob"))
print(h_plot)
plot_name <- paste0("H_coef_tree_plot_ptob.png")
ggsave(
  file.path(plot_dir,plot_name),
  plot = h_plot,
  bg = "transparent",
  width = 5,
  height = 10,
  dpi = 300
)

# # Gamma [1,] plots (old)  ----------------------------------------------------
# pred_len <- 100
# n_iter <- 1000
# preds <- lapply(mods, function(m){
#   
#   # Extract species and response
#   sp <- stringr::str_split_fixed(m, "_", n = 2)[1]
#   response <- stringr::str_split_fixed(m, "_", n = 2)[2]
#   
#   # Extract model coefficients
#   out <- out_list[[m]]
#   gamma_draws <- out$draws("gamma",format = "draws_matrix")
#   n_draws <- nrow(gamma_draws)
#   
#   #
#   vars <- gamma_bridge[[m]] %>% 
#     distinct(z_var,x_var) %>% 
#     filter(
#       z_var == "int",
#       x_var != "int"
#     ) %>% 
#     pull(x_var)
#   
#   pred_list <- lapply(1:length(vars), function(v) {
#     predictor <- vars[v]
#     x_data <- data_list[[m]]$x
#     x_mat <- matrix(
#       0,
#       nrow = pred_len,
#       ncol = ncol(x_data)
#     )
#     colnames(x_mat) <- colnames(x_data)
#     
#     x <- x_data[,predictor]
#     x_mat[,predictor] <- seq(min(x),max(x),length.out = pred_len)
#     x_mat[,"int"] <- 1
#     
#     pred_vec <- NULL
#     for(i in seq_len(n_iter)) {
#       iter <- sample(seq_len(n_draws),1,replace = T)
#       gamma <- as.data.frame(gamma_draws[iter,]) %>% 
#         pivot_longer( 
#           cols = everything(), 
#           names_to = "variable", 
#           values_to = "coef"
#         ) %>% 
#         left_join(
#           gamma_bridge[[m]],
#           by = join_by(variable)
#         ) %>% 
#         filter(z_var == "int") %>% 
#         pull(coef)
#       pred <- exp(x_mat %*% gamma)
#       pred_vec <- cbind(pred_vec,pred)
#     }
#     pred_df <- data.frame(
#       pred_md = apply(pred_vec, 1, quantile, probs=0.5),
#       pred_up = apply(pred_vec, 1, quantile, probs=0.975),
#       pred_lo = apply(pred_vec, 1, quantile, probs=0.025),
#       x_range=x_mat[,predictor]
#     )
#     
#     pred_df$sp <- sp
#     pred_df$predictor <- predictor
#     pred_df$response <- response
#     pred_df
#     
#   })
#   names(pred_list) <- vars
#   a <- bind_rows(pred_list)
# }
# ) %>% 
#   bind_rows()
# 
# lapply(unique(preds$predictor), function(p){
#   r <- "ptob"
#   plot_df <- preds %>% 
#     filter(
#       predictor == p,
#       response == r) 
#   
#   plot <- ggplot(
#     plot_df,
#     aes(
#       x = x_range,
#       y = pred_md,
#       group = sp,
#       colour = sp,
#       fill = sp
#     )
#   )+
#     geom_smooth()+
#     # geom_ribbon(
#     #   aes(ymin = pred_lo,ymax =pred_up,group = sp),
#     #   alpha=0.2,
#     # )+
#     xlab(p)+
#     ylab(r);print(plot)
# })
# # Coefficient plots (old)  ------------------------------------------------------
# # Plot
# sp_colors <- 
#   rev(c("black","#b5a331","#339d38","#c26a77","#8c6d3f","#2f2585","#2b695c"))
# pd <- position_dodge(width = 0.5)
# upr_preds <- unique(coef_df$z_var)
# upr_preds <- expand.grid(unique(coef_df$z_var),unique(coef_df$response))
# lapply(1:nrow(upr_preds),function(i){
#   u_p <- upr_preds[i,1]
#   r <- upr_preds[i,2]
#   plot_df <-coef_df %>% filter(z_var == u_p,response == r)
#   if(u_p == "int") plot_df <-plot_df %>% filter(coef != "int")
#   
#   plot <- ggplot(
#     data=plot_df, 
#     aes(
#       y=coef,
#       x=mean,
#       color = species,
#       alpha = overlap0)
#   ) +
#     geom_vline(xintercept = 0, color = "red",linewidth =2) + 
#     geom_errorbarh(
#       aes(xmin = lwr, xmax = upr),
#       height =0,
#       position = pd,
#       linewidth = 2,
#     )+
#     geom_point(position = pd,size = 5)+ 
#     scale_alpha_manual(
#       values = c(`TRUE` = 0.3, `FALSE` = 1),
#       guide = "none"
#     ) +
#     scale_color_manual(values =sp_colors)+
#     xlab("")+
#     theme(legend.position="none")+
#     theme(
#       axis.title.y = element_blank(),
#       axis.text.y  = element_blank(),
#       axis.ticks.y = element_blank(),
#       panel.grid.major.y = element_blank(),
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       panel.border = element_blank(),
#       plot.border  = element_blank(),
#       axis.line.x = element_line(color = "black", linewidth = 2),
#       panel.background = element_rect(fill = "transparent", color = NA),
#       plot.background  = element_rect(fill = "transparent", color = NA),
#       axis.text.x = element_text(size = 24)
#     ) #+
#   # labs(title = paste(r,u_p,sep="_"))
#   # print(plot)
#   plot_name <- paste0("coef_tree_plot_",r,"_",u_p,".png")
#   plot_height <- 10
#   if(u_p != "int") plot_height <- plot_height*1.25
#   ggsave(
#     file.path(
#       plot_dir,
#       "coef_tree",
#       plot_name
#     ),
#     plot = plot,
#     bg = "transparent",
#     width = 5,
#     height = plot_height,
#     dpi = 300
#   )
# })
# # CSI Plots (v1)  --------------------------
# 
# minmax_slopes <-predicted_slopes %>% 
#   group_by(
#     response,
#     species,
#     x_var,
#     z_var
#     ) %>% 
#   mutate(
#     z_range = case_when(
#       z_range == min(z_range) ~ "min",
#       z_range == max(z_range) ~ "max"
#     ),
#     overlap0 = case_when(
#       pred_up*pred_lo > 0 ~ F,
#       T ~T
#     )
#   ) %>% 
#   ungroup() %>% 
#   filter(
#     !is.na(z_range),
#     x_var != "int"
#     )
# 
# sp_colors <- 
#   rev(c("black","#b5a331","#339d38","#c26a77","#8c6d3f","#2f2585","#2b695c"))
# pd <- position_dodge(width = 0.5)
# 
# responses <- c("sample_den","biomass_mean","production_mean","ptob")
# lapply(responses, function(r) {
#   plot_df <- minmax_slopes %>% 
#     filter(response == r, z_var == "wet_sum_365day")
#   
#   plot <- ggplot(
#     data=plot_df, 
#     aes(
#       y=x_var,
#       x=pred_md,
#       color = species,
#       alpha = overlap0,
#       shape = factor(z_range))
#   ) +
#     geom_vline(xintercept = 0, color = "red",linewidth =2) + 
#     geom_errorbarh(
#       aes(xmin = pred_lo, xmax = pred_up),
#       height =0,
#       position = pd,
#       linewidth = 2,
#     )+
#     geom_point(position = pd,size = 5)+ 
#     scale_alpha_manual(
#       values = c(`TRUE` = 0.3, `FALSE` = 1),
#       guide = "none"
#     ) +
#     scale_color_manual(values =sp_colors)+
#     xlab("")+
#     theme(legend.position="none")+
#     theme(
#       axis.title.y = element_blank(),
#       # axis.text.y  = element_blank(),
#       axis.ticks.y = element_blank(),
#       panel.grid.major.y = element_blank(),
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       panel.border = element_blank(),
#       plot.border  = element_blank(),
#       axis.line.x = element_line(color = "black", linewidth = 2),
#       panel.background = element_rect(fill = "transparent", color = NA),
#       plot.background  = element_rect(fill = "transparent", color = NA),
#       axis.text.x = element_text(size = 24)
#     );print(plot)
# })

# # CSI Plots (v2)  --------------------------
# 
# # Find only significant interactions
# sig_df <- coef_df %>% 
#   filter(
#     z_var != "int",
#     x_var != "int",
#     !overlap0
#     ) %>% 
#   mutate(sig = T) %>% 
#   distinct(
#     species,
#     response,
#     x_var,
#     z_var,
#     sig
#     )
# 
# 
# minmax_slopes <-predicted_slopes %>% 
#   group_by(
#     response,
#     species,
#     x_var,
#     z_var
#   ) %>% 
#   mutate(
#     z_range = case_when(
#       z_range == min(z_range) ~ "min",
#       z_range == max(z_range) ~ "max"
#     ),
#     overlap0 = case_when(
#       pred_up*pred_lo > 0 ~ F,
#       T ~T
#     )
#   ) %>% 
#   ungroup() %>% 
#   filter(
#     !is.na(z_range),
#     x_var != "int"
#   ) %>% 
#   select(-pred_up,-pred_lo) %>% 
#   pivot_wider(
#     names_from = z_range,
#     values_from = c(pred_md,overlap0)
#   ) %>% 
#   mutate(
#     overlap0_all =case_when(
#       overlap0_min & overlap0_max ~ T,
#       T~F
#     )) %>% 
#   left_join(
#     sig_df,
#     by = join_by(response,species, x_var, z_var)
#   ) %>% 
#   mutate(
#     pred_md_min = case_when(sig ~ pred_md_min),
#     pred_md_max = case_when(sig ~ pred_md_max),
#     overlap0_min = case_when(sig ~ overlap0_min),
#     overlap0_max = case_when(sig ~ overlap0_max),
#     overlap0_all = case_when(sig ~ overlap0_all)
#   )
# 
# sp_colors <- 
#   c("all" = "black",
#     "FUNCHR" = "#b5a331",
#     "GAMHOL" = "#339d38",
#     "HETFOR" = "#c26a77",
#     "JORFLO" = "#8c6d3f",
#     "LUCGOO" = "#2f2585",
#     "POELAT" = "#2b695c"
#     )
# 
# # sp_colors <- 
# #   rev(c("black","#b5a331","#339d38","#c26a77","#8c6d3f","#2f2585","#2b695c"))
# pd <- position_dodge(width = .8)
# 
# responses <- c("sample_den","biomass_mean","production_mean","ptob")
# lapply(responses, function(r) {
#   
# 
#   
#   plot_df <- minmax_slopes %>% 
#     filter(response == r #,
#           # z_var == "pisc_index"
#            # z_var == "wet_sum_365day"
#            ) %>% 
#     mutate(
#       x_var = factor(
#         x_var,
#         levels= rev(c("int","depth","dsldd_int","plt_cov_int","peri_vol_int"))
#       ),
#       species = forcats::fct_rev(species),
#       z_var = factor(
#         z_var,
#         levels = c("wet_sum_365day", "pisc_index")
#       ),
#       dodge_group = interaction(
#         species,
#         z_var,
#         sep = "_",
#         lex.order = TRUE
#       )
#     )
#   
#   
#   plot <- ggplot(
#     data=plot_df, 
#     aes(
#       y=x_var,
#       x=pred_md_max,
#       color = species,
#       alpha = overlap0_all,
#       group = dodge_group
#       )
#   ) +
#     geom_vline(xintercept = 0, color = "red",linewidth =2) + 
#     geom_errorbarh(
#       aes(
#         xmin = pred_md_min, 
#         xmax = pred_md_max,
#         linetype = z_var
#         ),
#       height =0,
#       position = pd,
#       linewidth = 2,
#     )+
#     geom_point(
#       aes(
#         y=x_var,
#         x=pred_md_max,
#         color = species,
#         alpha = overlap0_max)
#       ,
#       position = pd,
#       size = 5,
#       shape = 17)+ 
#     geom_point(
#       aes(
#         y=x_var,
#         x=pred_md_min,
#         color = species,
#         alpha = overlap0_min
#         ),
#       position = pd,
#       size = 5,
#       shape = 15
#       )+ 
#     scale_alpha_manual(
#       values = c(`TRUE` = 0.3, `FALSE` = 1),
#       guide = "none"
#     ) +
#     scale_color_manual(values =sp_colors)+
#     xlab("")+
#     theme(legend.position="none")+
#     theme(
#       axis.title.y = element_blank(),
#       axis.text.y  = element_blank(),
#       axis.ticks.y = element_blank(),
#       panel.grid.major.y = element_blank(),
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       panel.border = element_blank(),
#       plot.border  = element_blank(),
#       axis.line.x = element_line(color = "black", linewidth = 2),
#       panel.background = element_rect(fill = "transparent", color = NA),
#       plot.background  = element_rect(fill = "transparent", color = NA),
#       axis.text.x = element_text(size = 24)
#     )#+
#     # labs(title = r);#print(plot)
#   
#   plot_name <- paste0("coef_csi_plot_",r,".png")
#   plot_height <- 10
#   ggsave(
#     file.path(
#       plot_dir,
#       "csi_coef",
#       plot_name
#     ),
#     plot = plot,
#     bg = "transparent",
#     width = 5,
#     height = plot_height,
#     dpi = 300
#   )
# })


# # CSI Plots (v3)  --------------------------
# 
# # Find only significant interactions
# sig_df <- coef_df %>% 
#   filter(
#     z_var != "int",
#     x_var != "int",
#     !overlap0
#   ) %>% 
#   mutate(sig = T) %>% 
#   distinct(
#     species,
#     response,
#     x_var,
#     z_var,
#     sig
#   )
# 
# # Population level coefficents
# pop_coef <- coef_df %>% 
#   filter(
#     z_var == "int",
#     x_var != "int"
#   ) %>% 
#   rename(
#     overlap0_all = overlap0
#   ) %>% 
#   select(
#     species,
#     response,
#     x_var,
#     z_var,
#     mean,
#     lwr,
#     upr,
#     overlap0_all
#   ) %>% 
#   anti_join(sig_df %>% distinct(species,response,x_var))
# 
# 
# minmax_slopes <-predicted_slopes %>% 
#   group_by(
#     response,
#     species,
#     x_var,
#     z_var
#   ) %>% 
#   mutate(
#     z_range = case_when(
#       z_range == min(z_range) ~ "min",
#       z_range == max(z_range) ~ "max"
#     ),
#     overlap0 = case_when(
#       pred_up*pred_lo > 0 ~ F,
#       T ~T
#     )
#   ) %>% 
#   ungroup() %>% 
#   filter(
#     !is.na(z_range),
#     x_var != "int"
#   ) %>% 
#   select(-pred_up,-pred_lo) %>% 
#   pivot_wider(
#     names_from = z_range,
#     values_from = c(pred_md,overlap0)
#   ) %>% 
#   mutate(
#     overlap0_all =case_when(
#       overlap0_min & overlap0_max ~ T,
#       T~F
#     )) %>% 
#   right_join(
#     sig_df,
#     by = join_by(response,species, x_var, z_var)
#   ) %>% 
#   # mutate(
#   #   pred_md_min = case_when(sig ~ pred_md_min),
#   #   pred_md_max = case_when(sig ~ pred_md_max),
#   #   overlap0_min = case_when(sig ~ overlap0_min),
#   #   overlap0_max = case_when(sig ~ overlap0_max),
#   #   overlap0_all = case_when(sig ~ overlap0_all)
#   # ) %>% 
#   bind_rows(pop_coef)
# 
# sp_colors <- 
#   c("all" = "black",
#     "FUNCHR" = "#b5a331",
#     "GAMHOL" = "#339d38",
#     "HETFOR" = "#c26a77",
#     "JORFLO" = "#8c6d3f",
#     "LUCGOO" = "#2f2585",
#     "POELAT" = "#2b695c"
#   )
# 
# # sp_colors <- 
# #   rev(c("black","#b5a331","#339d38","#c26a77","#8c6d3f","#2f2585","#2b695c"))
# pd <- position_dodge(width = .8)
# 
# responses <- c("sample_den","biomass_mean","production_mean","ptob")
# lapply(responses, function(r) {
#   
#   
#   
#   plot_df <- minmax_slopes %>% 
#     filter(response == r) %>% 
#     mutate(
#       x_var = factor(
#         x_var,
#         levels= rev(c("int","depth","dsldd_int","plt_cov_int","peri_vol_int"))
#       ),
#       species = forcats::fct_rev(species),
#       z_var = factor(
#         z_var,
#         levels = rev(c("int","wet_sum_365day", "pisc_index"))
#       ),
#       dodge_group = interaction(
#         species,
#         z_var,
#         sep = "_",
#         lex.order = TRUE
#       )
#     )
#   
#   
#   plot <- ggplot(
#     data=plot_df, 
#     aes(
#       y=x_var,
#       x=pred_md_max,
#       color = species,
#       alpha = overlap0_all,
#       group = dodge_group,
#       linewidth = species,
#       size = species
#     )
#   ) +
#     geom_vline(xintercept = 0, color = "red",linewidth =2) + 
#     geom_errorbarh(
#       aes(
#         xmin = upr, 
#         xmax = lwr,
#       ),
#       height =0,
#       position = pd,
#       # linewidth = 2,
#     )+
#     geom_errorbarh(
#       aes(
#         xmin = pred_md_min, 
#         xmax = pred_md_max,
#         linetype = z_var
#       ),
#       height =0,
#       position = pd,
#       # linewidth = 2,
#     )+
#     geom_point(
#       aes(
#         y=x_var,
#         x=mean,
#         color = species
#       ),
#       position = pd,
#       # size = 5
#     )+ 
#     geom_point(
#       aes(
#         y=x_var,
#         x=pred_md_max,
#         color = species,
#         alpha = overlap0_max)
#       ,
#       position = pd,
#       # size = 5,
#       shape = 17)+ 
#     geom_point(
#       aes(
#         y=x_var,
#         x=pred_md_min,
#         color = species,
#         alpha = overlap0_min
#       ),
#       position = pd,
#       # size = 5,
#       shape = 15
#     )+ 
#     scale_alpha_manual(
#       values = c(`TRUE` = 0.3, `FALSE` = 1),
#       guide = "none"
#     ) +
#     scale_linetype_manual(values = c("41","11","solid"))+
#     scale_linewidth_manual(values = rev(c(3,rep(2,6))))+
#     scale_size_manual(values = rev(c(8,rep(5,6))))+
#     scale_color_manual(values =sp_colors)+
#     xlab("")+
#     theme(legend.position="none")+
#     theme(
#       axis.title.y = element_blank(),
#       axis.text.y  = element_blank(),
#       axis.ticks.y = element_blank(),
#       panel.grid.major.y = element_blank(),
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       panel.border = element_blank(),
#       plot.border  = element_blank(),
#       axis.line.x = element_line(color = "black", linewidth = 2),
#       panel.background = element_rect(fill = "transparent", color = NA),
#       plot.background  = element_rect(fill = "transparent", color = NA),
#       axis.text.x = element_text(size = 24)
#     )#+
#   #labs(title = r);print(plot)
#   
#   plot_name <- paste0("coef_csi_plot_",r,"_v3.png")
#   plot_height <- 13
#   ggsave(
#     file.path(
#       plot_dir,
#       "csi_coef",
#       plot_name
#     ),
#     plot = plot,
#     bg = "transparent",
#     width = 5,
#     height = plot_height,
#     dpi = 300
#   )
# })
# # Pred plot (old) ------------------------------------------------------------
# pred_plots <- lapply(mods, function(m){
#   response <- stringr::str_split_fixed(m, "_", n = 2)[2]
#   pred_df <- pred_list[[m]] %>% 
#     
#     # Remove unsignifanct CSIs
#     inner_join(
#       sig_df,
#       by = join_by(species, response, x_var, z_var)
#     ) %>% 
#     filter(z_range != "mean")
#   
#   # Change one set of mean paramters (doesn't matter which, both are equal) 
#   # to intercept. THese will be used in absence of sig CSI
#   for_df <- pred_list[[m]] %>% 
#     filter(
#       z_var == "pisc_index",
#       z_range == "mean"
#     ) %>% 
#     mutate(z_var = "int") %>% 
#     anti_join(
#       sig_df,
#       by = join_by(species, response, x_var)
#     ) %>% 
#     bind_rows(pred_df)
#   
#   # Set maximum y axis values
#   y_max <- max(for_df$pred_md, for_df$pred_lo, for_df$pred_up)*1.01
#   y_max <- ceiling(y_max/5)*5
#   
#   # Create plots for each combination of variables
#   var_combo <- for_df %>% 
#     distinct(x_var,z_var)
#   
#   plot_list <- lapply(1:nrow(var_combo), function(v){
#     
#     # Extract indicated variables
#     x_v <- var_combo[v,"x_var"]
#     z_v <- var_combo[v,"z_var"]
#     
#     # Filter data 
#     plot_df <- for_df %>% 
#       filter(
#         x_var == x_v,
#         z_var == z_v
#       )
#     
#     if (z_v == "int") z_colors <- "black"
#     if (z_v == "wet_sum_365day") z_colors <- c(
#       "min" = "#FFD700",
#       "mean" = "#1ABC9C",
#       "max" = "#1B4F72"
#     )
#     
#     if (z_v == "pisc_index") z_colors <- c(
#       "max"    = "#E69F00",
#       "mean" = "#D55E00",
#       "min"   = "#0072B2"
#     )
#     
#     plot <- ggplot(
#       plot_df,
#       aes(
#         x = x_range,
#         y = log(pred_md),
#         colour = z_range,
#         fill = z_range,
#         group = z_range
#       ),
#     ) +
#       geom_smooth()+
#       geom_ribbon(
#         aes(
#           ymin = log(pred_lo),
#           ymax = log(pred_up)
#         ),
#         alpha=0.2,
#         colour = NA
#       )+
#       ylim(NA,log(y_max)) +
#       ylab("")+
#       xlab("")+
#       scale_color_manual(values = z_colors) +
#       scale_fill_manual(values = z_colors) +
#       theme_classic()+
#       theme(
#         legend.position = "none",
#         panel.border =  element_rect(color = "black", fill = NA, size = 1),
#         axis.title.y = element_blank(),
#         axis.title.x = element_blank(),
#         plot.margin = margin(0, 0, 0, 0)
#       )
#     if(x_v == "depth"){
#       plot <- plot + theme(axis.text.y = element_text(size = 18))
#     } else {
#       plot <- plot + theme(axis.text.y  = element_text(size = 0))
#     }
#     if(response == "ptob"){
#       plot <- plot + theme(axis.text.x = element_text(size = 18))
#     } else{
#       plot <- plot + theme(axis.text.x  = element_text(size = 0))
#     } 
#     plot
#   })
#   # names(plot_list) <- paste(var_combo$z_var,var_combo$x_var,sep = "_")
#   names(plot_list) <- var_combo$x_var
#   plot_list
# })