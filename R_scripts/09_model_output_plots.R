
# House keeping  ---------------------------------------------------------------
rm(list = ls())

# Packages
library(abind)

# Directories
out_dir <- "stan_outputs"
data_dir <- "hpc/data"
plot_dir <- "figures"


# Load in model out puts
# out_files <- list.files(out_dir,"M100.rds$")
out_files <- list.files(out_dir,".rds$")
out_list <- lapply(out_files, function(x) readRDS(file.path(out_dir,x)))
names(out_list) <- gsub("_stan_out_M100.rds|_stan_out.rds","",out_files)

# load in data
data_files <- list.files(data_dir)
data_list <- lapply(data_files, function(x) readRDS(file.path(data_dir,x)))
names(data_list) <- gsub("_input_data.rds","",data_files)


# Create gamma bridge  ---------------------------------------------------------

mods <- names(out_list)
# summarie model coefifcents across all species
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

mods <- names(out_list)
# summarie model coefifcents across all species
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

# Beta plots   -----------------------------------------------------------------

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
z_list <- lapply(mods, function(m){
  z_data <- data_list[[m]]$z
  ### TEMP. FIX THIS DATA PREP SCRIPT
  dimnames(z_data)[[1]] <- c("SRS","TSL","WCA")
  dimnames(z_data)[[2]] <- 1995:2023
  
  # CHange into a single data.frame with columns for region and wateryear
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
predicted_slopes <- lapply(mods, function(m){
  
  # Extract species and response
  sp <- stringr::str_split_fixed(m, "_", n = 2)[1]
  response <- stringr::str_split_fixed(m, "_", n = 2)[2]
  
  # Extract model coefficients
  out <- out_list[[m]]
  gamma_draws <- out$draws("gamma",format = "draws_matrix")
  n_draws <- nrow(gamma_draws)
  
  #
  var_combo <- gamma_bridge[[m]] %>% 
    distinct(z_var,x_var) %>% 
    filter(z_var != "int")
  
  lapply(1:nrow(var_combo), function(v) {
    vars <- var_combo[v,]
    predictor <- vars$z_var
    z_data <- z_list[[m]]
    z_mat <- matrix(
      0,
      nrow = pred_len,
      ncol = ncol(z_data %>% select(-region,-wateryear))
    )
    colnames(z_mat) <- colnames(z_data %>% select(-region,-wateryear))
    
    z <- z_data[,predictor]
    z_mat[,predictor] <- seq(min(z),max(z),length.out = pred_len)
    z_mat[,"int"] <- 1
    
    pred_vec <- NULL
    for(i in seq_len(n_iter)) {
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
      pred <- z_mat %*% gamma
      pred_vec <- cbind(pred_vec,pred)
    }
    pred_df <- data.frame(
      response = response,
      species = sp,
      x_var = vars$x_var,
      z_var = predictor,
      z_range=z_mat[,predictor],
      pred_md = apply(pred_vec, 1, quantile, probs=0.5),
      pred_up = apply(pred_vec, 1, quantile, probs=0.975),
      pred_lo = apply(pred_vec, 1, quantile, probs=0.025)
      
    )
  }
  ) #%>% 
  #bind_rows()
}) %>% 
  bind_rows()


# Gamma plots  -----------------------------------------------------------------  
gamma_plots <- lapply(mods,function(m){ 
  
  # Load in parametes and data
  z_data <- z_list[[m]]
  beta <- beta_list[[m]]
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
    z_df <- z_data
    z_df[,"z"] <- z_data[[predictor]]
    
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
      aes(x = z_range,y = pred_md)
    )+
      geom_smooth(color="black")+
      geom_ribbon(
        aes(ymin = pred_lo,ymax =pred_up),
        alpha=0.2,
        fill = "grey"
      )+
      geom_point(
        data=slope_df,
        aes(x=z,y=mean, color = overlap0),
        size = 2.5
      )+
      geom_errorbar(
        data=slope_df,
        aes(x=z,ymin =lwr,ymax = upr,color = overlap0),
        linewidth = 1,
        inherit.aes = FALSE
      )+
      scale_color_manual(values = c("black","grey"))+
      geom_abline(slope = 0,intercept = 0, color = "red",linewidth = 1.5)+
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
plotter("all","production_mean","pisc_index","plt_cov_int")
plotter("all","production_mean","wet_sum_365day","depth")


# Coefficient plots  ------------------------------------------------------------
mods <- names(out_list)
# summarie model coefifcents across all species
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

# Plot
sp_colors <- 
  rev(c("black","#b5a331","#339d38","#c26a77","#8c6d3f","#2f2585","#2b695c"))
pd <- position_dodge(width = 0.5)
upr_preds <- unique(coef_df$z_var)
upr_preds <- expand.grid(unique(coef_df$z_var),unique(coef_df$response))
lapply(1:nrow(upr_preds),function(i){
  u_p <- upr_preds[i,1]
  r <- upr_preds[i,2]
  plot_df <-coef_df %>% filter(z_var == u_p,response == r)
  if(u_p == "int") plot_df <-plot_df %>% filter(coef != "int")
  
  plot <- ggplot(
    data=plot_df, 
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
    # labs(title = paste(r,u_p,sep="_"))
  # print(plot)
  plot_name <- paste0("coef_tree_plot_",r,"_",u_p,".png")
  plot_height <- 10
  if(u_p != "int") plot_height <- plot_height*1.25
  ggsave(
    file.path(
      plot_dir,
      "coef_tree",
      plot_name
      ),
    plot = plot,
    bg = "transparent",
    width = 5,
    height = plot_height,
    dpi = 300
  )
})



# Hurdle Coefficient plots  ------------------------------------------------------------
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


# Gamma [1,] plots  -----------------------------------------------------------------
pred_len <- 100
n_iter <- 1000
preds <- lapply(mods, function(m){
  
  # Extract species and response
  sp <- stringr::str_split_fixed(m, "_", n = 2)[1]
  response <- stringr::str_split_fixed(m, "_", n = 2)[2]
  
  # Extract model coefficients
  out <- out_list[[m]]
  gamma_draws <- out$draws("gamma",format = "draws_matrix")
  n_draws <- nrow(gamma_draws)
  
  #
  vars <- gamma_bridge[[m]] %>% 
    distinct(z_var,x_var) %>% 
    filter(
      z_var == "int",
      x_var != "int"
    ) %>% 
    pull(x_var)
  
  pred_list <- lapply(1:length(vars), function(v) {
    predictor <- vars[v]
    x_data <- data_list[[m]]$x
    x_mat <- matrix(
      0,
      nrow = pred_len,
      ncol = ncol(x_data)
    )
    colnames(x_mat) <- colnames(x_data)
    
    x <- x_data[,predictor]
    x_mat[,predictor] <- seq(min(x),max(x),length.out = pred_len)
    x_mat[,"int"] <- 1
    
    pred_vec <- NULL
    for(i in seq_len(n_iter)) {
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
        filter(z_var == "int") %>% 
        pull(coef)
      pred <- exp(x_mat %*% gamma)
      pred_vec <- cbind(pred_vec,pred)
    }
    pred_df <- data.frame(
      pred_md = apply(pred_vec, 1, quantile, probs=0.5),
      pred_up = apply(pred_vec, 1, quantile, probs=0.975),
      pred_lo = apply(pred_vec, 1, quantile, probs=0.025),
      x_range=x_mat[,predictor]
    )
    
    pred_df$sp <- sp
    pred_df$predictor <- predictor
    pred_df$response <- response
    pred_df
    
  })
  names(pred_list) <- vars
  a <- bind_rows(pred_list)
}
) %>% 
  bind_rows()

lapply(unique(preds$predictor), function(p){
  r <- "ptob"
  plot_df <- preds %>% 
    filter(
      predictor == p,
      response == r) 
  
  plot <- ggplot(
    plot_df,
    aes(
      x = x_range,
      y = pred_md,
      group = sp,
      colour = sp,
      fill = sp
    )
  )+
    geom_smooth()+
    # geom_ribbon(
    #   aes(ymin = pred_lo,ymax =pred_up,group = sp),
    #   alpha=0.2,
    # )+
    xlab(p)+
    ylab(r);print(plot)
})




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
  
  # Remove nonsignficant CSI's and ...
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