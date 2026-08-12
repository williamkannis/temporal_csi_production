#-------------------------------------------------------------------------------
#
#   
#
#-------------------------------------------------------------------------------

# AUTHOR: William K. Annis

# CREATED: July 13, 2026

# DESCRIPTION: 


# Housekeeping  ----------------------------------------------------------------
rm(list = ls())

# Load in packages
library(dplyr)
library(purrr)
library(abind)
library(rstan)

# directories
mod_dir <- "stan_scripts"
input_dir <- "prod_data"
export_dir <- "mod_out"
plot_dir <- "figures"


# Data
prod_df <- 
  readRDS(file.path(input_dir,"fsprod_igr_2026-07-13.rds"))
phy_site <- 
  readRDS(file.path(input_dir,"phys_site_predictors_2026-08-12.rds"))
phy_year <- 
  readRDS(file.path(input_dir,"phys_year_predictors_2026-08-12.rds"))
phy_reg_year <- 
  readRDS(file.path(input_dir,"phys_regionyear_predictors_2026-08-12.rds"))

# Prep data  -------------------------------------------------------------------



# Remove NA production estimates
prod_for <- prod_df %>% filter(!is.na(production_mean))

# Create a composite production measure using all species
prod_all <- prod_for %>% 
  group_by(site,cum,date,area,interval) %>% 
  summarise(
    across(contains("sample_den"),sum),
    across(contains("biomass"),sum),
    across(contains("production"),sum)
  ) %>% 
  mutate(species = "all") %>% 
  bind_rows(prod_for)


# Add sample info to production data
samp_df <- phy_site %>% 
  distinct(wateryear,region,site,cum) %>% 
  filter(wateryear != 2024) %>%   ## TEMPORARY ASK NATE FOR NEWEST SHARK RIVER DATA (2025)
  left_join(prod_all, by = join_by(site,cum)) %>% 
  filter(!is.na(production_mean))

# Create a stan input data list for each species
sp <- unique(samp_df$species)
input_list <- lapply(sp, function(s){
  
  # # Data.frame to add year and region to production data
  # samp_df <- phy_site %>% 
  #   distinct(wateryear,region,site,cum)
  # 
  # if(s == "all"){
  #   df <- prod_df %>% 
  #     filter(!is.na(production_med)) %>% 
  #     summarise(
  #       across(contains("biomass"),sum),
  #       across(contains("production"),sum),
  #       .by = c(cum,site)
  #     )
  # } else {
  #   df <- prod_df %>% 
  #     filter(
  #       species == s,
  #       !is.na(production_med)
  #       )
  # }
  
  ## Site level data prep  ##
  # site_df <- df %>% 
  #   
  #   # Create site and year id for random effects
  #   left_join(samp_df,by=join_by(site,cum)) %>% 
  #   filter(wateryear != 2024) %>%   ## TEMPORARY ASK NATE FOR NEWEST SHARK RIVER DATA (2025)
  #   arrange(region,site,cum) %>% 
  #   mutate(site_id = cur_group_id(),.by = c(region,site)) %>% 
  #   mutate(year_id = cur_group_id(),.by = c(wateryear)) %>% 
  #   # mutate(regyear_id = cur_group_id(),.by = c(region,wateryear)) %>% 
  #   mutate(reg_id = cur_group_id(),.by = c(region)) %>% 
  #   
  # # merge in site level predictors
  #   left_join(phy_site,by = join_by(wateryear,region,site,cum))
  
  # Extract data for selected specices
  site_df <- samp_df %>% 
    filter(species == s) %>% 
    
    # Create site and year id for random effects
    arrange(region,site,cum) %>% 
    mutate(site_id = cur_group_id(),.by = c(region,site)) %>% 
    mutate(year_id = cur_group_id(),.by = c(wateryear)) %>% 
    # mutate(regyear_id = cur_group_id(),.by = c(region,wateryear)) %>% 
    mutate(reg_id = cur_group_id(),.by = c(region)) %>% 
    
    # merge in site level predictors
    left_join(phy_site,by = join_by(wateryear,region,site,cum))
  
  ## Region-site bridge
  reg_bridge <- site_df %>% 
    distinct(reg_id,site_id) %>% 
    pull(reg_id)
  
  ## Year-level data prep  ##
  
  # Add year id to year level predictors
  regyear_df <- site_df %>% 
    # distinct(region,wateryear,regyear_id) %>% 
    distinct(region,wateryear,reg_id,year_id) %>% 
    # arrange(regyear_id) %>% 
    arrange(reg_id,year_id) %>% 
    left_join(phy_reg_year,by = join_by(region,wateryear))
  
  
  ## Predictor prep  ##
  
  # Select 1st-level predictors (create column of 1 for intercept)
  x_df <- site_df %>% 
    mutate(int = 1) %>% 
    select(
      int,
      depth,
      dsldd_int,
      plt_cov_int,
      peri_vol_int
    ) %>% 
    
    # scale and center data
    mutate(across(!int,~as.numeric(scale(.x))))
  
  # Select 2nd level predictors (create column of 1 for intercept) and create
  # array for region specific
  # z_df <- regyear_df %>% 
  #   mutate(int = 1) %>% 
  #   select(
  #     int,
  #     wet_sum_365day
  #   ) %>% 
  # 
  # # scale and center data
  # mutate(across(!int,~as.numeric(scale(.x))))
  
  z_list <- lapply(1:max(reg_bridge), function(r){
    z_df <- regyear_df %>% 
      filter(reg_id == r) %>% 
      mutate(int = 1) %>% 
      select(
        int,
        wet_sum_365day
      ) %>% 
        
      # scale and center data
      mutate(across(!int,~as.numeric(scale(.x))))
  })
  z_bind <- abind(z_list,along = 3)
  z_data <- aperm(z_bind,c(3,1,2))
  
  
  ## Stan list  ##
  # stan <-list(
  #   N = nrow(site_df),
  #   L = n_distinct(site_df$year_id),
  #   O = n_distinct(site_df$site_id),
  #   P = n_distinct(site_df$regyear_id),
  #   K = ncol(x_df),
  #   J = ncol(z_df),
  #   y= log(1000*site_df$production_mean+1),
  #   # y=site_df$production_mean*1000,
  #   pp = site_df$regyear_id,
  #   ll = site_df$year_id,
  #   oo = site_df$site_id,
  #   x = x_df,
  #   z = z_df
  # )
  
  stan_data <- list(
    N = nrow(site_df),
    `T` = n_distinct(site_df$year_id),
    S = n_distinct(site_df$site_id),
    R = n_distinct(site_df$reg_id),
    K = ncol(x_df),
    L = ncol(z_data[1,,]),
    y= log(1000*site_df$production_mean+1),
    yr = site_df$year_id,
    st = site_df$site_id,
    rg = reg_bridge,
    x = x_df,
    z = z_data
  )
  
  
  list(stan_data = stan_data, xbridge = site_df, zbridge = regyear_df)
  
}
)
names(input_list) <- sp
input_list_t <- transpose(input_list)
stan_list <- input_list_t$stan_data
xbridge_list <-input_list_t$xbridge
zbridge_list <-input_list_t$zbridge

sapply(stan_list,function(ls) sum(ls$y == 0)/length(ls$y)*100)


# Run models  ------------------------------------------------------------------
## CHANGE TO STEM COUNT FOR PLT NOT COVER
s <-"LUCGOO"

out_list <- lapply(sp[1],function(s) {
  stan_data <- stan_list[[s]]
  
  stan(
    file = file.path(mod_dir,"mvn_second_level_regyear_effects.stan"),
    data = stan_data,
    iter = 5000,
    warmup = 1000,
    chains =4,
    # control = list(adapt_delta = .97),  
    cores = 4
  )
})
names(out_list) <- sp

lapply(1:length(out_list),function(o){
  print(names(out_list)[o])
  check_hmc_diagnostics(out_list[[o]])
  })

lapply(1:length(out_list),function(o){
  print(sp[o])
  print(
    out_list[[o]],
    pars = c(
      "gamma",
      "tau",
      "tau_s",
      # "tau1",
      # "tau2",
      "sigma",
      # "cor_mat",
      # "cor_mat1",
      # "cor_mat_s",
      "mean_gt",
      "sd_gt",
      "min_gt",
      "max_gt",
      "neg_rep")
    )
}
)
out_sum <- summary(out)[[1]]


# Model fit --------------------------------------------------------------------


s <- sp[1]
print(s)
out <- out_list[[s]]
# loo::loo(out)
post <- extract(out,c("y_rep","mu","residuals"))


# Predicitve density
bayesplot::ppc_ecdf_overlay(
  stan_list[[s]]$y,
  post$y_rep[1:100,]
)
a <-bayesplot::ppc_dens_overlay(
  stan_list[[s]]$y,
  post$y_rep[1:100,]
  )
ggsave("total_ppc_dens.png",a,dpi = 300, height = 4,width = 3)
bayesplot::ppc_stat(
  stan_list[[s]]$y,
  post$y_rep,
  stat="min"
)
bayesplot::ppc_stat(
  stan_list[[s]]$y,
  post$y_rep,
  stat="max"
)
bayesplot::ppc_stat(
  stan_list[[s]]$y,
  post$y_rep,
  stat="median"
)
bayesplot::ppc_stat(
  stan_list[[s]]$y,
  post$y_rep,
  stat="sd"
)

# Residuals
resid <- apply(post$residuals,2,mean)
mu_hat <-apply(post$mu,2,mean)


plot(mu_hat,resid)
abline(h=0,col="red")
qqnorm(resid)
qqline(resid)

# variance
plot(mu_hat,resid^2)

# Bayesian Q-Q plot
log_lik <- extract_log_lik(out, parameter_name = "log_lik")
psis <- psis(log_lik)
a <-ppc_loo_pit_qq(
  y = stan_list[[s]]$y,
  yrep = post$y_rep, 
  lw = weights(psis)
  )
ggsave("bayesQ-Q.png",a,dpi = 300, height = 4,width = 3)

# Predrction 
plot(mu_hat,stan_list[[s]]$y)
a <- ggplot(mapping = aes(x = stan_list[[s]]$y,y=mu_hat))+
  geom_point()+
  geom_abline(intercept = 0, slope = 1)+
  theme_classic()+
  xlab("actual")+
  ylab("predicted")
summary(lm(mu_hat~stan_list[[s]]$y))
abline(0,1,col="red")
ggsave("pred_actual.png",a,dpi = 300, height = 4,width = 3)


# Export  ----------------------------------------------------------------------
saveRDS(out_list,file.path(
  export_dir,
  paste0("csi_model_outlist_regyear",Sys.Date(),".rds")
  )
  )
saveRDS(list(stan_list,xbridge_list,zbridge_list),file.path(
  export_dir,
  paste0("bridge_csi_model_outlist_regyear",Sys.Date(),".rds")
)
)


# Plotting  --------------------------------------------------------------------
s <- sp[7]
print(s)
out <- out_list[[s]]
# Prediction input
wet_scaled <- scale(phy_reg_year$wet_sum_365day)
wet <- seq(
  min(wet_scaled),
  max(wet_scaled),
  length.out = 500
)
wet_range <-seq(
  min(phy_reg_year$wet_sum_365day),
  max(phy_reg_year$wet_sum_365day),
  length.out=500
)

# Extract posterior distributions
gamma_out <- extract(out,"gamma")[[1]]
dimnames(gamma_out) <- list(
  NULL,
  # colnames(stan_list[[s]]$z),
  colnames(stan_list[[s]]$z[1,,]),
  colnames(stan_list[[s]]$x)
)

beta_out <-extract(out,"beta")[[1]]
# dimnames(beta_out) <- list(
#   NULL,
#   unique(stan_list[[s]]$pp),
#   colnames(stan_list[[s]]$x)
#   )
dimnames(beta_out) <- list(
  NULL,
  unique(stan_list[[s]]$rg),
  colnames(stan_list[[s]]$x),
  unique(stan_list[[s]]$yr)
  
)


# Summarize posterior distributions
# beta_mean <- apply(beta_out, c(2,3), mean)
# beta_lwr <- apply(beta_out, c(2,3), quantile,.025)
# beta_upr <- apply(beta_out, c(2,3), quantile,.975)
beta_mean <- apply(beta_out, c(2,3,4), mean)
beta_lwr <- apply(beta_out, c(2,3,4), quantile,.025)
beta_upr <- apply(beta_out, c(2,3,4), quantile,.975)

gamma_mean <- apply(gamma_out[,"int",],2,mean)
gamma_lwr <- apply(gamma_out[,"int",],2,quantile,.025)
gamma_upr <- apply(gamma_out[,"int",],2,quantile,.975)


# Create plots for each predictor
preds <- colnames(beta_mean)
lapply(preds, function(p){

  ### Beta plot  ###

  # Create data.frame for each first level predictor
  # beta <- data.frame(
  #   regyear_id = (row.names(beta_mean)),
  #   mean = beta_mean[,p],
  #   lwr = beta_lwr[,p],
  #   upr = beta_upr[,p]
  # ) %>% 
  #   mutate(overlap0 = case_when(
  #     lwr*upr >0 ~ F,
  #     T ~ T
  #   ))
  beta <- lapply(1:dim(beta_mean)[1],function(r){
    data.frame(
      year_id = labels(beta_mean)[[3]],
      mean = beta_mean[r,p,],
      lwr = beta_lwr[r,p,],
      upr = beta_upr[r,p,]
    ) %>% 
      mutate(
        reg_id = r,
        overlap0 = case_when(
          lwr*upr >0 ~ F,
          T ~ T
        )
        )
  }) %>% 
    bind_rows()

  # Add popuation effect
  overall <- data.frame(
    # regyear_id = "overall",
    year_id = "overall",
    reg_id = "overall",
    mean = gamma_mean[p],
    lwr = gamma_lwr[p],
    upr = gamma_upr[p]
  ) %>% 
    mutate(overlap0 = case_when(
      lwr*upr >0 ~ F,
      T ~ T
    ))
  
  # Add in labels for plots
  beta_plot_df <- rbind(overall,beta) %>% 
    left_join(
      zbridge_list[[s]] %>% 
        mutate(year_id = as.character(year_id),reg_id = as.character(reg_id)),
      by = join_by(reg_id,year_id)
      # zbridge_list[[s]] %>% mutate(regyear_id = as.character(regyear_id)),
      # by = "regyear_id"
      ) %>% 
    mutate(wateryear = case_when(
      # regyear_id == "overall" ~ "overall",
      year_id == "overall" ~ "overall",
      T ~ as.character(wateryear)
    ))

  # Create plot
  pd <- position_dodge(width = 0.5)
  beta_plot <- ggplot(beta_plot_df,
                      aes(x = factor(wateryear),
                          y = mean,
                          color = region,
                          alpha = overlap0)) +
    geom_point(position = pd, size = 2) +
    geom_errorbar(aes(ymin = lwr, ymax = upr),
                  position = pd,
                  width = 0.2) +
    geom_hline(yintercept = 0, color = "red") +
    coord_flip() +
    scale_alpha_manual(values = c(`TRUE` = 0.35,
                                  `FALSE` = 1)) +
    theme_classic() +
    theme(
      panel.border = element_rect(color = "black",
                                  fill = NA,
                                  linewidth = 1)
    ) +
    labs(
      x = "Year",
      y = paste0("Slope of annual ", p, "-production relationship"),
      color = "Region",
      alpha = "Overlap"
    )
  print(beta_plot)


  ### Predicted slope plot  ###

  # Predict 1st-level slopes across range of 2nd-level predictors
  pred.all <-NULL
  for(y in 1:1000){
    x <- sample(1:nrow(gamma_out),1)
    pred <- gamma_out[x,"int",p] +
      gamma_out[x,"wet_sum_365day",p]*wet
    pred.all <- rbind(pred.all,pred)
  }
  
  # Summarize predictions
  pred.md <- apply(pred.all, 2, quantile, probs=0.5)
  pred.up <- apply(pred.all, 2, quantile, probs=0.975)
  pred.lo <- apply(pred.all, 2, quantile, probs=0.025)
  pred.df <- data.frame(
    pred.md = pred.md,
    pred.up = pred.up,
    pred.lo = pred.lo,
    wet_range=wet_range
    )
  
  # Add in real data points (reg/year slopes)
  # real_df <- beta %>%
  #   mutate(regyear_id = as.numeric(regyear_id)) %>% 
  #   left_join(zbridge_list[[s]], by ="regyear_id")
  real_df <- beta %>%
    mutate(
      year_id = as.numeric(year_id),
      reg_id = as.numeric(reg_id)
      ) %>% 
    left_join(
      zbridge_list[[s]], 
      by = join_by(reg_id,year_id)
      )
  
  # Create plot
  g_plot <- ggplot(
    data=pred.df,
    aes(x = wet_range,y = pred.md)
    )+
    geom_smooth(color="black")+
    geom_ribbon(
      aes(ymin = pred.lo,ymax =pred.up),
      alpha=0.2,
      fill = "grey"
      )+
    geom_point(
      data=real_df,
      aes(x=wet_sum_365day,y=mean, color = overlap0),
      size = 2.5
      )+
    geom_errorbar(
      data=real_df,
      aes(x=wet_sum_365day,ymin =lwr,ymax = upr,color = overlap0),
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
  print(g_plot)
  
  g_plot_name <- paste0("gamma_plots/gamma_plot_",p,"_",s,".png")
  ggsave(
    file.path(plot_dir,g_plot_name),
    plot = g_plot,
    width = 8,
    height = 6,
    dpi = 300
  )
  
})

s <- sp[1]
print(s)
out <- out_list[[s]]
n.iter = 1000
input.len <- 100
lapply(preds[-1], function(p){
  
  # create input using full range of predictor
  var <- stan_list[[s]]$x[,p]
  input <- seq(min(var),max(var),length.out = input.len)
  lab <- xbridge_list[[s]][,p]
  input_lab <- seq(min(lab),max(lab),length.out = input.len)
  
  # Create a design matrix containing input data for selected variable, and
  # mean values (i.e., 0's) for all but intercept (1)
  X <- matrix(0,nrow=input.len,ncol = ncol(stan_list[[s]]$x))
  colnames(X) <- names(stan_list[[s]]$x)
  X[,"int"] <- 1
  X[,p] <- input
  
  # How many year and regions are there?
  rg_len <- dim(beta_out)[[2]]
  yr_len <- dim(beta_out)[[4]]
  
  # Create year/region combinations
  rg_vec <- rep(seq(rg_len),each = yr_len)
  yr_vec <- rep(seq(yr_len),rg_len)
  
  # Random posterior samples
  iter <-sample(1:dim(beta_out)[1],n.iter)
  
  # Predict production for each region/site along range of predictor.
  ry_list <- map2(rg_vec,yr_vec, function(r,y){
    
    # Create production predictions
    pred_list <- lapply(iter, function(i) X %*% beta_out[i,r,,y])
    # pred_list <- lapply(iter, function(i) exp(X %*% beta_out[i,r,,y])-1)
    pred_bind <-do.call(cbind,pred_list)
    
    # Summarize production predictions
    data.frame(
      reg_id = r,
      year_id =y,
      input = input_lab,
      mean =  apply(pred_bind,1,mean),
      upr = apply(pred_bind,1,quantile,.975),
      lower = apply(pred_bind,1,quantile,.025)
    ) 
  })
  
  # Overall slopes
  overall_list <- lapply(iter, function(i) X %*% gamma_out[i,"int",])
  # overall_list <- lapply(iter, function(i) exp(X %*% gamma_out[i,"int",])+1)
  overall_bind <-do.call(cbind,overall_list)
  overal_df <-data.frame(
    reg_id = NA,
    year_id =NA,
    input = input_lab,
    mean =  apply(overall_bind,1,mean),
    upr = apply(overall_bind,1,quantile,.975),
    lower = apply(overall_bind,1,quantile,.025)
  ) 
  
  # Combine all predictions and prepare for plotting
  pred_df <- bind_rows(ry_list) %>% 
    mutate(regyear_id = cur_group_id(),.by = c(reg_id,year_id)) %>% 
    left_join(zbridge_list[[s]])
  
  # Prepare actual production data for plotting
  real_df <- samp_df %>% 
    filter(species == s) %>% 
    left_join(xbridge_list[[s]]) %>% 
    mutate(production_mean = log(production_mean*1000+1)) %>% 
    # mutate(production_mean = production_mean*1000) %>% 
    left_join(zbridge_list[[s]])
  # real_df$input <- scale(real_df[[p]])[,1]
  real_df$input <- real_df[[p]]
  real_df<- real_df[,c("production_mean","input","wet_sum_365day")] 
  

  b_plot <- ggplot(
    real_df, 
    aes(x = input,y=production_mean)
  )+
    geom_point(colour = "darkgrey",size = .5) +
    geom_line(
      data=overal_df,
      aes(x = input,y=mean),
      inherit.aes = F,
      color = rev(sp_colors)[sp == s],
      linewidth = 1
    )+
    geom_ribbon(
      data=overal_df,
      mapping = aes(x=input,ymin = lower,ymax =upr),
      inherit.aes = F,
      alpha=0.2,
      fill = rev(sp_colors)[sp == s]
    )+
    theme_classic()+
    theme(
      axis.text.x = element_text(size = 18),  
      axis.text.y = element_text(size = 18),
      legend.position = "none",
      panel.border =  element_rect(color = "black", fill = NA, size = 1)
    )+
    ylab("")+
    xlab("")
  print(b_plot)
  b_plot_name <- paste0("beta_plots/beta_plot_",p,"_",s,".png")
  ggsave(
    file.path(plot_dir,b_plot_name),
    plot = b_plot,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  b_plot_all <- ggplot(
    real_df,
    aes(x = input,y=production_mean, colour = wet_sum_365day)
  )+
    geom_point(size = .5) +
    geom_line(
      data = pred_df,
      aes(
        x=input,
        y=mean, 
        group = regyear_id, 
        color = wet_sum_365day),
      inherit.aes = F
    )+
    scale_colour_gradientn(
      colours = c( "gold","forestgreen", "dodgerblue3")
    )+
    theme_classic()+
    theme(
      axis.text.x = element_text(size = 18),  
      axis.text.y = element_text(size = 18),
      legend.position = "none",
      panel.border =  element_rect(color = "black", fill = NA, size = 1)
    )+
    ylab("")+
    xlab("")
   
  print(plot)
  b_plot_all_name <- paste0("beta_plots/beta_plot_all",p,"_",s,".png")
  ggsave(
    file.path(plot_dir,b_plot_all_name),
    plot = b_plot_all,
    width = 8,
    height = 6,
    dpi = 300
  )
  
})


# COefficnet plots  ------------------------------------------------------------

# summarie model coefifcents across all species
coef_df <- lapply(sp, function(s){
  
  # Extract model coefficients
  out <- out_list[[s]]
  gamma_out <- extract(out,"gamma")[[1]]
  dimnames(gamma_out) <- list(
    NULL,
    colnames(stan_list[[s]]$z[1,,]),
    colnames(stan_list[[s]]$x)
  )
  
  # Summarize mean and credible intervals into table
  upr_preds <- dimnames(gamma_out)[[2]]
  lapply(upr_preds,function(u_p){
    data.frame(
      species = s,
      upr_coef = u_p,
      coef = dimnames(gamma_out)[[3]],
      mean = apply(gamma_out[,u_p,],2,mean),
      lwr = apply(gamma_out[,u_p,],2,quantile,.025),
      upr = apply(gamma_out[,u_p,],2,quantile,.975)
    )
  }) %>% bind_rows()%>% 
    mutate(overlap0 = case_when(
      lwr*upr >0 ~ F,
      T ~ T
    )
    )
}
) %>% bind_rows() %>% 
  mutate(
    coef = factor(
      coef,
      levels= rev(c("int","depth","dsldd_int","plt_cov_int","peri_vol_int"))
      ),
    species = forcats::fct_rev(species)
    )
row.names(coef_df) <- NULL


# Plot
sp_colors <- 
  rev(c("black","#b5a331","#339d38","#c26a77","#8c6d3f","#2f2585","#2b695c"))
pd <- position_dodge(width = 0.5)

lapply(upr_preds,function(u_p){
  plot_df <-coef_df %>% filter(upr_coef == u_p)
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
    )
  print(plot)
  plot_name <- paste0("coef_tree_plot_",u_p,".png")
  ggsave(
    file.path(plot_dir,plot_name),
    plot = plot,
    bg = "transparent",
    width = 5,
    height = 10,
    dpi = 300
  )
})




