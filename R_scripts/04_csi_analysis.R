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
library(rstan)

# directories
mod_dir <- "stan_scripts"
input_dir <- "prod_data"
export_dir <- "mod_out"


# Data
prod_df <- readRDS(file.path(input_dir,"fsprod_igr_2026-07-13.rds"))
phy_site <- readRDS(file.path(input_dir,"phys_site_predictors_2026-07-13.rds"))
phy_year <- readRDS(file.path(input_dir,"phys_year_predictors_2026-07-13.rds"))
phy_reg_year <- readRDS(file.path(input_dir,"phys_regionyear_predictors_2026-07-14.rds"))

# Prep data  -------------------------------------------------------------------
sp <- unique(prod_df$species)

input_list <- lapply(sp, function(s){
  
  # Data.frame to add year and region to production data
  samp_df <- phy_site %>% 
    distinct(wateryear,region,site,cum)
  
  ## Site level data prep  ##
  site_df <- prod_df %>% 
    
    # Select sampling events with production estimates for selected species
    filter(
      species == s,
      !is.na(production_med)
      ) %>% 
    
    # Create site and year id for random effects
    left_join(samp_df,by=join_by(site,cum)) %>% 
    filter(wateryear != 2024) %>%  ## TEMPORARY ASK NATE FOR NEWEST SHARK RIVER DATA (2025)
    arrange(region,site,cum) %>% 
    mutate(site_id = cur_group_id(),.by = c(region,site)) %>% 
    mutate(year_id = cur_group_id(),.by = c(wateryear)) %>% 
    
    # merge in site level predictors
    left_join(phy_site,by = join_by(wateryear,region,site,cum)) 
  
  
  ## Year-level data prep  ##
  
  # Add year id to year level predictors
  year_df <- site_df %>% 
    distinct(wateryear,year_id) %>% 
    arrange(year_id) %>% 
    right_join(phy_year,by = join_by(wateryear))
  
  
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
  
  # Select 2nd level predictors (create column of 1 for intercept)
  z_df <- year_df %>% 
    mutate(int = 1) %>% 
    select(
      int,
      wet_sum_365day
    ) %>% 
    
    # scale and center data
    mutate(across(!int,~as.numeric(scale(.x))))
  
  
  ## Stan list  ##
  list(
    N = nrow(site_df),
    L = n_distinct(site_df$year_id),
    O = n_distinct(site_df$site_id),
    K = ncol(x_df),
    J = ncol(z_df),
    y= log(1000*site_df$production_mean+1),
    # y=site_df$production_mean,
    ll = site_df$year_id,
    oo = site_df$site_id,
    x = x_df,
    z1 = z_df
  )
  
}
)
names(input_list) <- sp




# Run models  ------------------------------------------------------------------
## CHANGE TO STEM COUNT FOR PLT NOT COVER
stan_data <- input_list$LUCGOO

out <- stan(
  # file = file.path(mod_dir,"mvn_second_level_effects_site_effect.stan"),
  file = file.path(mod_dir,"mvn_second_level_regyear_effects.stan"),
  data = stan_data,
  iter = 5000,
  warmup = 1000,
  chains =4,
  # control = list(adapt_delta = .97),  
  cores = 4
)
print(out,pars = c("gamma1","tau1","tau2","sigma","cor_mat1"))
print(out,pars = c("gamma","tau","tau1","tau2","sigma","cor_mat1"))
out_sum <- summary(out)[[1]]


# Beta plot  -------------------------------------------------------------------



wet_scaled <- scale(phy_year$wet_sum_365day)
wet <- seq(
  min(wet_scaled),
  max(wet_scaled),
  length.out = 500
)
wet_range <-seq(
  min(phy_year$wet_sum_365day),
  max(phy_year$wet_sum_365day),
  length.out=500
)

gamma_out <- extract(out,"gamma1")[[1]]
dimnames(gamma_out) <- list(
  NULL,
  colnames(input_list[[s]]$z1),
  colnames(input_list[[s]]$x)
)


beta_out <-extract(out,"beta1")[[1]]

dimnames(beta_out) <- list(
  NULL,
  1995:2024,
  colnames(input_list[[s]]$x)
  )

beta_mean <- apply(beta_out, c(2,3), mean)
beta_lwr <- apply(beta_out, c(2,3), quantile,.025)
beta_upr <- apply(beta_out, c(2,3), quantile,.975)

gamma_mean <- apply(gamma_out[,"int",],2,mean)
gamma_lwr <- apply(gamma_out[,"int",],2,quantile,.025)
gamma_upr <- apply(gamma_out[,"int",],2,quantile,.975)

preds <- colnames(beta_mean)



lapply(preds, function(p){
  
  # Create data.frame for each first level predictor
  beta <- data.frame(
    year = (row.names(beta_mean)),
    mean = beta_mean[,p],
    lwr = beta_lwr[,p],
    upr = beta_upr[,p]
  ) %>% 
    mutate(overlap0 = case_when(
      lwr*upr >0 ~ F,
      T ~ T
    ))
  
  overall <- data.frame(
    year = "overall",
    mean = gamma_mean[p],
    lwr = gamma_lwr[p],
    upr = gamma_upr[p]
  ) %>% 
    mutate(overlap0 = case_when(
      lwr*upr >0 ~ F,
      T ~ T
    ))
  
  beta_plot_df <- rbind(overall,beta)


  beta_plot <- ggplot(beta_plot_df,aes(x = year,y = mean,color = overlap0))+
    scale_color_manual(values = c("black","grey"))+
    geom_point()+
    geom_errorbar(aes(ymin =lwr,ymax = upr))+
    geom_abline(slope = 0,intercept = 0, color = "red")+
    coord_flip()+
    theme_classic()+
    theme(
      legend.position = "none",
      panel.border =  element_rect(color = "black", fill = NA, size = 1))+
    ylab(paste0("Slope of annual ",p,"-production relationship"))+
    xlab("Year")
  print(beta_plot)


# Predicted slope plot  --------------------------------------------------------


  pred.all <-NULL
  for(y in 1:1000){
    x <- sample(1:nrow(gamma_out),1)
    pred <- gamma_out[x,"int",p] +
      gamma_out[x,"wet_sum_365day",p]*wet
    pred.all <- rbind(pred.all,pred)
  }
  pred.md <- apply(pred.all, 2, quantile, probs=0.5)
  pred.up <- apply(pred.all, 2, quantile, probs=0.975)
  pred.lo <- apply(pred.all, 2, quantile, probs=0.025)
  
  pred.df <- data.frame(
    pred.md = pred.md,
    pred.up = pred.up,
    pred.lo = pred.lo,
    wet_range=wet_range
    )
  # Second level
  d <- beta %>% left_join(phy_year %>% mutate(year = as.factor(wateryear)))
  
  plot <- ggplot(data=pred.df,aes(x = wet_range,y = pred.md))+
    geom_smooth(color="black")+
    geom_ribbon(aes(ymin = pred.lo,ymax =pred.up),alpha=0.2,fill = "grey")+
    geom_point(data=d,aes(x=wet_sum_365day,y=mean, color = overlap0))+
    geom_errorbar(data=d,aes(x=wet_sum_365day,ymin =lwr,ymax = upr,color = overlap0,),inherit.aes = FALSE)+
    scale_color_manual(values = c("black","grey"))+
    geom_abline(slope = 0,intercept = 0, color = "red")+
    theme_classic()+
    theme(legend.position = "none",panel.border =  element_rect(color = "black", fill = NA, size = 1))+
    ylab(paste0("Slope of annual ",p,"-production relationship"))+
    # xlab("Depth - 365 day average")
    xlab("Wetted days")
  print(plot)
  
})





