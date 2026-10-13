#-------------------------------------------------------------------------------
#
#  Driver summary and plotting 
#
#-------------------------------------------------------------------------------

# AUTHOR: William K. Annis

# CREATED: October 8, 2026

# DESCRIPTION: 

## TODO
# 1. color by hydroperiod and predator abundance (maybe not now)

# House keeping  ---------------------------------------------------------------
rm(list=ls())

# Load in R packages
library(ggplot2)
library(dplyr)

# Directories
input_dir <- "input_data"
prod_dir <- "prod_data"
plot_dir <- "figures"

# Load in data
phy_site <- 
  readRDS(file.path(prod_dir,"phys_site_predictors_2026-08-17.rds"))
phy_reg_year <- 
  readRDS(file.path(prod_dir,"phys_regionyear_predictors_2026-08-17.rds"))
prod_df <- 
  readRDS(file.path(prod_dir,"fsprod_formatted.rds"))
len_df <-
  readRDS(file.path(input_dir,"fslen_imputed_2026-07-09.rds"))


# Data preparation  ------------------------------------------------------------

# seasonal drivers - annual variation
seasonal_drivers_year <- phy_site %>% 
  left_join(len_df %>% distinct(cum,period,waterperiod)) %>% 
  right_join(prod_df %>% distinct(cum,site)) %>% 
  left_join(phy_reg_year) %>% 

  # Aggregate to sampling interval
  summarise(
    across(
      .cols = c(
        depth,
        dsldd_int,
        plt_cov_int,
        peri_vol_int,
        wet_sum_365day,
        pisc_index
      ),
      .fns = \(x) mean(x, na.rm = TRUE)
    ),
    .by = c(wateryear,waterperiod,cum)
  )


# seasonal drivers - site variation
seasonal_drivers_site <- phy_site %>% 
  left_join(len_df %>% distinct(cum,period,waterperiod)) %>% 
  right_join(prod_df %>% distinct(cum,site)) %>% 
  
  # Aggregate to sampling interval
  summarise(
    across(
      .cols = c(
        depth,
        dsldd_int,
        plt_cov_int,
        peri_vol_int
      ),
      .fns = \(x) mean(x, na.rm = TRUE)
    ),
    .by = c(region,site,waterperiod)
  )

# Annual drivers
annual_drivers <-  phy_reg_year %>% 
  right_join(prod_df %>% distinct(region,wateryear))
  

# Driver seasonal plot ---------------------------------------------------------

drivers <- c(
  "depth",
  "dsldd_int",
  "plt_cov_int",
  "peri_vol_int"
)

# annual variation
season_year_list <- lapply(drivers,function(r){
  
  
  # Max plotting value
  max_y <-max(seasonal_drivers_year[[r]], seasonal_drivers_site[[r]])
  max_y <-ceiling(max_y / 5) * 5
  
  # Filter plot data
  plot_df <-seasonal_drivers_year
  plot_df$y <- plot_df[[r]]
  plot <-ggplot(
    data = plot_df ,
    aes(
      x = waterperiod,
      y=y, 
      group = wateryear,
      # fill = wateryear
      fill = wet_sum_365day
    )
  )+
    ylim(0,max_y)+
    geom_line(aes(color = wateryear),linewidth = 1.5)+
    # scale_color_gradientn(
    #   colours = c("#5E3C99", "#CC79A7", "#E66101")
    # )+
    scale_color_gradientn(
      colours = c("#FFD700","#1ABC9C","#1B4F72")
    )+
    theme_classic()+
    theme(
      axis.text.x = element_blank(),  
      axis.text.y = element_text(size = 18),
      legend.position = "none",
      panel.border =  element_rect(color = "black", fill = NA, size = 1),
      panel.background = element_rect(fill = "transparent", colour = NA), 
      plot.background = element_rect(fill = "transparent", colour = NA)
    )+
    
    xlab("")+
    ylab("")
  # print(plot)
  plot
})
names(season_year_list) <- drivers

# site variation
season_site_list <- lapply(drivers,function(r){
  
  # Max plotting value
  max_y <-max(seasonal_drivers_year[[r]], seasonal_drivers_site[[r]])
  max_y <-ceiling(max_y / 5) * 5
  
  # Filter plot data
  plot_df <-seasonal_drivers_site
  plot_df$y <- plot_df[[r]]
  plot <-ggplot(
    data = plot_df ,
    aes(
      x = waterperiod,
      y=y, 
      group = site,
      color= region
    )
  )+
    ylim(0,max_y)+
    geom_line(linewidth = 1.5)+
    theme_classic()+
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 18),
      legend.position = "none",
      panel.border =  element_rect(color = "black", fill = NA, size = 1),
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background = element_rect(fill = "transparent", colour = NA)
    )+
    scale_color_manual(values = c(
      "TSL" = "#FFD700",
      "SRS" = "#1ABC9C",
      "WCA" = "#1B4F72")
      )+
    xlab("")+
    ylab("")
  # print(plot)
  plot
})
names(season_site_list) <- drivers

# Combined year ans site
season_list <- c(season_site_list,season_year_list)

season_plot <- cowplot::plot_grid(
  plotlist = season_list,
  ncol = 2,
  align = "vh",
  byrow = F
)

ggsave(
  file.path(
    plot_dir,
    "driver_trend",
    "_season_driver.png"
  ),
  plot = season_plot,
  width = 16,
  height = 12,
  dpi = 300
)


# Inter-annual response plots---------------------------------------------------
annual_driver <- c("wet_sum_365day","pisc_index")
year_list <- lapply(annual_driver, function(r){
  plot_df <-annual_drivers
  plot_df$y <- plot_df[[r]]
  
  theme(legend.position = "none") 
  
  plot<-ggplot(
    data = plot_df,
    aes(
      x = wateryear,
      y=y, 
      group = region,
      fill = region
    )
  )+
    geom_line(
      aes(
        color =  region,
        linewidth = region
      )
    )+
    scale_color_manual(values = c(
      "TSL" = "#FFD700",
      "SRS" = "#1ABC9C",
      "WCA" = "#1B4F72")
      )+
    scale_fill_manual(values = c(
      "TSL" = "#FFD700",
      "SRS" = "#1ABC9C",
      "WCA" = "#1B4F72")
    )+
    scale_linewidth_manual(values = c(2,rep(1,6)))+
    theme_classic()+
    theme(
      axis.text.x = element_text(size = 18),  
      axis.text.y = element_text(size = 18),
      legend.position = "none",
      panel.border =  element_rect(color = "black", fill = NA, size = 1)
    )+
    xlab("")+
    ylab("")
  # print(plot)
  plot
})
names(year_list) <- annual_driver


year_plot <- cowplot::plot_grid(
  plotlist = year_list,
  ncol = 1,
  align = "v"
)

ggsave(
  file.path(
    plot_dir,
    "driver_trend",
    "_annual_driver.png"
  ),
  plot = year_plot,
  width = 8,
  height = 6,
  dpi = 300
)
