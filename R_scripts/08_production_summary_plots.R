#-------------------------------------------------------------------------------
#
#  Production summary and plotting 
#
#-------------------------------------------------------------------------------

# AUTHOR: William K. Annis

# CREATED: June 25, 2026

# DESCRIPTION: 


# Housekeeping  ----------------------------------------------------------------
rm(list = ls())

# Load in packages
library(dplyr)
library(ggplot2)
library(tidyr)

# Directories
prod_dir <- "prod_data"
input_dir <- "input_data"
plot_dir <- "figures"

# Data
prod_df <- 
  readRDS(file.path(prod_dir,"fsprod_formatted.rds"))
len_df <-
  readRDS(file.path(input_dir,"fslen_imputed_2026-07-09.rds"))


# Data preparation  ------------------------------------------------------------
# Change zero lengths and P:B into NAs. These were changed to zeros as stan
# hurdle models need zeros not NAs
prod_for <- prod_df %>% 
  mutate(
    across(
      .cols = c(
        mean_wt,
        mean_length,
        interval_mean_wt,
        interval_mean_length,
        ptob),
      .fns = \(x) case_when(
        x == 0 ~ NA,
        x > 0 ~ x
      )
    )
  )

# seasonal response variables
seasonal_prod <- prod_for %>% 
  left_join(len_df %>% distinct(cum,period,waterperiod)) %>% 
  


  
  # Aggergate to sampling interval
  summarise(
    across(
      .cols = c(
        interval_mean_wt,
        interval_density,
        interval_biomass_mean,
        production_mean,
        ptob
        ),
      .fns = mean
    ),
    .by = c(wateryear,waterperiod,cum,species)
  )


year_prod <- prod_for %>% 
  
  # Change daily production to interval production
  mutate(
    across(
      .cols = c(
        production_mean,production_lwr,production_upr
      ),
      .fns = \(x) x*interval
    )
  ) %>% 
  
  # Aggregate production estimate to annual scale
  summarise(
    across(
      .cols = c(
        production_mean,
        interval
      ),
      .fns = sum
    ),
    across(
      .cols = c(
        mean_wt,
        sample_den,
        biomass_mean
      ),
      .fns = mean
    ),
    .by = c(wateryear,region,site,species)
  ) %>% 
  
  # Standardized values to 365 (not all annual intervals are the same)
  mutate(
    across(
      .cols = c(production_mean),
      .fns = \(x) x*(365/interval)
    )#,
    # ptob = production_mean/biomass_mean
  ) %>% 
  
  # Summarized across all sites
  summarize(
    across(
      .cols = c(
        mean_wt,
        sample_den,
        biomass_mean,
        production_mean,
        # ptob
        ),
      .fns = mean, na.rm=T
    ),
    .by = c(wateryear,species)
  ) %>% 
  mutate(ptob = production_mean/biomass_mean)


# General plot details  --------------------------------------------------------
sp <- unique(prod_for$species)
sp <- sp[order(sp)]
sp_colors <- 
  c("black","#b5a331","#339d38","#c26a77","#8c6d3f","#2f2585","#2b695c")

annual_response <- c(
  "mean_wt",
  "sample_den",
  "biomass_mean",
  "production_mean",
  "ptob"
  )
response <- c(
  "interval_mean_wt",
  "interval_density",
  "interval_biomass_mean",
  "production_mean",
  "ptob"
  )

# Response summary  ------------------------------------------------------------


# Total response summary
prod_for %>% 
  filter(
    species == "all",
    ) %>% 
  select(
    interval_mean_wt,
    interval_density,
    interval_biomass_mean,
    production_mean,
    ptob
    ) %>% 
  summary()

# Total PtoB range
prod_for %>% 
  select(ptob) %>% 
  summary()

# Species comparison
prod_for %>% 
  select(
    species,
    interval_mean_wt,
    interval_density,
    interval_biomass_mean,
    production_mean,
    ptob
    ) %>% 
  group_by(species) %>% 
  summarise(across(everything(),~median(.x,na.rm = T))) %>% 
  arrange(interval_biomass_mean)
  
  
# Response boxplots  -----------------------------------------------------------

bar_list <- lapply(response, function(r){
  plot_df <-prod_for 
  plot_df$y <- plot_df[[r]]
  plot <- ggplot(
    data = plot_df,
    aes(
      x = species,
      y= y,
      color = species,
      fill = species
    ))+
    geom_boxplot(fatten = NULL)+
    stat_summary(
      fun = median, 
      geom = "crossbar", 
      fun.min = median, 
      fun.max = median, 
      width = 0.75,       
      color = "white",      
      fatten = 1          
    )+
    theme_classic()+
    theme(
      axis.text.x  = element_blank(),
      axis.text.y = element_text(size = 18),
      legend.position = "none",
      panel.border =  element_rect(
        color = "black", 
        fill = NA, 
        size = 1
        )
    )+
    scale_fill_manual(values = sp_colors)+
    scale_color_manual(values = sp_colors)+
    xlab("")+
    ylab("")
  # print(plot)
  plot
}
)


barplot <- cowplot::plot_grid(
  plotlist = bar_list,
  ncol = 1,
  align = "v"
)

ggsave(
  file.path(
    plot_dir,
    "response_trend",
    "_response_summary.png"
  ),
  plot = barplot,
  width = 8,
  height = 12,
  dpi = 300
)


# Inter-annual response plots---------------------------------------------------

year_list <- lapply(annual_response[annual_response != "ptob"], function(r){
  plot_df <-year_prod
  plot_df$y <- plot_df[[r]]

    theme(legend.position = "none") 
    
  plot<-ggplot(
    data = plot_df%>% filter(species != "all"),
    aes(
      x = wateryear,
      y = y,
      fill = species
    )
  ) +
    geom_col(
      position = "stack",
      width = 0.8
      )+
    scale_fill_manual(
      values = sp_colors[sp!="all"]
    ) +
    geom_line(
      data = plot_df %>% filter(species == "all"),
      aes(
        x = wateryear,
        y = y
      ),
      inherit.aes = F,
      color = "black",
      linewidth = 1
    )+
    theme_classic()+
    theme(
      axis.text.x = element_text(size = 18),  
      axis.text.y = element_text(size = 18),
      legend.position = "none",
      panel.border =  element_rect(color = "black", fill = NA, size = 1)
    )+
    xlab("")+
    ylab("")
  print(plot)
  plot
})
names(year_list) <- response[response != "ptob"]

year_list$ptob <- ggplot(
  data = year_prod,
  aes(
    x = wateryear,
    y=ptob, 
    group = species,
    fill = species
  )
)+
  geom_line(
    aes(
      color =  species,
      linewidth = species
      )
    )+
  # geom_ribbon(
  #   mapping = aes(ymin = production_lwr,ymax =production_upr),
  #   alpha=0.2,
  # )+
  scale_color_manual(values =sp_colors)+
  scale_fill_manual(values =sp_colors)+
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

year_plot <- cowplot::plot_grid(
  plotlist = year_list,
  ncol = 1,
  align = "v"
)

ggsave(
  file.path(
    plot_dir,
    "response_trend",
    "_annual_response.png"
  ),
  plot = year_plot,
  width = 8,
  height = 12,
  dpi = 300
)
  



# Intra-annual response plots  -------------------------------------------------
season_list <- lapply(response,function(r){
  plot_df <-seasonal_prod %>% filter(species == "all")
  plot_df$y <- plot_df[[r]]
  plot <-ggplot(
    data = plot_df ,
    aes(
      x = waterperiod,
      y=y, 
      group = wateryear,
      fill = wateryear
    )
  )+
    geom_line(aes(color = wateryear),linewidth = 1.5)+
    scale_color_gradientn(
      colours = c("#5E3C99", "#CC79A7", "#E66101")
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
  print(plot)
  plot
})
names(season_list) <- response

season_plot <- cowplot::plot_grid(
  plotlist = season_list,
  ncol = 1,
  align = "v"
)

ggsave(
  file.path(
    plot_dir,
    "response_trend",
    "_season_response.png"
  ),
  plot = season_plot,
  width = 8,
  height = 12,
  dpi = 300
)


# Species seasonal plot  -------------------------------------------------------

avg_seasonal_prod <- prod_for %>% 
  left_join(len_df %>% distinct(cum,period,waterperiod)) %>% 
  
  # Aggergate to sampling interval
  summarise(
    across(
      .cols = c(
        # sample_den,
        interval_density,
        # biomass_mean,biomass_lwr,biomass_upr,
        interval_biomass_mean,interval_biomass_lwr,interval_biomass_upr,
        production_mean,production_lwr,production_upr,
        ptob,
        interval_mean_wt),
      .fns = mean
    ),
    .by = c(waterperiod,species)
  )

sp_season_list <- lapply(response[response != "ptob"], function(r){
  plot_df <-avg_seasonal_prod
  plot_df$y <- plot_df[[r]]
  
  theme(legend.position = "none") 
  
  plot<-ggplot(
    data = plot_df%>% filter(species != "all"),
    aes(
      x = waterperiod,
      y = y,
      fill = species
    )
  ) +
    geom_col(
      position = "stack",
      width = 0.8
    )+
    scale_fill_manual(
      values = sp_colors[sp!="all"]
    ) +
    geom_line(
      data = plot_df %>% filter(species == "all"),
      aes(
        x = waterperiod,
        y = y
      ),
      inherit.aes = F,
      color = "black",
      linewidth = 1
    )+
    theme_classic()+
    theme(
      axis.text.x = element_blank(),    
      axis.text.y = element_text(size = 18),
      legend.position = "none",
      panel.border =  element_rect(color = "black", fill = NA, size = 1)
    )+
    xlab("")+
    ylab("")
  print(plot)
  plot
})
names(sp_season_list) <- response[response != "ptob"]

sp_season_list$ptob <- ggplot(
  data = avg_seasonal_prod ,
  aes(
    x = waterperiod,
    y=ptob, 
    group = species,
    fill = species
  )
)+
  geom_line(aes(color = species),linewidth = 1.5)+
  scale_color_manual(values = sp_colors)+
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

sp_season_plot <- cowplot::plot_grid(
  plotlist = sp_season_list,
  ncol = 1,
  align = "v"
)

ggsave(
  file.path(
    plot_dir,
    "response_trend",
    "_sp_season_response.png"
  ),
  plot = sp_season_plot,
  width = 8,
  height = 12,
  dpi = 300
)
