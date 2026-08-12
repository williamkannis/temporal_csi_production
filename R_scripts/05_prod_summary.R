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
  readRDS(file.path(prod_dir,"fsprod_igr_2026-07-13.rds"))
phy_site <- 
  readRDS(file.path(prod_dir,"phys_site_predictors_2026-07-13.rds"))
phy_year <- 
  readRDS(file.path(prod_dir,"phys_year_predictors_2026-07-13.rds"))
phy_reg_year <- 
  readRDS(file.path(prod_dir,"phys_regionyear_predictors_2026-07-14.rds"))
len_df <- 
  readRDS(file.path(input_dir,"fslen_imputed_2026-07-09.rds"))

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
  mutate(
    species = "all",
    PtoB_mean = production_mean/biomass_mean
    ) %>% 
  bind_rows(prod_for)


# Add sample info to production data
samp_df <- phy_site %>% 
  distinct(wateryear,region,site,cum) %>% 
  filter(wateryear != 2024) %>%   ## TEMPORARY ASK NATE FOR NEWEST SHARK RIVER DATA (2025)
  left_join(prod_all, by = join_by(site,cum)) %>% 
  filter(!is.na(production_mean)) %>% 
  
  # add watr period
  left_join(
    len_df %>% distinct(cum,waterperiod),
    by = join_by(cum)
    )


# Total summary boxplots  ------------------------------------------------------

tot_bio <- ggplot(
  data = samp_df %>% filter(species == "all"),
  aes(
    x = species,
    y= biomass_mean,
  ))+
  geom_boxplot(fatten = NULL,fill="black")+
  stat_summary(
    fun = median, 
    geom = "crossbar", 
    fun.min = median, 
    fun.max = median, 
    width = 0.75,       
    color = "white",      
    fatten = 2          
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
tot_bio_name <- paste0("prod_plots/tot_bio.png")
ggsave(
  file.path(plot_dir,tot_bio_name),
  plot = tot_bio,
  width = 4.6,
  height = 8,
  dpi = 300
)

tot_prod <- ggplot(
  data = samp_df %>% filter(species == "all"),
  aes(
    x = species,
    y= production_mean,
  ))+
  geom_boxplot(fatten = NULL,fill="black")+
  stat_summary(
    fun = median, 
    geom = "crossbar", 
    fun.min = median, 
    fun.max = median, 
    width = 0.75,       
    color = "white",      
    fatten = 2          
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
tot_prod_name <- paste0("prod_plots/tot_prod.png")
ggsave(
  file.path(plot_dir,tot_prod_name),
  plot = tot_prod,
  width = 4.6,
  height = 8,
  dpi = 300
)

tot_bio <- ggplot(
  data = samp_df %>% filter(species == "all"),
  aes(
    x = species,
    y= biomass_mean,
  ))+
  geom_boxplot(fatten = NULL,fill="black")+
  stat_summary(
    fun = median, 
    geom = "crossbar", 
    fun.min = median, 
    fun.max = median, 
    width = 0.75,       
    color = "white",      
    fatten = 2          
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
tot_bio_name <- paste0("prod_plots/tot_bio.png")
ggsave(
  file.path(plot_dir,tot_bio_name),
  plot = tot_bio,
  width = 4.6,
  height = 8,
  dpi = 300
)

tot_p2b <- ggplot(
  data = samp_df %>% filter(species == "all"),
  aes(
    x = species,
    y= PtoB_mean,
  ))+
  geom_boxplot(fatten = NULL,fill="black")+
  stat_summary(
    fun = median, 
    geom = "crossbar", 
    fun.min = median, 
    fun.max = median, 
    width = 0.75,       
    color = "white",      
    fatten = 2          
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
tot_p2b_name <- paste0("prod_plots/tot_p2b.png")
ggsave(
  file.path(plot_dir,tot_p2b_name),
  plot = tot_p2b,
  width = 4.6,
  height = 8,
  dpi = 300
)




# Species composition ----------------------------------------------------------
sp_df <-samp_df %>% 
  select(wateryear,cum,region,site,species,production_mean) %>% 
  pivot_wider(
    names_from = species,
    values_from = production_mean
    ) %>% 
  mutate(across(
    .cols = c(FUNCHR,GAMHOL,HETFOR,JORFLO,LUCGOO,POELAT),
    .fns = \(x) (x/all)*100
  )) %>% 
  
  # Sites with no production cant have species compostion
  filter(all > 0) %>% 
  select(-all) %>% 
  pivot_longer(
    cols=c(FUNCHR,GAMHOL,HETFOR,JORFLO,LUCGOO,POELAT),
    names_to = "species",
    values_to = "prop_prod"
    ) %>% 
  mutate(
    species = factor(
      species,
      c("LUCGOO","GAMHOL","FUNCHR","JORFLO","HETFOR","POELAT")
      )
    )

# number o samples with production values
sp_df %>% 
  distinct(cum,site) %>% 
  nrow()
# 2573

box_col <- c("#2f2585","#339d38","#b5a331","#8c6d3f","#c26a77","#2b695c")
sp_plot <- ggplot(
  data = sp_df,
  aes(
    x = species,
    y= prop_prod,
    fill = species,
    color = species
  ))+
  geom_boxplot(fatten = NULL)+
  stat_summary(
    fun = median, 
    geom = "crossbar", 
    fun.min = median, 
    fun.max = median, 
    width = 0.75,       
    color = "white",      
    fatten = 2          
  )+
  scale_color_manual(values =box_col)+
  scale_fill_manual(values =box_col)+
  theme_classic()+
  theme(
    axis.text.x = element_text(size = 18),  
    axis.text.y = element_text(size = 18),
    legend.position = "none",
    panel.border =  element_rect(color = "black", fill = NA, size = 1)
  )+
  xlab("")+
  ylab("")
sp_plot_name <- paste0("prod_plots/sp_comp.png")
ggsave(
  file.path(plot_dir,sp_plot_name),
  plot = sp_plot,
  width = 18,
  height = 8,
  dpi = 300
)






# Species comparison  ----------------------------------------------------------
sp_prod <- samp_df %>% 
  
  # Aggergate to sampling interval
  summarise(
    across(
      .cols = c(production_mean,production_lwr,production_upr),
      .fns = mean
    ),
    .by = c(wateryear,waterperiod,cum,species)
  ) %>% 
  mutate(
    across(
      .cols = c(production_mean,production_lwr,production_upr),
      .fns = \(x) x*1000
    ))


ggplot(
  data = sp_prod %>% filter(species == "all"),
  aes(
    x = cum,
    y=production_mean 
    # group = region,
    # fill = region
  )
)+
  geom_line()+
  geom_ribbon(
    mapping = aes(ymin = production_lwr,ymax =production_upr),
    alpha=0.2,
  )+
  # scale_color_manual(values =reg_colors)+
  # scale_fill_manual(values =reg_colors)+
  theme_classic()+
  xlab("")+
  ylab("")


# Seasonal changes  
season_plot <-ggplot(
  data = sp_prod %>% filter(species == "all"),
  aes(
    x = waterperiod,
    y=production_mean, 
    group = wateryear,
    fill = wateryear
  )
)+
  geom_line(aes(color = wateryear),linewidth = 1.5)+
  # geom_ribbon(
  #   mapping = aes(ymin = production_lwr,ymax =production_upr),
  #   alpha=0.2,
  # )+
  scale_color_gradientn(
    colours = c("#5E3C99", "#CC79A7", "#E66101")
    )+
  theme_classic()+
  theme(
    axis.text.x = element_text(size = 18),  
    axis.text.y = element_text(size = 18),
    legend.position = "none",
    panel.border =  element_rect(color = "black", fill = NA, size = 1),
    panel.background = element_rect(fill = "transparent", colour = NA),  # Transparent panel
    plot.background = element_rect(fill = "transparent", colour = NA)
  )+

  xlab("")+
  ylab("")
s_plot_name <- paste0("prod_plots/total_seasonal.png")
ggsave(
  file.path(plot_dir,s_plot_name),
  plot = season_plot,
  bg = "transparent",
  width = 14,
  height = 8,
  dpi = 300
)

season_phy <- 
  phy_site %>% 
  left_join(samp_df %>% distinct(cum,waterperiod)) %>% 
  summarise(
    depth = mean(depth,na.rm=T),
    plt_cov_int = mean(plt_cov_int,na.rm=T),
    .by = c(wateryear,waterperiod)
    ) %>% 
  left_join(phy_year)

ggplot(
  data = season_phy,
  aes(
    x = waterperiod,
    y=plt_cov_int, 
    group = wateryear,
    fill = wateryear
  )
)+
  geom_line(aes(color = wateryear),linewidth = 1.5)+
  # geom_ribbon(
  #   mapping = aes(ymin = production_lwr,ymax =production_upr),
  #   alpha=0.2,
  # )+
  scale_color_gradientn(
    colours = c("#5E3C99", "#CC79A7", "#E66101")
  )+
  theme_classic()+
  theme(
    axis.text.x = element_text(size = 18),  
    axis.text.y = element_text(size = 18),
    # legend.position = "none",
    panel.border =  element_rect(color = "black", fill = NA, size = 1),
    panel.background = element_rect(fill = "transparent", colour = NA),  # Transparent panel
    plot.background = element_rect(fill = "transparent", colour = NA)
  )+
  
  xlab("")+
  ylab("")

d_plot <-ggplot(phy_year,aes(x=wateryear,y=depth_ave_365day)) + 
  geom_line(colour="blue",linewidth = 2) +
  theme_classic()+
  theme(
    axis.text.x = element_text(size = 30),  
    axis.text.y = element_text(size = 30),
    legend.position = "none",
    panel.border =  element_rect(color = "black", fill = NA, size = 1),
    panel.background = element_rect(fill = "transparent", colour = NA),  # Transparent panel
    plot.background = element_rect(fill = "transparent", colour = NA)
  )+
  
  xlab("")+
  ylab("")
d_plot_name <- paste0("prod_plots/annual_depth.png")
ggsave(
  file.path(plot_dir,d_plot_name),
  plot = d_plot,
  bg = "transparent",
  width = 26,
  height = 8,
  dpi = 300
)


# Total region:year  -----------------------------------------------------------

year_prod <- samp_df %>% 
  # filter(species == "all") %>% 
  
  # Change daily production to interval production
  mutate(
    across(
      .cols = c(production_mean,production_lwr,production_upr),
      .fns = \(x) x*interval
      )
    ) %>% 
  
  # Aggreage production estiate to annual scale
  summarise(
    across(
      .cols = c(production_mean,production_lwr,production_upr,interval),
      .fns = sum
      ),
    .by = c(wateryear,region,site,species)
  ) %>% 
  
  # Standardized values to 365 (not all annual intervals are the same)
  mutate(
    across(
      .cols = c(production_mean,production_lwr,production_upr),
      .fns = \(x) x*(365/interval)
    )
  ) 

# region:year level
regyear_prod <-year_prod %>% 
  summarise(
    across(
      .cols = c(production_mean,production_lwr,production_upr),
      .fns = mean
    ),
    .by = c(wateryear,region,species)
  )

# Plot
reg_colors <- c("forestgreen","gold", "dodgerblue3")
ry_plot <- ggplot(
  data = regyear_prod %>% filter(species == "all"),
  aes(
    x = wateryear,
    y=production_mean, 
    group = region,
    fill = region
    )
  )+
  geom_line(aes(color = region))+
  geom_ribbon(
    mapping = aes(ymin = production_lwr,ymax =production_upr),
    alpha=0.2,
  )+
  scale_color_manual(values =reg_colors)+
  scale_fill_manual(values =reg_colors)+
  theme_classic()+
  theme(
    axis.text.x = element_text(size = 18),  
    axis.text.y = element_text(size = 18),
    legend.position = "none",
    panel.border =  element_rect(color = "black", fill = NA, size = 1)
  )+
  xlab("")+
  ylab("")
b_plot_all_name <- paste0("prod_plots/total_reg_year.png")
ggsave(
  file.path(plot_dir,b_plot_all_name),
  plot = ry_plot,
  width = 14,
  height = 8,
  dpi = 300
)


ggplot(
  data = phy_reg_year %>% filter(wateryear != 2024),
  aes(
    x = wateryear,
    y=depth_ave_365day, 
    group = region,
    color = region
    )
  )+
  geom_line() 

# Specie year level
sp_year <- year_prod %>% 
  summarise(
    across(
      .cols = c(production_mean,production_lwr,production_upr),
      .fns = mean
    ),
    .by = c(wateryear,species)
  )
sp_colors <- 
  c("#b5a331","#339d38","#c26a77","#8c6d3f","#2f2585","#2b695c")
ggplot(
  data = sp_year %>% filter(species != "all"),
  aes(
    x = wateryear,
    y=production_mean, 
    group = species,
    fill = species
  )
)+
  geom_line(aes(color = species))+
  geom_ribbon(
    mapping = aes(ymin = production_lwr,ymax =production_upr),
    alpha=0.2,
  )+
  scale_color_manual(values =sp_colors)+
  scale_fill_manual(values =sp_colors)+
  theme_classic()

# Phy data ----
a <-phy_site %>% 
  left_join(len_df %>% distinct(cum,waterperiod)) %>% 
  summarise(
    depth = mean(plt_cov_int,na.rm=T),
    .by = c(wateryear,waterperiod,cum)
  ) %>% 
  left_join(phy_year)

ggplot(
  data = a,
  aes(
    x = waterperiod,
    y=depth, 
    group = wet_sum_365day,
    fill = wet_sum_365day
  )
)+
  geom_line(aes(color = wet_sum_365day),linewidth = 1.5)


a <-phy_site %>% 
  left_join(len_df %>% distinct(cum,waterperiod)) %>% 
  # summarise(
  #   depth = mean(plt_cov_int,na.rm=T),
  #   .by = c(wateryear,region)
  # ) %>% 
  left_join(phy_year) %>% 
  mutate(id = cur_group_id(),.by = c(region,wateryear))

plot(a$wet_sum_365day,a$plt_cov_int)

dep_plot <- ggplot(
  data = a,
  aes(
    x = wet_sum_365day,
    y = depth,
    group = wet_sum_365day
    # fill = wet_sum_365day,
    # color = wet_sum_365day
    )
  )+
  geom_boxplot(
    width = 0.9,
    fatten = NULL,
    color = "darkblue",
    fill = "darkblue")+
  stat_summary(
    fun = median, 
    geom = "crossbar", 
    fun.min = median, 
    fun.max = median, 
    width = 0.75,       
    color = "white",      
    fatten = 2          
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
dep_plot_name <- paste0("prod_plots/dep_bar.png")
ggsave(
  file.path(plot_dir,dep_plot_name),
  plot = dep_plot,
  width = 18,
  height = 5,
  dpi = 300
)

veg_plot <-ggplot(data = a,aes(x = wet_sum_365day,y = plt_cov_int,group = wet_sum_365day))+
  geom_boxplot(
    width = 0.9,
    fatten = NULL,
    color = "darkgreen",
    fill = "darkgreen")+
  stat_summary(
    fun = median, 
    geom = "crossbar", 
    fun.min = median, 
    fun.max = median, 
    width = 0.75,       
    color = "white",      
    fatten = 2          
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
veg_plot_name <- paste0("prod_plots/veg_bar.png")
ggsave(
  file.path(plot_dir,veg_plot_name),
  plot = veg_plot,
  width = 18,
  height = 5,
  dpi = 300
)
