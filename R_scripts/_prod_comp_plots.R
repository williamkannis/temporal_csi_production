

plot_dir <- "figures"

pred_prop_change <-data_list$predictions_prop_change

library(dplyr)
library(ggplot2)



## ADD PLOTS WITH TOTA PRODUCTION AND SPECIES


# Prepare prop change data for plots
prop_plot_df <-pred_prop_change %>% 
  
  # Filter to only include credible CSIs, or global values in absence of CSI
  filter(
    x_var %in% c("peri_vol_int","dsldd_int") & 
      z_var == "wet_sum_365day" &
      z_cat == "mean" |
      
      x_var == "depth"  & 
      z_var == "wet_sum_365day" &
      z_cat != "mean" |
      
      x_var == "plt_cov_int"  & 
      z_var == "pisc_index" &
      z_cat != "mean" |
      
      x_var %in% c("wet_sum_365day","pisc_index")
    ) %>% 
  
  # Create factor to allow custom plot spacing
  mutate(
    cat = paste(x_var,z_cat,sep = "_"),
    x_pos = case_when(
      cat == "depth_min"         ~ 1,
      cat == "depth_max"         ~ 3,
      cat == "plt_cov_int_min"   ~ 6,
      cat == "plt_cov_int_max"   ~ 8,
      cat == "dsldd_int_mean"    ~ 11,
      cat == "peri_vol_int_mean" ~ 14,
      cat == "wet_sum_365day_mean" ~ 17,
      cat == "pisc_index_mean"   ~ 20
    ),
    # x_pos = x_pos *-1
    )


# Separate data into positive and negative
bar_data <-  prop_plot_df %>% 
  filter(response %in% c("interval_mean_wt","interval_density","ptob"))
bar_pos <- bar_data %>% filter(prop_change_md >= 0)
bar_neg <- bar_data %>% filter(prop_change_md < 0)
bar_pos$response <- factor(
  bar_pos$response, 
  levels = rev(c("interval_density","interval_mean_wt","ptob"))
  )
bar_neg$response <- factor(
  bar_neg$response, 
  levels = rev(c("interval_density","interval_mean_wt","ptob"))
  )

# Prepare production data
prod_data <- prop_plot_df %>% 
  filter(response == "production_mean")
net_data <- a %>% 
  group_by(x_pos) %>% 
  summarise(
    prop_change_md = sum(prop_change_md),
    prop_change_up = sum(prop_change_up),
    prop_change_lo = sum(prop_change_lo)
  )

# Create plot
comp_plot <- ggplot() +
  
  # Postive stacked bars
  geom_col(
    data = bar_pos, 
    aes(
      x = x_pos, 
      y = prop_change_md, 
      fill = response,
      alpha = overlap0
      ), 
    position = "stack" 
    ) +
  
  # Negative stacked bars
  geom_col(
    data = bar_neg, 
    aes(
      x = x_pos, 
      y = prop_change_md, 
      fill = response,
      alpha = overlap0
      ), 
    position = "stack"
    ) +
  
  # Production mean
  geom_errorbar(
    data = prod_data,
    aes(
      x = x_pos, 
      ymin = prop_change_lo, 
      ymax = prop_change_up,
      alpha = overlap0
    ),
    color = "black",
    size = .8,
    width = 0,
    inherit.aes = F
  ) +
  
  # Production CI
  geom_point(
    data = prod_data, 
    aes(
      x = x_pos, 
      y = prop_change_md,
      alpha = overlap0
    ),
    color = "black", 
    size = 3,
    inherit.aes = F
  ) +
  
  # Zero line
  geom_hline(
    yintercept = 0, 
    color = "black", 
    linetype = "dashed"
    ) +
  
  # Plot aestetics
  scale_fill_brewer(palette = "Set2") + 
  scale_alpha_manual(values = c("TRUE" = 0.3, "FALSE" = 1.0))+
  ylab("Proportion change")+
  theme_classic()+
  theme(
    legend.position="none",
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    panel.border =  element_rect(color = "black", fill = NA, size = 1),
  )

# Export plot
ggsave(
  file.path(
    plot_dir,
    "prod_components",
    "_prod_components_wide.png"
  ),
  comp_plot,
  dpi = 300,
  height = 8,
  width = 10
)  
  
