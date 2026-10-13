pred_prop_change <-data_list$predictions_prop_change

library(dplyr)
library(ggplot2)



## ADD PLOTS WITH TOTA PRODUCTION AND SPECIES



a <-pred_prop_change %>% 
  filter(
    response %in% c("interval_mean_wt","interval_density","ptob"),
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


net_data <- a %>% 
  group_by(x_pos) %>% 
  summarise(change_md = sum(change_md))




# Separate data into positive and negative
df_pos <- a %>% filter(change_md >= 0)
df_neg <- a %>% filter(change_md < 0)
df_pos$response <- factor(df_pos$response, levels = rev(c("interval_density","interval_mean_wt","ptob")))
df_neg$response <- factor(df_neg$response, levels = rev(c("interval_density","interval_mean_wt","ptob")))

comp_plot <- ggplot() +
  geom_col(
    data = df_pos, 
    aes(x = x_pos, y = change_md, fill = response), 
    position = "stack" 
    ) +
  geom_col(
    data = df_neg, 
    aes(x = x_pos, y = change_md, fill = response), 
    position = "stack"
    ) +
  scale_fill_brewer(palette = "Set2") + 
  geom_hline(
    yintercept = 0, 
    color = "black", 
    linetype = "dashed"
    ) +
  # geom_errorbar(
  #   data = net_data, 
  #   aes(x = x_pos, ymin = change_md, ymax = change_md), 
  #   color = "black", 
  #   size = .8, 
  #   width = 1.8,
  #   inherit.aes = F
  #   ) +
  geom_point(
    data = net_data, 
    aes(x = x_pos, y = change_md),
    color = "black", 
    size = 3,
    inherit.aes = F
    ) +
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
  
