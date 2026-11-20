
# House keeping  ---------------------------------------------------------------

# Packages
library(dplyr)
library(ggplot2)

# Directories
plot_dir <- "figures"

# Data
coef_df <- data_list$coef_gamma
theta <- data_list$coef_theta
predictions_theta <- data_list$predictions_theta
pred_response <- data_list$predictions_response
pred_prop_change <-data_list$predictions_prop_change


# CSI Plots  -------------------------------------------------------------------

coef_df %>% 
  group_by(response,species) %>% 
  summarise(n=n()) %>% 
  arrange(species)

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
minmax_slopes <-predictions_theta %>% 
  group_by(
    response,
    species,
    x_var,
    z_var
  ) %>% 
  mutate(
    overlap0 = case_when(
      pred_up*pred_lo > 0 ~ F,
      T ~T
    )
  ) %>% 
  ungroup() %>% 
  filter(
    !is.na(z_cat),
    x_var != "int"
  ) %>% 
  select(-pred_up,-pred_lo,-z_raw_range,-z_range) %>% 
  tidyr::pivot_wider(
    names_from = z_cat,
    values_from = c(pred_md,overlap0)
  ) %>% 
  # mutate(
  #   overlap0_all =case_when(
  #     overlap0_min & overlap0_max ~ T,
  #     T~F
  #   )) %>% 
  
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
# responses <- c("sample_den","biomass_mean","production_mean","ptob")
responses <- c("interval_density","interval_biomass_mean","production_mean","ptob")

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
  
  # # Export plot
  # plot_name <- paste0("coef_csi_plot_",r,"_v4.png")
  # plot_height <- 13*1.5
  # ggsave(
  #   file.path(
  #     plot_dir,
  #     "csi_coef",
  #     plot_name
  #   ),
  #   plot = plot,
  #   bg = "transparent",
  #   width = 5,
  #   height = plot_height,
  #   dpi = 300
  # )
})

# Gamma plots  ----------------------------------------------------------------- 

theta_cb <- predictions_theta %>% 
  distinct(species,response,z_var,x_var)

theta_plots <- lapply(1:nrow(gamma_cb), function(i){
  
  # Extract predictions for selected variable combination
  selection <- theta_cb[i,]
  pred_df <- predictions_theta %>% 
    semi_join(
      selection,
      by = join_by(response, species, x_var, z_var)
      )
  
  # Extract actual theta
  act_df <- theta %>% 
    semi_join(
      selection,
      by = join_by(response, species, x_var)
    )
  
  # Assign z data based on predictor
  act_df$z <- act_df[[selection$z_var]]
    
  # x axis name
  if(predictor == "wet_sum_365day") x_name <- "Hydroperiods (days)"
  if(predictor == "pisc_index") x_name <- "Predator abundance (CPUE)"
  
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
      data=act_df,
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
    ylab("")+
    xlab(x_name)+
    theme_classic()+
    theme(
      axis.title.x = element_text(size = 18),
      axis.text.x = element_text(size = 18),  
      axis.text.y = element_text(size = 18),
      legend.position = "none",
      panel.border =  element_rect(color = "black", fill = NA, size = 1),
      plot.margin = margin(r = 20, unit = "pt")
    )
  # print(g_plot)
  g_plot
  
  
})
names(theta_plots) <- paste(
  theta_combo$species,
  theta_combo$response,
  theta_combo$z_var,
  var_combo$x_var,
  sep = "_"
  )
plot_list
  
plotter <- function(species,response,z,x) {
  outer_name <- paste(species,response,sep="_")
  inner_name <- paste(z,x,sep="_")
  gamma_plots[[outer_name]][[inner_name]]
}
plotter("POELAT","ptob","pisc_index","peri_vol_int")

p <-list(
  plotter("all","interval_density","wet_sum_365day","depth"),
  plotter("all","interval_biomass_mean","wet_sum_365day","depth"),
  plotter("all","production_mean","wet_sum_365day","depth"),
  NULL,
  plotter("all","interval_density","pisc_index","plt_cov_int"),
  plotter("all","interval_biomass_mean","pisc_index","plt_cov_int"),
  plotter("all","production_mean","pisc_index","plt_cov_int"),
  plotter("all","ptob","wet_sum_365day","plt_cov_int")
  
)

a <-cowplot::plot_grid(
  plotlist = p,
  ncol = 4,
  nrow = 2,
  align = "hv",
  axis = "none",
  byrow = TRUE
)

ggsave(
  file.path(
    plot_dir,
    "csi_plot",
    "csi_plot_all.png"
  ),
  dpi = 600,
  height = 8,
  width = 16
)

# First level plot with CSI  ---------------------------------------------------

## THINK IF YOU WANT ON LOG SCALE

# Extract drivers with at least one credible CSI
# Remove non-significance CSIs
csi_pred_df <- pred_response %>% 
  semi_join(
    sig_df,
    by = join_by(species, response, x_var, z_var)
  ) %>% 
  filter(z_cat != "mean")


# Extract drivers without credible CSIs
mean_pred_df <- pred_response %>% 
  
  # Change one set of mean parameters (doesn't matter which, both are equal) 
  # to intercept. These will be used in absence of sig CSI
  filter(
    z_var == "pisc_index",
    z_cat == "mean"
  ) %>% 
  mutate(z_cat = "int") %>% 
  
  # Remove signficant CSI
  anti_join(
    sig_df,
    by = join_by(species, response, x_var)
  )

# Cmobine genral and CSI effects
pred_df <- mean_pred_df %>% 
  bind_rows(csi_pred_df)

# All combination of spcies, plots, and drivers to plot
response_cb <- pred_response %>% 
  group_by(species, response, x_var) %>% 
  summarise(
    y_max = max(pred_md, pred_lo, pred_up)*1.01,
    y_min <- min(for_df$pred_md, for_df$pred_lo, for_df$pred_up)*1.01
  )
  distinct(species, response, x_var)

# Set min and max y value for each species and resonse
response_minmax <- pred_response %>% 
  group_by(species, response) %>% 
  summarise(
    y_max = max(pred_md, pred_lo, pred_up)*1.01,
    y_min = min(pred_md, pred_lo, pred_up)*1.01
  )


pred_plots <- lapply(1:nrow(response_cb), function(i){
  
  # Model details
  selection <- response_cb[i,]
  response <- selection$response
  x_v <- selection$x_var
  
  # Extract response prediction for selected driver
  plot_df <- pred_df %>% 
    semi_join(
      selection,
      by = join_by(species,response,x_var)
    )  %>% 
    
    # Create grouping variable based on z variable and range
    mutate(
      group_idx = paste(z_cat,z_var,sep="_")
    )
  
  # Extract min and max y values
  minmax_select <- response_minmax %>% 
    semi_join(
      selection,
      by = join_by(species,response)
    ) 
  y_min <- minmax_select$y_min
  y_max <- minmax_select$y_max
    
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
names(pred_plots) <- paste(
  response_cb$species,
  response_cb$response,
  response_cb$x_var
)



# Make into grid for each species
species <- unique(stringr::str_split_fixed(mods, "_", n = 2)[,1])
lapply(species, function(s){
  
  # ordering for grid
  resp_order <- c(
    "interval_mean_wt",
    "interval_density",
    "interval_biomass_mean",
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
    nrow = 5,
    # nrow = 4,
    align = "hv",
    axis = "none",
    byrow = FALSE
  )
  
  # Export
  plot_name <- paste0("_response_pred_",s,".png")
  
  ggsave(
    filename = file.path(
      plot_dir,
      "response_prediction",
      plot_name
    ),
    plot = plot,
    width = 12,
    # height = 8,
    height = 10,
    dpi = 600
  )
  
})


#  Production components  ------------------------------------------------------

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

# Prepare total production data
prod_data <- prop_plot_df %>% 
  filter(
    species == "all",
    response == "production_mean"
  )


# Separate data into positive and negative
comp_data <-  prop_plot_df %>% 
  filter(
    species == "all",
    response %in% c("interval_mean_wt","interval_density","ptob")
  )
comp_pos <- comp_data %>% filter(prop_change_md >= 0)
comp_neg <- comp_data %>% filter(prop_change_md < 0)
comp_pos$response <- factor(
  comp_pos$response, 
  levels = rev(c("interval_density","interval_mean_wt","ptob"))
)
comp_neg$response <- factor(
  comp_neg$response, 
  levels = rev(c("interval_density","interval_mean_wt","ptob"))
)


# net_data <- a %>% 
#   group_by(x_pos) %>% 
#   summarise(
#     prop_change_md = sum(prop_change_md),
#     prop_change_up = sum(prop_change_up),
#     prop_change_lo = sum(prop_change_lo)
#   )

# Create plot
comp_plot <- ggplot() +
  
  # Postive stacked bars
  geom_col(
    data = comp_pos, 
    aes(
      x = x_pos, 
      y = prop_change_md, 
      fill = response,
      alpha = prop_change_overlap0
    ), 
    position = "stack" 
  ) +
  
  # Negative stacked bars
  geom_col(
    data = comp_neg, 
    aes(
      x = x_pos, 
      y = prop_change_md, 
      fill = response,
      alpha = prop_change_overlap0
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
      alpha = prop_change_overlap0
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
      alpha = prop_change_overlap0
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


# Species contribution plots  --------------------------------------------------


plot_cats <- data.frame(
  x_var = c("depth","depth","plt_cov_int","plt_cov_int",
            "peri_vol_int","dsldd_int"),
  z_var = c("wet_sum_365day","wet_sum_365day","pisc_index", "pisc_index",
            "pisc_index", "pisc_index"),
  z_cat = c("min","max","min","max","mean","mean")
)

plot_cats <- data.frame(
  x_var = c("depth","plt_cov_int",
            "peri_vol_int","dsldd_int"),
  z_var = c("pisc_index", "pisc_index",
            "pisc_index", "pisc_index"),
  z_cat = c("mean","mean","mean","mean")
)


lapply(1:nrow(plot_cats),function(i){
  
  x <- plot_cats$x_var[i]
  z <- plot_cats$z_var[i]
  c <- plot_cats$z_cat[i]
  
  plot_df <- pred_response %>% 
    filter(
      response == "production_mean",
      z_var == z,
      z_cat == c,
      x_var == x,
      species != "all"
    )
  
  p <-ggplot() +
    
    # Postive stacked bars
    geom_col(
      data = plot_df, 
      aes(
        x = x_raw_range, 
        y = exp(pred_md), 
        fill = species
      ), 
      position = "stack" 
    )+
    labs(title = paste(x,z,c))
  
  
})

# Annual predictors
annual_df <- theta %>% 
  filter(
    x_var == "int",
    response == "production_mean",
    species != "all",
    # z_var == "pisc_index"
    z_var == "wet_sum_365day"
  )

ggplot() +
  
  # Postive stacked bars
  geom_col(
    data = annual_df, 
    aes(
      x = z_raw_range, 
      y = exp(pred_md), 
      fill = species
    ), 
    position = "stack" 
  )



