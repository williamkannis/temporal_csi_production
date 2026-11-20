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