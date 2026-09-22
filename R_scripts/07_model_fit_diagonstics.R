
# House keeping  ---------------------------------------------------------------
rm(list = ls())

# Packages


# Directories
out_dir <- "stan_outputs"
data_dir <- "hpc/data"


# Load in model out puts
out_files <- list.files(out_dir,"M100.rds$")
out_list <- lapply(out_files, function(x) readRDS(file.path(out_dir,x)))
names(out_list) <- gsub("_stan_out_M100.rds","",out_files)

# load in data
data_files <- list.files(data_dir)
data_list <- lapply(data_files, function(x) readRDS(file.path(data_dir,x)))
names(data_list) <- gsub("_input_data.rds","",data_files)


for (i in 1:length(out_list)) {
  
  out <- out_list[[i]]
  name <- names(out_list)[[i]]
  data <- data_list[[name]]
  y_rep <- out$draws("y_rep", format = "matrix")
  
 print( bayesplot::ppc_ecdf_overlay(
    data$y,
    y_rep[1:100,]
  ) +ggplot2::ggtitle(name))
  
  print(bayesplot::ppc_dens_overlay(
    data$y,
    y_rep[1:100,]
  ) +ggplot2::ggtitle(name)) #+  ggplot2::coord_cartesian(xlim = c(0, 500))
  
print(  bayesplot::ppc_stat(
    data$y,
    y_rep,
    stat = function(x) mean(x == 0)
  ) +ggplot2::ggtitle(name))
  
print(  bayesplot::ppc_stat(
    data$y,
    y_rep,
    stat = mean
  ) +ggplot2::ggtitle(name))

print(  bayesplot::ppc_stat(
    data$y,
    y_rep,
    stat = var
  ) +ggplot2::ggtitle(name))
  
print(  bayesplot::ppc_stat(
    data$y,
    y_rep,
    stat = max
  ) +ggplot2::ggtitle(name))
  
  
}
param <- c(
  "gamma",
  "tau",
  "tau_s",
  "phi",
  "theta"
      )
a <- lapply(out_list, function(x){
  x$summary(
    param,
    "mean",
    "rhat",
    quantile, 
    .args = list(probs = c(0.025, 0.975))
    ) %>% 
    mutate(
      overlap0 = case_when(
        `2.5%` < 0 & `97.5%` > 0 ~ T,
        T ~ F
      )
      )
  
})

lapply(1:length(out_list), function(x){
  out <- out_list[[x]]
  name <- names(out_list)[[x]]
  out$summary(c("phi","theta"),mean) %>% 
    mutate(mod = name)
  }
  ) %>% 
  bind_rows() %>% 
  tidyr::pivot_wider(names_from = variable, values_from = mean)

b <- lapply(1:length(out_list), function(x){
  out <- out_list[[x]]
  name <- names(out_list)[[x]]
  out$summary("gamma",mean,quantile, 
              .args = list(probs = c(0.025, 0.975))
  ) %>% 
    mutate(
      overlap0 = case_when(
        `2.5%` < 0 & `97.5%` > 0 ~ T,
        T ~ F
      )
    ) %>% 
    mutate(
      mod = name,
      `mean` = round(`mean`,2),
      mean_ci = case_when(
        overlap0 == T ~as.character(`mean`),
        overlap0 == F ~ paste0(`mean`,"*"),
        T ~ NA
      )
      ) %>% 
    select(mod,variable,mean_ci)
}
) %>% 
  bind_rows() %>% 
  tidyr::pivot_wider(names_from = variable, values_from = mean_ci)


for (i in 1:length(out_list)) {
lambda_draws <- posterior::as_draws_matrix(out_list[[i]]$draws("lambda"))

max_lambda <- quantile(lambda_draws, 0.999)

M_required <- qpois(
  1 - 1e-10,
  lambda = max_lambda
)

lambda_max <- quantile(lambda_draws, 0.999)
print("==================================================")
print(names(out_list)[i])
print("==================================================")
print(data.frame(
  tolerance = c(1e-6, 1e-8, 1e-10, 1e-12),
  M = c(
    qpois(1 - 1e-6, lambda_max),
    qpois(1 - 1e-8, lambda_max),
    qpois(1 - 1e-10, lambda_max),
    qpois(1 - 1e-12, lambda_max)
  )
))
M=60
poisson_tail <- 1 - ppois(
  M,
  lambda = lambda_draws
)

print(quantile(
  poisson_tail,
  probs = c(0.5, 0.9, 0.99, 0.999, 1)
))
}


for(i in 1:length(out_list)) {
  out <- out_list[[i]]
  name <- names(out_list)[[i]]
  
  # convergence
  conv <- out$summary(NULL, c("rhat", "ess_bulk"))
  high_rhat <- sum(conv$rhat > 1.1, na.rm = TRUE)
  low_ess   <- sum(conv$ess_bulk < 400, na.rm = TRUE)
  
  # sampling
  diag <- out$diagnostic_summary()
  div <- sum(diag$num_divergent)
  tree <- sum(diag$num_max_treedepth)
  
  # combined issues
  issues <- div + high_rhat + low_ess
  
  if (issues == 0) {
    message("============================================================")
    message(name)
    message("No issues dectected:")
    message("  Rhat > 1.1: ", high_rhat)
    message( "Bulk ESS < 400: ", low_ess)
    message("  Divergent transitions: ", div)
    message(" Max tree depth exceeded: ", tree)
    message("============================================================")
  } else {
    message("============================================================")
    message("WARNING")
    message(name)
    message("Issues detected:")
    message("  Rhat > 1.1: ", high_rhat)
    message("  Divergent transitions: ", div)
    message(" Max tree depth exceeded: ", tree)
    message("============================================================")
  }
}


b <-out$draws("gamma")
c <- posterior::as_draws_rvars(b)

diag <- tweedie_diagnostics(
  fit = out_list$all_production_mean,
  y = data_list$all_production_mean$y,
  st = data_list$all_production_mean$st,
  yr = data_list$all_production_mean$yr,
  rg = data_list$all_production_mean$rg
)

hist(diag$group_ppc$summary$p_zero)
diag$plots$group_mean
diag$plots$zero_probability_group

diag_list <- parallel::mclapply(
  mc.cores = 7,
  X = seq_along(out_list),
  FUN = function(x){
    
    out <- out_list[[x]]
    data <- data_list[[names(out_list)[x]]]
    
    tweedie_diagnostics(
      fit = out,
      y = data$y,
      st = data$st,
      yr = data$yr,
      rg = data$rg
    )
  }
  )
diag$plots$rqr_qq
seq_along(diag_list)
for(i in seq_along(diag_list)) {
  # print(diag_list[[i]]$plots$rqr_qq + ggplot2::ggtitle(names(out_list)[i]))
  print(diag_list[[i]]$plots$rqr_vs_fitted + ggplot2::ggtitle(names(out_list)[i]))
  # print(diag_list[[i]]$plots$residual_order + ggplot2::ggtitle(names(out_list)[i]))
  # print(diag_list[[i]]$plots$zero_probability + ggplot2::ggtitle(names(out_list)[i]))
  # print(diag_list[[i]]$plots$group_mean + ggplot2::ggtitle(names(out_list)[i]))
  # print(diag_list[[i]]$plots$zero_probability_group + ggplot2::ggtitle(names(out_list)[i]))
}


st_raw <- posterior::as_draws_matrix(
  out$draws("st_raw")
)

st_summary <- data.frame(
  
  site = seq_len(ncol(st_raw)),
  
  mean = apply(
    st_raw,
    2,
    mean
  ),
  
  sd = apply(
    st_raw,
    2,
    sd
  ),
  
  q025 = apply(
    st_raw,
    2,
    quantile,
    0.025
  ),
  
  q975 = apply(
    st_raw,
    2,
    quantile,
    0.975
  )
)

ggplot(
  st_summary,
  aes(
    x = mean,
    y = reorder(site, mean)
  )
) +
  
  geom_errorbarh(
    aes(
      xmin = q025,
      xmax = q975
    )
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  
  labs(
    x = "Standardized site effect",
    y = "Site"
  ) +
  
  theme_classic()


beta_raw <- posterior::as_draws_matrix(
  out$draws("beta_raw")
)
beta_raw[1:10,]