rm(list = ls())
library(MASS)
library(purrr)
library(dplyr)
library(abind)
library(cmdstanr)
library(clusterGeneration)
stan_dir <- "stan_scripts"
stan_dir <- "hpc/stan_scripts"


# Hurdle predictors  -----------------------------------------------------------
rm(list=ls())
# Simulation input
group_size = 5
nest_size = 7
K_T <- 30
K_R <- 3
K_S <- nest_size*K_R
n = K_T*K_S*group_size
gamma <- list(  # first value in each vector is beta intercept
  c(-.3,.3,0.04),
  c(0,-0.07,0.02),
  c(.3,.08,-0.1),
  c(-.04,-.3,.1),
  c(.1,.1,0.02))
beta_hurdle <- c(-1.5,.1,-.3,.7,-.6)
tau_T <- c(.41,.15,.46,.21,.21)
tau_S <- .69

# theta <- .6
# sd <- .001
shape <- 2

# COV matrix
cor.mat = matrix(
  c(1,.3,.8,.2,
    .3,1,.5,.05,
    .8,.5,1,.3,
    .2,.05,.3,1),
  4, 
  4
)
cor.mat = rcorrmatrix(d = 5)
cov.mat <-diag(tau_T) %*% cor.mat  %*%  diag(tau_T)


# number of predictors
n_beta <- length(gamma)
n_beta_h <- length(beta_hurdle)
n_gamma <- length(gamma[[1]])

# Create groupings
group_t <- rep(seq_len(K_T), each = group_size*K_S)
group_s <- rep(seq_len(K_S), group_size*K_T)
group_r <- rep(seq_len(K_R),each = nest_size)
reg_bridge <- group_r[group_s]

# Simulate random predictor data
X <- cbind(1,replicate(n_beta-1, rnorm(n)))

# Simulate random hurdle predictor data
X_hurdle <- cbind(1,replicate(n_beta_h-1, rnorm(n)))

# Simulate second level predictor data for each region
Z_list <- replicate(K_R,cbind(1,replicate(n_gamma-1,rnorm(K_T))),simplify=F)

# Simulate region-year specific beta
beta_tr_group <- lapply(Z_list, function(Z){
  beta_mu <- sapply(gamma,function(g) Z%*%g)
  beta_error <- mvrnorm(K_T,mu=rep(0,n_beta),Sigma =cov.mat)
  beta_mu + beta_error
})


# Expand beta for each row in X
beta <- do.call(
  rbind,
  map2(
    reg_bridge,
    group_t,
    function(r,t) beta_tr_group[[r]][t,]
  )
)

# Create site error
site_error <- rnorm(K_S,0,tau_S)
site_error_full <- site_error[group_s]

# Simulate zeros
theta_logit <- X_hurdle %*% beta_hurdle
theta <- plogis(theta_logit)
hurdle <- rbernoulli(n,p =theta)

# Simulate response data
mu <- rowSums(X*beta) + site_error_full
# scale <- exp(mu)/shape
rate <- shape/exp(mu)
y <- rep(0,n)
# y[hurdle] <- rlnorm(sum(hurdle),meanlog = mu[hurdle], sdlog = sd)
y[!hurdle] <- rgamma(
  n = sum(!hurdle),
  shape = shape,
  # scale = scale,
  rate = rate[!hurdle]
)


# First level data
data <- data.frame(X)
names(data) <- c("intercept",sapply(1:(n_beta-1),function(i) paste0("x",i)))
data$y <- y
data$group_t <- group_t
data$group_s <- group_s

# # filter out zeros
# data <- data[hurdle,]
# data$group_r <-reg_bridge[hurdle]
# regyear_no0 <- data %>% distinct(group_t,group_r)


# second level data
data2_list <- lapply(Z_list, function(Z) {
  df <- data.frame(Z)
  names(df) <- c("intercept",sapply(1:(n_gamma-1),function(i) paste0("z",i)))
  df
}
)
# ## REMOVE ZERO
# data2_list <- lapply(1:3, function(z){
#   df <- data2_list[[z]]
#   years <- regyear_no0$group_t[regyear_no0$group_r == z]
#   years <- years[order(years)]
#   df[years,]
# })

z_bind <- abind(data2_list,along=3)
z_data <- aperm(z_bind,c(3,1,2))


x_data <- data %>% 
  # mutate(intercept = 1) %>% 
  dplyr::select(intercept,x1,x2,x3,x4)


stan_data <- list(
  N = nrow(data),
  `T` = n_distinct(data$group_t),
  S = n_distinct(data$group_s),
  R = n_distinct(group_r),
  K = ncol(x_data),
  H = ncol(X_hurdle),
  L = ncol(z_data[1,,]),
  y = data$y,
  yr = data$group_t,
  st = data$group_s,
  rg = group_r,
  x = x_data,
  x_hurdle = X_hurdle,
  z = z_data
)

mod <- cmdstan_model(
  "stan_scripts/hurdle_mvn_second_level_regyear_effects_h_pred.stan"
)
# mod <- cmdstan_model(
#   "stan_scripts/gamma_mvn_second_level_regyear_effects.stan"
# )

out <- mod$sample(
  data = stan_data,
  iter_sampling = 2000,
  iter_warmup = 1000,
  chains = 4,
  parallel_chains = 4,
  init = 0.5
)

out$diagnostic_summary()


y_rep <- out$draws("y_rep", format = "matrix")
bayesplot::ppc_dens_overlay(
  stan_data$y,
  y_rep[1:100,]
) +  ggplot2::coord_cartesian(xlim = c(0, 20))

# loo::loo(loo::extract_log_lik(out))


real <- c(do.call(cbind,gamma),beta_hurdle,tau_T,tau_S,shape)
# real <- c(do.call(cbind,gamma),tau_T,tau_S,shape)
means <- out$summary(
  c("gamma","beta_hurdle","tau","tau_s","shape"),
  "mean",
  quantile, 
  .args = list(probs = c(0.025, 0.975))
)
cbind(means,real) %>% 
  mutate(overlap = case_when(
    real >= `2.5%` & real <= `97.5%` ~ T,
    T~F
  ))



# MVN region-year effects  -----------------------------------------------------
# Simulation input
group_size = 5
nest_size = 7
K_T <- 30
K_R <- 3
K_S <- nest_size*K_R
n = K_T*K_S*group_size
gamma <- list(  # first value in each vector is beta intercept
  c(-3,.3,0.04),
  c(0,-0.07,0.02),
  c(.3,.08,-0.1),
  c(-.04,-.3,.1),
  c(.1,.1,0.02))
tau_T <- c(.41,.15,.46,.21,.21)
tau_S <- .69

theta <- .6
# sd <- .001
shape <- 30

# COV matrix
cor.mat = matrix(
  c(1,.3,.8,.2,
    .3,1,.5,.05,
    .8,.5,1,.3,
    .2,.05,.3,1),
  4, 
  4
)
cor.mat = rcorrmatrix(d = 5)
cov.mat <-diag(tau_T) %*% cor.mat  %*%  diag(tau_T)


# number of predictors
n_beta <- length(gamma)
n_gamma <- length(gamma[[1]])

# Create groupings
group_t <- rep(seq_len(K_T), each = group_size*K_S)
group_s <- rep(seq_len(K_S), group_size*K_T)
group_r <- rep(seq_len(K_R),each = nest_size)
reg_bridge <- group_r[group_s]

# Simulate random predictor data
X <- cbind(1,replicate(n_beta-1, rnorm(n)))

# Simulate second level predictor data for each region
Z_list <- replicate(K_R,cbind(1,replicate(n_gamma-1,rnorm(K_T))),simplify=F)

# Simulate region-year specific beta
beta_tr_group <- lapply(Z_list, function(Z){
  beta_mu <- sapply(gamma,function(g) Z%*%g)
  beta_error <- mvrnorm(K_T,mu=rep(0,n_beta),Sigma =cov.mat)
  beta_mu + beta_error
})


# Expand beta for each row in X
beta <- do.call(
  rbind,
  map2(
    reg_bridge,
    group_t,
    function(r,t) beta_tr_group[[r]][t,]
  )
)

# Create site error
site_error <- rnorm(K_S,0,tau_S)
site_error_full <- site_error[group_s]

# Simuate zeros
hurdle <- rbernoulli(n,p =theta)

# Simulate response data
mu <- rowSums(X*beta) + site_error_full
# scale <- exp(mu)/shape
rate <- shape/exp(mu)
y <- rep(0,n)
# y[hurdle] <- rlnorm(sum(hurdle),meanlog = mu[hurdle], sdlog = sd)
y[hurdle] <- rgamma(
  n = sum(hurdle),
  shape = shape,
  # scale = scale,
  rate = rate[hurdle]
  )


# First level data
data <- data.frame(X)
names(data) <- c("intercept",sapply(1:(n_beta-1),function(i) paste0("x",i)))
data$y <- y
data$group_t <- group_t
data$group_s <- group_s

# # filter out zeros
# data <- data[hurdle,]
# data$group_r <-reg_bridge[hurdle]
# regyear_no0 <- data %>% distinct(group_t,group_r)


# second level data
data2_list <- lapply(Z_list, function(Z) {
  df <- data.frame(Z)
  names(df) <- c("intercept",sapply(1:(n_gamma-1),function(i) paste0("z",i)))
  df
}
)
# ## REMOVE ZERO
# data2_list <- lapply(1:3, function(z){
#   df <- data2_list[[z]]
#   years <- regyear_no0$group_t[regyear_no0$group_r == z]
#   years <- years[order(years)]
#   df[years,]
# })

z_bind <- abind(data2_list,along=3)
z_data <- aperm(z_bind,c(3,1,2))


x_data <- data %>% 
  # mutate(intercept = 1) %>% 
  dplyr::select(intercept,x1,x2,x3,x4)





stan_data <- list(
  N = nrow(data),
  `T` = n_distinct(data$group_t),
  S = n_distinct(data$group_s),
  R = n_distinct(group_r),
  K = ncol(x_data),
  L = ncol(z_data[1,,]),
  y = data$y,
  yr = data$group_t,
  st = data$group_s,
  rg = group_r,
  x = x_data,
  z = z_data
)

mod <- cmdstan_model(
  "stan_scripts/hurdle_mvn_second_level_regyear_effects.stan"
  )
# mod <- cmdstan_model(
#   "stan_scripts/gamma_mvn_second_level_regyear_effects.stan"
# )

out <- mod$sample(
  data = stan_data,
  iter_sampling = 2000,
  iter_warmup = 1000,
  chains = 4,
  parallel_chains = 4,
  init = 0.5
)

out$diagnostic_summary()


y_rep <- out$draws("y_rep", format = "matrix")
bayesplot::ppc_dens_overlay(
  stan_data$y,
  y_rep[1:100,]
) +  ggplot2::coord_cartesian(xlim = c(0, 20))

# loo::loo(loo::extract_log_lik(out))






real <- c(do.call(cbind,gamma),tau_T,tau_S,shape,1-theta)
# real <- c(do.call(cbind,gamma),tau_T,tau_S,shape)
means <- out$summary(
  # c("gamma","tau","tau_s","shape"),
  c("gamma","tau","tau_s","shape","hu"),
  "mean",
  quantile, 
  .args = list(probs = c(0.025, 0.975))
  )
cbind(means,real) %>% 
  mutate(overlap = case_when(
    real >= `2.5%` & real <= `97.5%` ~ T,
    T~F
  ))


### Simple model  --------------------------------------------------------------
# rm(list=ls())
# Simulation input
n = 1000
gamma <- c(-0.30,  0.00,  0.30, -0.04,  0.10)
theta <- .6
shape <- 2




# number of predictors
n_gamma <- length(gamma)



# Simulate random predictor data
X <- cbind(1,replicate(n_gamma-1, rnorm(n)))

# Simulate zeros
hurdle <- rbernoulli(n,p =theta)

# Simulate response data
mu <- X%*%gamma
# scale <- exp(mu)/shape
rate <- shape/exp(mu)
y <- rep(0,n)
# y[hurdle] <- rlnorm(sum(hurdle),meanlog = mu[hurdle], sdlog = sd)
y[hurdle] <- rgamma(
  n = sum(hurdle),
  shape = shape,
  # scale = scale,
  rate = rate[hurdle]
)




# First level data
data <- data.frame(X)
names(data) <- c("intercept",sapply(1:(n_gamma-1),function(i) paste0("x",i)))
data$y <- y


# # filter out zeros
# data <- data[hurdle,]
# data$group_r <-reg_bridge[hurdle]
# regyear_no0 <- data %>% distinct(group_t,group_r)

x_data <- data %>% 
  # mutate(intercept = 1) %>% 
  dplyr::select(intercept,x1,x2,x3,x4)


stan_data <- list(
  N = nrow(data),
  K = ncol(x_data),
  y = data$y,
  x = x_data
)

mod <- cmdstan_model(
  "stan_scripts/hurdle_simple.stan"
  )
# mod <- cmdstan_model(
#   "stan_scripts/gamma_mvn_second_level_regyear_effects.stan"
# )

out <- mod$sample(
  data = stan_data,
  iter_sampling = 2000,
  iter_warmup = 1000,
  chains = 4,
  parallel_chains = 4,
  init = 0.5
)

out$diagnostic_summary()


# y_rep <- out$draws("y_rep", format = "matrix")
# bayesplot::ppc_dens_overlay(
#   stan_data$y,
#   y_rep[1:100,]
# )
# 
# loo::loo(loo::extract_log_lik(out))






# real <- c(do.call(cbind,gamma),tau_T,tau_S,shape,theta)
real <- c(gamma,shape,theta)
means <- out$summary(
  c("gamma","shape","hu"),
  # c("gamma","tau","tau_s","shape","hu"),
  "mean",
  quantile, 
  .args = list(probs = c(0.025, 0.975))
)
cbind(means,real) %>% 
  mutate(overlap = case_when(
    real >= `2.5%` & real <= `97.5%` ~ T,
    T~F
  ))

