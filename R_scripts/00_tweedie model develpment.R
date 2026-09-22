

## TO DO 
# 1. TEST AT MORE VALUES OF PHI AND POWER
# 2. come up with diagonsitc for M
# 3. try and use better priors to reduce processing time

rm(list = ls())
library(tweedie)
library(MASS)
library(purrr)
library(dplyr)
library(abind)
library(rstan)
library(clusterGeneration)
stan_dir <- "stan_scripts"
stan_dir <- "hpc/stan_scripts"

hist(tweedie::rtweedie(10000,power=1.1, mu=.75, phi=9),breaks = 40)
# ?rtweedie
# 
# hist(rpois(1000,10000000))


# MVN region-year effects  -----------------------------------------------------
# Simulation input
group_size = 5
nest_size = 7
K_T <- 10
K_R <- 3
K_S <- nest_size*K_R
n = K_T*K_S*group_size
gamma <- list(  # first value in each vector is beta intercept
  c(5.27,.3,0.04),
  c(0,-0.07,0.02),
  c(.3,.08,-0.1),
  c(-.04,-.3,.1),
  c(.1,.1,0.02))
tau_T <- c(.41,.15,.46,.21,.21)
tau_S <- .69

power <- 1.5
phi <- 8

power <- 1.4
phi <- .26

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

# Simulate response data
mu <- exp(rowSums(X*beta) + site_error_full)
y <- rtweedie(n,mu = mu, phi = phi, power = power)

# First level data
data <- data.frame(X)
names(data) <- c("intercept",sapply(1:(n_beta-1),function(i) paste0("x",i)))
data$y <- y
data$group_t <- group_t
data$group_s <- group_s


# second level data
data2_list <- lapply(Z_list, function(Z) {
  df <- data.frame(Z)
  names(df) <- c("intercept",sapply(1:(n_gamma-1),function(i) paste0("z",i)))
  df
}
)
z_bind <- abind(data2_list,along=3)
z_data <- aperm(z_bind,c(3,1,2))

x_data <- data %>% 
  # mutate(intercept = 1) %>% 
  dplyr::select(intercept,x1,x2,x3,x4)

stan_data <- list(
  M = 100,
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

library(cmdstanr)
mod <- cmdstan_model(file.path(stan_dir,"tweedie_mvn_regyear_effects.stan"))

out <- mod$sample(
  data = stan_data,
  iter_sampling = 2000,
  iter_warmup = 1000,
  chains = 4,
  parallel_chains = 4,
  init = 0.5
)

out110$summary()
out110$diagnostic_summary()

# out110i <- stan(
#   file = file.path(stan_dir,"tweedie_mvn_second_level_regyear_effects_opt.stan"),
#   data = stan_data,
#   iter = 3000,
#   warmup = 1000,
#   chains =4,
#   # control = list(adapt_delta = .97),  
#   cores = 4
# )
# new version ran in 1106.94 secs
# new version ran in 1790.61 secs
print(out,pars = c("gamma","tau","tau_s","phi","theta"))
print(out2,pars = c("gamma","tau","tau_s","phi","theta"))

out110$summary(NULL,"mean")
post <- out110$draws(c("max_lambda","max_eta")) %>% posterior::as_draws_df()
post <- out110$draws(c("y_rep","mu","residuals")) %>% posterior::as_draws_df()
post <- extract(out110i,c("max_lambda","max_eta"))
quantile(post$max_lambda,
         c(0, .5, .9, .95, .99, .999, 1))

quantile(post$max_eta,
         c(0, .5, .9, .95, .99, .999, 1))

post <- extract(out110i,c("y_rep","mu","residuals"))
bayesplot::ppc_ecdf_overlay(
  stan_data$y,
  post$y_rep[1:100,]
)
y_rep <- out110$draws("y_rep", format = "matrix")
bayesplot::ppc_dens_overlay(
  stan_data$y,
  y_rep[1:100,]
)

loo::loo(loo::extract_log_lik(out110))


lambda_draws <- extract(out110i,"lambda")[[1]]

max_lambda <- quantile(lambda_draws, 0.999)

M_required <- qpois(
  1 - 1e-10,
  lambda = max_lambda
)

lambda_max <- quantile(lambda_draws, 0.999)

data.frame(
  tolerance = c(1e-6, 1e-8, 1e-10, 1e-12),
  M = c(
    qpois(1 - 1e-6, lambda_max),
    qpois(1 - 1e-8, lambda_max),
    qpois(1 - 1e-10, lambda_max),
    qpois(1 - 1e-12, lambda_max)
  )
)
M=110
poisson_tail <- 1 - ppois(
  M,
  lambda = lambda_draws
)

quantile(
  poisson_tail,
  probs = c(0.5, 0.9, 0.99, 0.999, 1)
)

real <- c(do.call(rbind,gamma),tau_T,tau_S,phi,power)
means1 <-rstan::get_posterior_mean(out, c("gamma","tau","tau_s","phi","theta"))[,5]
means2 <-rstan::get_posterior_mean(out2, c("gamma","tau","tau_s","phi","theta"))[,5]

rbind(means1,means2,real)
rbind(means1,real)

real <- c(do.call(rbind,gamma),tau_T,tau_S,phi,power)
means30 <-rstan::get_posterior_mean(out30, c("gamma","tau","tau_s","phi","theta"))[,5]
means40 <-rstan::get_posterior_mean(out40, c("gamma","tau","tau_s","phi","theta"))[,5]
means50 <-rstan::get_posterior_mean(out50, c("gamma","tau","tau_s","phi","theta"))[,5]
means60 <-rstan::get_posterior_mean(out60, c("gamma","tau","tau_s","phi","theta"))[,5]
means80 <-rstan::get_posterior_mean(out80, c("gamma","tau","tau_s","phi","theta"))[,5]
means110 <-rstan::get_posterior_mean(out110i, c("gamma","tau","tau_s","phi","theta"))[,5]
means130 <-rstan::get_posterior_mean(out130, c("gamma","tau","tau_s","phi","theta"))[,5]
rbind(means30,means40,means50,means60,means80,means110,means130,real)
rbind(means110,real)

a <- out$summary(c("gamma","tau","tau_s","phi","theta"),"mean")
b <- a$mean
names(b) <- a$variable

rbind(b,real)

# simple model  ----------------------------------------------------------------

# Simulation inputs
n <- 2500
n <- 1000
power <- 1.3
phi <- 1
beta <- c(3,0,.2,-.1)


# number of predictors
n_beta <- length(beta)

# Simulate random predictor data
X <- cbind(1,replicate(n_beta-1, rnorm(n)))

# Simulate response data
mu <- exp(X%*%beta) 
Y <- rtweedie(n,mu=mu,power = power,phi=phi)

# Prepare data
stan_data <- list(
  M = 30,
  N = nrow(Y),
  K = ncol(X),
  X = X,
  Y= Y
)

out <- stan(
  file = file.path(stan_dir,"first_level_tweedie.stan"),
  data = stan_data,
  iter = 3000,
  warmup = 1000,
  chains =4,
  # control = list(adapt_delta = .97),  
  cores = 4
)

print(out,pars = c("theta","phi","beta"))


