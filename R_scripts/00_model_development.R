
rm(list = ls())

library(purrr)
library(lme4)
library(MASS)
library(dplyr)
library(rstan)
library(abind)

stan_dir <- "stan_scripts"

# one level multiple regression  -----------------------------------------------

# Simulation inputs
n <- 10000
beta <- c(3,0,2,-1)
sd <- 1

# number of predictors
n_beta <- length(beta)

# Simulate random predictor data
X <- cbind(1,replicate(n_beta-1, rnorm(n)))

# Simulate response data
y <- X%*%beta +rnorm(n,sd=sd)

# Prepare data
data <- data.frame(X)
names(data) <- c("intercept",sapply(1:(n_beta-1),function(i) paste0("x",i)))
data$y <- y[,1]

# Test data
summary(lm(y~x1+x2+x3,data))

pred_data <- data %>% select(intercept,x1,x2,x3)
stan_data <- list(
  N = nrow(data),
  K = ncol(pred_data),
  x = pred_data,
  y= data$y
)

out <- stan(
  file = file.path(stan_dir,"first_level_effects.stan"),
  data = stan_data,
  iter = 5000,
  warmup = 1000,
  chains =4,
  # control = list(adapt_delta = .97),  
  cores = 4
)
print(out,pars = c("beta_mu","tau","sigma"))


# with random effect - no correlation ------------------------------------------

# Simulation inputs
group_size = 100
K <- 30
n = K*group_size
beta_mu <- c(1,0,2,-1)
tau <- c(.3,.1,.7,.4)
sd <- 1

# number of predictors
n_beta <- length(beta_mu)

# Create groupings 
group <- rep(seq_len(K), each = group_size)

# Simulate random predictor data
X <- cbind(1,replicate(n_beta-1, rnorm(n)))


# Simulate group-specific betas
beta_group <- do.call(cbind,map2(beta_mu,tau,rnorm,n = K))

# Expand betas for each row in X
beta <- beta_group[group,]

# Simulate response data
y <- rowSums(X*beta) +rnorm(n,sd=sd)

# Prepare data
data <- data.frame(X)
names(data) <- c("intercept",sapply(1:(n_beta-1),function(i) paste0("x",i)))
data$y <- y
data$group <- group

# Test data
summary(lmer(y~ x1 + x2+x3+ (1 +x1 + x2+x3 || group),data = data))

pred_data <- data %>% select(intercept,x1,x2,x3)
stan_data <- list(
  N = nrow(data),
  L = n_distinct(data$group),
  K = ncol(pred_data),
  y= data$y,
  ll = data$group,
  x = pred_data
)

out <- stan(
  file = file.path(stan_dir,"random_effects.stan"),
  data = stan_data,
  iter = 5000,
  warmup = 1000,
  chains =4,
  # control = list(adapt_delta = .97),  
  cores = 4
)
print(out,pars = c("beta_mu","tau","sigma"))


# 2nd level effects - no correlation -------------------------------------------

# Simulation input
group_size = 100
K <- 30
n = K*group_size
int <- 3
gamma <- list(  # first value in each vector is beta intercept
  c(1,0,0.5),
  c(0,1.1,-0.5),
  c(2,-1.9,-0.2),
  c(-1,2,0))
tau <- c(.3,.1,.7,.4)
sd <- 1

# number of predictors
n_beta <- length(beta_mu)
n_gamma <- length(gamma[[1]])

# Create groupings
group <- rep(seq_len(K), each = group_size)

# Simulate random predictor data
X <- cbind(1,replicate(n_beta-1, rnorm(n)))

# Simulate second level predictor data
Z <- cbind(1,replicate(n_gamma-1,rnorm(K)))

# Simulate group specific beta
beta_mu <- sapply(gamma,function(g) Z%*%g)
beta_error <- sapply(tau, rnorm, n= K, mean = 0)
beta_group <- beta_mu + beta_error

# Expand beta for each row in X
beta <- beta_group[group,]

# Simulate response data
y <- rowSums(X*beta) +rnorm(n,sd=sd)

# First level data
data <- data.frame(X)
names(data) <- c("intercept",sapply(1:(n_beta-1),function(i) paste0("x",i)))
data$y <- y
data$group <- group

# second level data
data2 <- data.frame(Z)
names(data2) <- c("intercept",sapply(1:(n_gamma-1),function(i) paste0("z",i)))
data2$group <- 1:K

# Test data
summary(lmer(y~ x1 + x2+x3+ (1 +x1 + x2+x3 || group),data = data))

pred_data <- data %>% 
  # mutate(intercept = 1) %>% 
  select(intercept,x1,x2,x3)
pred_data2 <- data2 %>% 
  # mutate(intercept = 1) %>% 
  select(intercept,z1,z2)
stan_data <- list(
  N = nrow(data),
  L = nrow(pred_data2),
  K = ncol(pred_data),
  J = ncol(pred_data2),
  y= data$y,
  ll = data$group,
  x = pred_data,
  z = pred_data2
)

out <- stan(
  file = file.path(stan_dir,"second_level_effects.stan"),
  data = stan_data,
  iter = 5000,
  warmup = 1000,
  chains =4,
  control = list(adapt_delta = .97),  
  cores = 4
)
print(out,pars = c("gamma","tau","sigma"))


# Second-level year effects, site intercept  -----------------------------------

# Simulation input
group_size = 5
K <- 30
K2 <- 20
n = K*K2*group_size
gamma <- list(  # first value in each vector is beta intercept
  c(1,0,0.5),
  c(0,1.1,-0.5),
  c(2,-1.9,-0.2),
  c(-1,2,0))
tau <- c(.3,.1,.7,.4)
tau2 <- .5
sd <- 1

# number of predictors
n_beta <- length(gamma)
n_gamma <- length(gamma[[1]])

# Create groupings
group <- rep(seq_len(K), each = group_size*K2)
group2 <- rep(seq_len(K2), group_size*K)

# Simulate random predictor data
X <- cbind(1,replicate(n_beta-1, rnorm(n)))

# Simulate second level predictor data
Z <- cbind(1,replicate(n_gamma-1,rnorm(K)))

# Simulate group specific beta
beta_mu1 <- sapply(gamma,function(g) Z%*%g)
beta_error1 <- sapply(tau, rnorm, n= K, mean = 0)
beta_error2 <- rnorm(K2,sd = tau2)
beta_group1 <- beta_mu1 + beta_error1

# Expand beta for each row in X
beta1 <- beta_group1[group,]
beta2 <- beta_error2[group2]
beta <- beta1
beta[,1] <- beta[,1] + beta2

# Simulate response data
y <- rowSums(X*beta) +rnorm(n,sd=sd)

# First level data
data <- data.frame(X)
names(data) <- c("intercept",sapply(1:(n_beta-1),function(i) paste0("x",i)))
data$y <- y
data$group1 <- group
data$group2 <- group2

# second level data
data2 <- data.frame(Z)
names(data2) <- c("intercept",sapply(1:(n_gamma-1),function(i) paste0("z",i)))
data2$group <- 1:K

# Test data
summary(lmer(y~ x1 + x2+x3+ (1 +x1 + x2+x3 || group),data = data))

pred_data <- data %>% 
  # mutate(intercept = 1) %>% 
  select(intercept,x1,x2,x3)
pred_data2 <- data2 %>% 
  # mutate(intercept = 1) %>% 
  select(intercept,z1,z2)
stan_data <- list(
  N = nrow(data),
  L = nrow(pred_data2),
  O = n_distinct(data$group2),
  K = ncol(pred_data),
  J = ncol(pred_data2),
  y= data$y,
  ll = data$group1,
  oo = data$group2,
  x = pred_data,
  z1 = pred_data2
)

out <- stan(
  file = file.path(stan_dir,"second_level_effects_site_effect.stan"),
  data = stan_data,
  iter = 5000,
  warmup = 1000,
  chains =4,
  # control = list(adapt_delta = .97),  
  cores = 4
)
print(out,pars = c("gamma1","tau1","tau2","sigma"))

# MVN Second-level year effects, site intercept  -------------------------------

# Simulation input
group_size = 5
K <- 30
K2 <- 20
n = K*K2*group_size
gamma <- list(  # first value in each vector is beta intercept
  c(1,0,0.5),
  c(0,1.1,-0.5),
  c(2,-1.9,-0.2),
  c(-1,2,0))
tau1 <- c(.3,.1,.7,.4)
tau2 <- .5
sd <- 1

# COV matrix
cor.mat = matrix(
  c(1,.3,.8,.2,
    .3,1,.5,.05,
    .8,.5,1,.3,
    .2,.05,.3,1),
  4, 
  4
)
cov.mat <-diag(tau1) %*% cor.mat  %*%  diag(tau1)

# number of predictors
n_beta <- length(gamma)
n_gamma <- length(gamma[[1]])

# Create groupings
group <- rep(seq_len(K), each = group_size*K2)
group2 <- rep(seq_len(K2), group_size*K)

# Simulate random predictor data
X <- cbind(1,replicate(n_beta-1, rnorm(n)))

# Simulate second level predictor data
Z <- cbind(1,replicate(n_gamma-1,rnorm(K)))

# Simulate group specific beta
beta_mu1 <- sapply(gamma,function(g) Z%*%g)
beta_error1 <- mvrnorm(K,mu=rep(0,n_beta),Sigma =cov.mat)
beta_error2 <- rnorm(K2,sd = tau2)
beta_group1 <- beta_mu1 + beta_error1

# Expand beta for each row in X
beta1 <- beta_group1[group,]
beta2 <- beta_error2[group2]
beta <- beta1
beta[,1] <- beta[,1] + beta2

# Simulate response data
y <- rowSums(X*beta) +rnorm(n,sd=sd)

# First level data
data <- data.frame(X)
names(data) <- c("intercept",sapply(1:(n_beta-1),function(i) paste0("x",i)))
data$y <- y
data$group1 <- group
data$group2 <- group2

# second level data
data2 <- data.frame(Z)
names(data2) <- c("intercept",sapply(1:(n_gamma-1),function(i) paste0("z",i)))
data2$group <- 1:K

# Test data
summary(lmer(y~ x1 + x2+x3+ (1 +x1 + x2+x3 || group),data = data))

pred_data <- data %>% 
  # mutate(intercept = 1) %>% 
  dplyr::select(intercept,x1,x2,x3)
pred_data2 <- data2 %>% 
  # mutate(intercept = 1) %>% 
  dplyr::select(intercept,z1,z2)
stan_data <- list(
  N = nrow(data),
  L = nrow(pred_data2),
  O = n_distinct(data$group2),
  K = ncol(pred_data),
  J = ncol(pred_data2),
  y= data$y,
  ll = data$group1,
  oo = data$group2,
  x = pred_data,
  z1 = pred_data2
)

out <- stan(
  file = file.path(stan_dir,"mvn_second_level_effects_site_effect.stan"),
  data = stan_data,
  iter = 5000,
  warmup = 1000,
  chains =4,
  # control = list(adapt_delta = .97),  
  cores = 4
)
print(out,pars = c("gamma1","tau1","tau2","sigma","cor_mat1"))


# MVN second level site-year effects -------------------------------------------
# Simulation input
group_size = 5
K1 <- 30
K2 <- 20
n = K1*K2*group_size
gamma <- list(  # first value in each vector is beta intercept
  c(1,0,0.5),
  c(0,1.1,-0.5),
  c(2,-1.9,-0.2),
  c(-1,2,0))
tau <- c(.3,.1,.7,.4)
tau1 <- c(.2,.5,.1,.2)
tau2 <- c(.5,.2,.3,.5)
sd <- 1

# COV matrix
cor.mat = matrix(
  c(1,.3,.8,.2,
    .3,1,.5,.05,
    .8,.5,1,.3,
    .2,.05,.3,1),
  4, 
  4
)
cov.mat <-diag(tau) %*% cor.mat  %*%  diag(tau)

# number of predictors
n_beta <- length(gamma)
n_gamma <- length(gamma[[1]])

# Create groupings
group1 <- rep(seq_len(K1), each = group_size*K2)
group2 <- rep(seq_len(K2), group_size*K1)
group <- data.frame(group1 = group1,group2 =group2) %>% 
  group_by(group1, group2) %>% 
  mutate(group = cur_group_id()) %>% 
  pull(group)
K = n_distinct(group)


# Simulate random predictor data
X <- cbind(1,replicate(n_beta-1, rnorm(n)))

# Simulate second level predictor data
Z <- cbind(1,replicate(n_gamma-1,rnorm(K)))

# Simulate group specific beta
beta_mu <- sapply(gamma,function(g) Z%*%g)
beta_error <- mvrnorm(K,mu=rep(0,n_beta),Sigma =cov.mat)
beta_group <- beta_mu + beta_error

beta_error1 <- sapply(tau1,rnorm,n=K1,mean = 0)
beta_error2 <- sapply(tau2,rnorm,n=K2,mean = 0)


# Expand beta for each row in X
beta_c <- beta_group[group,]
beta1 <- beta_error1[group1,]
beta2 <- beta_error2[group2,]
beta <- beta_c+beta1+beta2

# Simulate response data
y <- rowSums(X*beta) +rnorm(n,sd=sd)

# First level data
data <- data.frame(X)
names(data) <- c("intercept",sapply(1:(n_beta-1),function(i) paste0("x",i)))
data$y <- y
data$group <- group
data$group1 <- group1
data$group2 <- group2

# second level data
data2 <- data.frame(Z)
names(data2) <- c("intercept",sapply(1:(n_gamma-1),function(i) paste0("z",i)))
data2$group <- 1:K
data2 <-data2 %>% 
  left_join(
    data %>% dplyr::distinct(group,group1,group2)
  )

# Test data
summary(lmer(y~ x1 + x2+x3+ (1 +x1 + x2+x3 || group),data = data))

pred_data <- data %>% 
  # mutate(intercept = 1) %>% 
  dplyr::select(intercept,x1,x2,x3)
pred_data2 <- data2 %>% 
  # mutate(intercept = 1) %>% 
  dplyr::select(intercept,z1,z2)
stan_data <- list(
  N = nrow(data),
  L = n_distinct(data$group1),
  O = n_distinct(data$group2),
  P = n_distinct(data$group),
  K = ncol(pred_data),
  J = ncol(pred_data2),
  y= data$y,
  pp = data$group,
  ll = data2$group1,
  oo = data2$group2,
  x = pred_data,
  z = pred_data2
)

out <- stan(
  file = file.path(stan_dir,"mvn_second_level_siteyear_effects.stan"),
  data = stan_data,
  iter = 7000,
  warmup = 1000,
  chains =4,
  # control = list(adapt_delta = .97),  
  cores = 4
)
print(out,pars = c("gamma","tau","tau1","tau2","sigma","cor_mat1"))


# MVN region-year effects  -----------------------------------------------------
# Simulation input
group_size = 5
nest_size = 7
K1 <- 30
K3 <- 3
K2 <- nest_size*K3
n = K1*K2*group_size
gamma <- list(  # first value in each vector is beta intercept
  c(1,0,0.5),
  c(0,1.1,-0.5),
  c(2,-1.9,-0.2),
  c(-1,2,0))
tau <- c(.3,.1,.7,.4)
tau2 <- c(.5,.2,.3,.5)
sd <- 1

# COV matrix
cor.mat = matrix(
  c(1,.3,.8,.2,
    .3,1,.5,.05,
    .8,.5,1,.3,
    .2,.05,.3,1),
  4, 
  4
)
cov.mat <-diag(tau) %*% cor.mat  %*%  diag(tau)
# cor.mat2 = matrix(
#   c(1,.3,.4,.2,
#     .3,1,.8,.7,
#     .4,.8,1,.1,
#     .2,.7,.1,1),
#   4, 
#   4
# )
cov.mat2 <-diag(tau2) %*% cor.mat  %*%  diag(tau2)

# number of predictors
n_beta <- length(gamma)
n_gamma <- length(gamma[[1]])

# Create groupings
group1 <- rep(seq_len(K1), each = group_size*K2)
group2 <- rep(seq_len(K2), group_size*K1)
# nest_bridge <- cbind(
#   group2 =seq_len(K2),
#   group3 = rep(seq_len(K3),each = nest_size)
#   )
group3 <- rep(seq_len(K3),each = nest_size)
nest_bridge <- group3[group2]

# Simulate random predictor data
X <- cbind(1,replicate(n_beta-1, rnorm(n)))

# Simulate second level predictor data for each region
Z_list <- replicate(K3,cbind(1,replicate(n_gamma-1,rnorm(K1))),simplify=F)

# Simulate region-year specific beta
beta1_group <- lapply(Z_list, function(Z){
  beta_mu <- sapply(gamma,function(g) Z%*%g)
  beta_error <- mvrnorm(K1,mu=rep(0,n_beta),Sigma =cov.mat)
  beta_mu + beta_error
})
beta2_error <- mvrnorm(K2,mu=rep(0,n_beta),Sigma =cov.mat)




# Expand beta for each row in X
beta1 <- do.call(
  rbind,
  map2(
    nest_bridge,
    group1,
    function(g3,g1) beta1_group[[g3]][g1,]
    )
  )
beta2 <- beta2_error[group2,]
beta <- beta1+beta2

# Simulate response data
y <- rowSums(X*beta) +rnorm(n,sd=sd)

# First level data
data <- data.frame(X)
names(data) <- c("intercept",sapply(1:(n_beta-1),function(i) paste0("x",i)))
data$y <- y
data$group1 <- group1
data$group2 <- group2


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
  dplyr::select(intercept,x1,x2,x3)

stan_data <- list(
  N = nrow(data),
  `T` = n_distinct(data$group1),
  S = n_distinct(data$group2),
  R = n_distinct(group3),
  K = ncol(x_data),
  L = ncol(z_data[1,,]),
  y= data$y,
  yr = data$group1,
  st = data$group2,
  rg = group3,
  x = x_data,
  z = z_data
)

out <- stan(
  file = file.path(stan_dir,"mvn_second_level_regyear_effects.stan"),
  data = stan_data,
  iter = 5000,
  warmup = 1000,
  chains =4,
  # control = list(adapt_delta = .97),  
  cores = 4
)
print(out,pars = c("gamma","tau","tau_s","tau_r","sigma","cor_mat","cor_mat_s"))
