data {
  int<lower=0> N;                     // number of data items
  int<lower=1> L;                     // random effect group1 size
  int<lower=1> O;                     // random effect group2 size
  int<lower=1> K;                     // number of predictors
  int<lower=1> J;                     // number of group1 predictors
  array[N] real y;                    // outcome array
  array[N] int<lower=1, upper=L> ll;  // group1 array
  array[N] int<lower=1, upper=O> oo;  // group1 array
  array[N] row_vector[K] x;           // predictor array
  array[L] row_vector[J] z1;          // 2nd-level predictor array - group1
}
parameters {
  matrix[J,K] gamma1;                 // group1 2nd-level effect coefficents
  vector<lower=0>[K] tau1;            // group1 effect coefficents error
  real<lower=0> tau2;                 // group2 effect coefficents error
  array[L] vector[K] beta_raw1;       // standardized group1-level deviations
  array[O] real beta_raw2;            // standardized group2-level deviations
  real<lower=0> sigma;                // error scale
}

transformed parameters {
  array[L] vector[K] beta1;           // group1-level effect coefficents
  array[O] real beta2;                // group2-level effect coefficents
  
  
  // Group-level MVN random effects - noncentered parametrization
  for (l in 1:L) beta1[l] = (z1[l] * gamma1)' + tau1 .* beta_raw1[l];
  for (o in 1:O) beta2[o] = tau2 .* beta_raw2[o];
}

model {
  
  // hyperpriors
  for (l in 1:L) beta_raw1[l] ~ std_normal();
  for (o in 1:O) beta_raw2[o] ~ std_normal();
  to_vector(gamma1) ~ normal(0,10);
  tau1 ~ cauchy(0,1);
  tau2 ~ cauchy(0,1);
  
  // 1st-level priors
  sigma ~ cauchy(0,1);
  
  // Likelihood
  vector[N] ll_mu;                    // group1-specific means
  vector[N] oo_mu;                    // group2-specific means
  for (n in 1:N) {
    ll_mu[n] = x[n] * beta1[ll[n]];
    oo_mu[n] = beta2[oo[n]];
  }
  y ~ normal(ll_mu+oo_mu,sigma);
}

