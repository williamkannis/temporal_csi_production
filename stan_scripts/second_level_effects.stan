data {
  int<lower=0> N;                     // number of data items
  int<lower=1> L;                     // number of random effect groups
  int<lower=1> K;                     // number of predictors
  int<lower=1> J;                     // number of 2nd-level predictors
  array[N] real y;                    // outcome array
  array[N] int<lower=1, upper=L> ll;  // grouping array
  array[N] row_vector[K] x;           // predictor array
  array[L] row_vector[J] z;           // 2nd-level predictor array
}
parameters {
  matrix[J,K] gamma;                  // 2nd-level effect coefficents
  vector<lower=0>[K] tau;             // group-level effect coefficents error
  array[L] vector[K] beta_raw;        // standardized group-level deviations
  real<lower=0> sigma;                // error scale
}

transformed parameters {
  array[L] vector[K] beta;            // group-level effect coefficents

  // Site-level MVN random effects - noncentered parametrization
  for (l in 1:L)
    beta[l] = (z[l] * gamma)' + tau .* beta_raw[l];
}

model {
  
  // hyperpriors
  for (l in 1:L) beta_raw[l] ~ std_normal();
  to_vector(gamma) ~ normal(0,10);
  tau ~ cauchy(0,1);
  
  // 1st-level priors
  sigma ~ cauchy(0,1);
  
                                    
  // Likelihood
  vector[N] x_beta_ll;                // group-specific means
  for (n in 1:N) {
    x_beta_ll[n] = x[n] * beta[ll[n]];
  }
  y ~ normal(x_beta_ll,sigma);
}

