data {
  int<lower=0> N;                         // number of data items
  int<lower=1> L;                         // number of random effect groups
  int<lower=1> K;                         // number of predictors
  int<lower=1> J;                         // number of 2nd-level predictors
  array[N] real y;                        // outcome array
  array[N] int<lower=1, upper=L> ll;      // grouping array
  array[N] row_vector[K] x;               // predictor array
  array[L] row_vector[J] z;               // 2nd-level predictor array
}
parameters {
  matrix[J,K] gamma;                      // 2nd-level effect coefficents
  array[K] real<lower=0> tau;             // effect coefficents error
  array[L] vector[K] beta;                // group-level effect coefficents
  real<lower=0> sigma;                    // error scale
}
model {
  
  // hyperpriors
  to_vector(gamma) ~ normal(0,10);
  tau ~ cauchy(0,1);
  
  // 1st-level priors
  array[L] vector[K] beta_mu;
  for(l in 1:L){
    beta_mu[l] = (z[l] * gamma)';
    beta[l] ~ normal(beta_mu[l],tau);
  }
  sigma ~ cauchy(0,1);
  

  // Likelihood
  vector[N] x_beta_ll;
  for (n in 1:N) {
    x_beta_ll[n] = x[n] * beta[ll[n]];
  }
  y ~ normal(x_beta_ll,sigma);
}

