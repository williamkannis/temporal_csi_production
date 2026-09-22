data {
  int<lower=0> N;                     // number of data items
  int<lower=1> T;                     // number of years
  int<lower=1> S;                     // number of sites
  int<lower=1> R;                     // number of regions
  int<lower=1> K;                     // number of predictors
  int<lower=1> H;                     // number of hurdle predictors
  int<lower=1> L;                     // number of 2nd-level predictors
  
  array[N] int<lower=1, upper=T> yr;  // year index
  array[N] int<lower=1, upper=S> st;  // site 
  array[S] int<lower=1, upper=R> rg;  // region-site bridge
  array[N] real y;                    // outcome array
  array[N] row_vector[K] x;           // period-level predictor data
  matrix[N,H] x_hurdle;               // period-level hurdle predictor data
  array[R] matrix[T, L] z;            // year-level predictor data, by region
}

parameters {
  matrix[L,K] gamma;                  // year-level effect coefficents
  vector[H] beta_hurdle;              // effect coefficents for hurdle
  array[R] matrix[K,T] beta_raw;      // standardized year-level deviations
  vector[S] st_raw;                   // standardized site effect deviations
  cholesky_factor_corr[K] L_omega;    // Cholesky transformed correlation matrix
  vector<lower=0>[K] tau;             // year-level effect coefficents error
  real<lower=0> tau_s;                // site random effect error
  real<lower=0> shape;                // Gamma shape parameter

}

transformed parameters {

  // year-level MVN random effects - noncentered parametrization
  array[R] matrix[K,T] beta;          // year_region-level effect coefficents
  for (r in 1:R){
    matrix[K,T] beta_error;
    beta_error = diag_pre_multiply(tau,L_omega)*beta_raw[r];
    beta[r] = (z[r] * gamma)' + beta_error;
  }
  
  // site-level MVN random effects - noncentered parametrization
  vector[S] st_eff = tau_s * st_raw;
  
  // hurdle parameters 
  vector[N] logit_hu = x_hurdle * beta_hurdle;
  vector[N] hu = inv_logit(logit_hu);
  
  // Mu vector
  vector[N] eta;

  for (n in 1:N){
    int t = yr[n];
    int s = st[n];
    int r = rg[s];
    eta[n] = x[n] * beta[r][,t] + st_eff[s];
  }
  vector[N] mu = exp(eta);
    
}

model {

  // Hyperpriors
  to_vector(gamma) ~ normal(0, 5);
  beta_hurdle ~ normal(0, 5);

  for (r in 1:R)
    to_vector(beta_raw[r]) ~ std_normal();

  st_raw ~ std_normal();

  L_omega ~ lkj_corr_cholesky(2);

  tau ~ normal(0, 1);
  tau_s ~ normal(0, 1);

  shape ~ exponential(1);
  // shape ~cauchy(0,.5);

  // Likelihood
  for (n in 1:N) {

    if (y[n] == 0) {

      target += bernoulli_lpmf(1 | hu[n]);

    } else {
      
      // // Mu vector
      // int t = yr[n];
      // int s = st[n];
      // int r = rg[s];
      // real eta = x[n] * beta[r][,t] + st_eff[s];
      // real mu = exp(eta);
      
      target += bernoulli_lpmf(0 | hu[n]);
      // target += gamma_lpdf(y[n] | shape, shape / mu);
      target += gamma_lpdf(y[n] | shape, shape / mu[n]);
    
    }
  }
}

generated quantities {

  vector[N] log_lik;
  vector[N] y_rep;
  array[N] real residuals;

  for (n in 1:N) {

    // Log likelihood
    if (y[n] == 0) {

      log_lik[n] = bernoulli_lpmf(1 | hu[n]);

    } else {

      log_lik[n] = bernoulli_lpmf(0 | hu[n]) +
        gamma_lpdf(y[n] | shape, shape / mu[n]);

    }

    // residuals
    residuals[n] = y[n]-mu[n];
    
    // Posterior predictive draw
    if (bernoulli_rng(hu[n]) == 1) {
      y_rep[n] = 0;
    } else {
      y_rep[n] = gamma_rng(
        shape,
        shape / mu[n]
      );
    }
  }
}

