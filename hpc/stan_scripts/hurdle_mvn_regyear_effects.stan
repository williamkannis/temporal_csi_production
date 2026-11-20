//  CHANGE MODEL SYMBOLS TO MATCH MANUSCRIPT

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
  vector[H] delta;              // effect coefficents for hurdle
  array[R] matrix[K,T] theta_raw;      // standardized year-level deviations
  vector[S] st_raw;                   // standardized site effect deviations
  cholesky_factor_corr[K] L_omega;    // Cholesky transformed correlation matrix
  vector<lower=0>[K] sigma;             // year-level effect coefficents error
  real<lower=0> tau;                // site random effect error
  // vector<lower=0>[R] tau;           // site random effect error
  real<lower=0> alpha;                // Gamma shape parameter

}

transformed parameters {

  // year-level MVN random effects - noncentered parametrization
  array[R] matrix[K,T] theta;          // year_region-level effect coefficents
  for (r in 1:R){
    matrix[K,T] theta_error;
    theta_error = diag_pre_multiply(sigma,L_omega)*theta_raw[r];
    theta[r] = (z[r] * gamma)' + theta_error;
  }
  
  // site-level random effects - noncentered parametrization
  vector[S] st_eff = tau * st_raw;
  // vector[S] st_eff = tau[rg] .* st_raw;
  
  // hurdle parameters 
  vector[N] logit_hu = x_hurdle * delta;
  vector[N] hu = inv_logit(logit_hu);
  
  // Mu vector
  vector[N] eta;

  for (n in 1:N){
    int t = yr[n];
    int s = st[n];
    int r = rg[s];
    eta[n] = x[n] * theta[r][,t] + st_eff[s];
  }
  vector[N] mu = exp(eta);
    
}

model {

  // Hyperpriors
  to_vector(gamma) ~ normal(0, 1);
  delta ~ normal(0, 5);

  for (r in 1:R)
    to_vector(theta_raw[r]) ~ std_normal();

  st_raw ~ std_normal();

  L_omega ~ lkj_corr_cholesky(2);

  sigma ~ normal(0, 1);
  tau ~ normal(0, 1);

  alpha ~ exponential(1);
  // alpha ~cauchy(0,.5);

  // Likelihood
  for (n in 1:N) {

    if (y[n] == 0) {

      target += bernoulli_lpmf(1 | hu[n]);

    } else {
      
      // // Mu vector
      // int t = yr[n];
      // int s = st[n];
      // int r = rg[s];
      // real eta = x[n] * theta[r][,t] + st_eff[s];
      // real mu = exp(eta);
      
      target += bernoulli_lpmf(0 | hu[n]);
      // target += gamma_lpdf(y[n] | alpha, alpha / mu);
      target += gamma_lpdf(y[n] | alpha, alpha / mu[n]);
    
    }
  }
}

generated quantities {

  vector[N] log_lik;
  vector[N] y_rep;
  array[N] real raw_residuals;
  array[N] real pearson_residual_hurdle;
  array[N] real pearson_residual_gamma;

  for (n in 1:N) {

    // Log likelihood
    if (y[n] == 0) {

      log_lik[n] = bernoulli_lpmf(1 | hu[n]);

    } else {

      log_lik[n] = bernoulli_lpmf(0 | hu[n]) +
        gamma_lpdf(y[n] | alpha, alpha / mu[n]);

    }

    // Residuals
    raw_residuals[n] = y[n]-mu[n];
    

    {
      real is_zero = (y[n] == 0);

      pearson_residual_hurdle[n] =
        (is_zero - hu[n]) /
        sqrt(hu[n] * (1 - hu[n]));
    }

    if (y[n] > 0) {

      pearson_residual_gamma[n] =
        (y[n] - mu[n]) /
        (mu[n] / sqrt(alpha));

    } else {

      pearson_residual_gamma[n] = not_a_number();
    }
    
    // Posterior predictive draw
    if (bernoulli_rng(hu[n]) == 1) {
      y_rep[n] = 0;
    } else {
      y_rep[n] = gamma_rng(
        alpha,
        alpha / mu[n]
      );
    }
  }
}

