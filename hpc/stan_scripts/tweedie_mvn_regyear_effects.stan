data {
  int M;
  int<lower=0> N;                     // number of data items
  int<lower=1> T;                     // number of years
  int<lower=1> S;                     // number of sites
  int<lower=1> R;                     // number of regions
  int<lower=1> K;                     // number of predictors
  int<lower=1> L;                     // number of 2nd-level predictors
  
  array[N] int<lower=1, upper=T> yr;  // year index
  array[N] int<lower=1, upper=S> st;  // site 
  array[S] int<lower=1, upper=R> rg;  // region-site bridge
  array[N] real y;                    // outcome array
  array[N] row_vector[K] x;           // period-level predictor data
  array[R] matrix[T, L] z;            // year-level predictor data, by region
}

parameters {
  matrix[L,K] gamma;                  // year-level effect coefficents
  array[R] matrix[K,T] theta_raw;      // standardized year-level deviations
  vector[S] st_raw;                   // standardized site effect deviations
  cholesky_factor_corr[K] L_omega;    // Cholesky transformed correlation matrix
  vector<lower=0>[K] sigma;             // year-level effect coefficents error
  real<lower=0> tau;                // site random effect error
  // vector<lower=0>[R] tau;           // site random effect error
  real<lower=0> phi;
  real<lower=1, upper=2> power;
  // real<lower=1.01, upper=1.99> power;
}

transformed parameters {
  array[R] matrix[K,T] theta;          // year_region-level effect coefficents
  vector[S] st_eff;                   // site random effect

  // year-level MVN random effects - noncentered parametrization
  for (r in 1:R){
    matrix[K,T] theta_error;
    theta_error = diag_pre_multiply(sigma,L_omega)*theta_raw[r];
    theta[r] = (z[r] * gamma)' + theta_error;
  }
  
  // site-level MVN random effects - noncentered parametrization
  st_eff = tau * st_raw;
  // st_eff = tau[rg] .* st_raw;
  
// Mu vector
    vector[N] eta;    
    vector[N] mu;
    
    for (n in 1:N){
      int t = yr[n];
      int s = st[n];
      int r = rg[s];
      eta[n] = x[n] * theta[r][,t] + st_eff[s];
    }
    mu = exp(eta);
    
  // Tweetie parameters 
  vector[N] lambda;
  vector[N] beta_gamma;
  real alpha;
  vector[M] log_m_factorial;

  alpha = (2 - power) / (power - 1);

  lambda =
      exp((2 - power) * eta)
      / phi
      / (2 - power);

  beta_gamma =
      exp((1 - power) * eta)
      / phi
      / (power - 1);
      
  for (m in 1:M)
    log_m_factorial[m] = lgamma(m + 1);


}

model {

  // Hyperpriors
  to_vector(gamma) ~ normal(0, 1);

  for (r in 1:R)
    to_vector(theta_raw[r]) ~ std_normal();

  st_raw ~ std_normal();

  L_omega ~ lkj_corr_cholesky(2);

  sigma ~ normal(0, 1);
  tau ~ normal(0, 1);

  // Tweedie dispersion
  // phi ~ cauchy(0, 5);
  phi ~ lognormal(0,1);

  // Likelihood
  for (n in 1:N) {

    if (y[n] == 0) {

      target += -lambda[n];

    } else {

      vector[M] lp;

      real common;
      real log_term;

      common =
          -lambda[n]
          - beta_gamma[n] * y[n]
          - log(y[n]);

      log_term =
          log(lambda[n])
          + alpha * log(beta_gamma[n])
          + alpha * log(y[n]);

      for (m in 1:M) {

        lp[m] =
            common
            + m * log_term
            - lgamma(m * alpha)
            - lgamma(m + 1);
      }

      target += log_sum_exp(lp);
    }
  }
}

generated quantities {

  vector[N] log_lik;
  vector[N] y_rep;
  vector[N] p_zero;
  array[N] real raw_residuals;
  array[N] real pearson_residual;

  for (n in 1:N) {

    // Log likelihood
    if (y[n] == 0) {

      log_lik[n] = -lambda[n];

    } else {

      vector[M] lp;

      real common;
      real log_term;

      common =
          -lambda[n]
          - beta_gamma[n] * y[n]
          - log(y[n]);

      log_term =
          log(lambda[n])
          + alpha * log(beta_gamma[n])
          + alpha * log(y[n]);

      for (m in 1:M) {

        lp[m] =
            common
            + m * log_term
            - lgamma(m * alpha)
            - lgamma(m + 1);
      }

      log_lik[n] = log_sum_exp(lp);
    }

    // residuals
    raw_residuals[n] = y[n]-mu[n];
    pearson_residual[n] = raw_residuals[n] / sqrt(phi * pow(mu[n], power));

    // probability of zero
    p_zero[n] = exp(-lambda[n]);

    // Posterior predictive draw
    {
      int m_rep;

      m_rep = poisson_rng(lambda[n]);

      if (m_rep == 0) {

        y_rep[n] = 0;

      } else {

        y_rep[n] =
          gamma_rng(
            m_rep * alpha,
            beta_gamma[n]
          );
      }
    }
  }
}
