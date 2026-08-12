data {
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
  array[R] matrix[T, L] z;              // year-level predictor data, by region
}

parameters {
  matrix[L,K] gamma;                  // year-level effect coefficents
  array[R] matrix[K,T] beta_raw;      // standardized year-level deviations
  matrix[K,R] rg_raw;                 // standardized region effect deviations
  matrix[K,S] st_raw;                 // standardized site effect deviations
  matrix[K,T] yr_raw;                 // standardized year effect deviations
  cholesky_factor_corr[K] L_omega;    // Cholesky transformed correlation matrix
  cholesky_factor_corr[K] L_omega_r;  // Cholesky transformed correlation matrix
  cholesky_factor_corr[K] L_omega_s;  // Cholesky transformed correlation matrix
  cholesky_factor_corr[K] L_omega_y;  // Cholesky transformed correlation matrix
  vector<lower=0>[K] tau;             // year-level effect coefficents error
  vector<lower=0>[K] tau_r;           // region random effect error
  vector<lower=0>[K] tau_s;           // site random effect error
  vector<lower=0>[K] tau_y;           // year random effect error
  real<lower=0> sigma;                // error scale
}

transformed parameters {
  array[R] matrix[K,T] beta;          // year_region-level effect coefficents
  matrix[K,R] rg_eff;                 // region random effect
  matrix[K,S] st_eff;                 // site random effect
  matrix[K,T] yr_eff;                 // year random effect

  // MVN random effects - noncentered parametrization

  // region-year-level
  for (r in 1:R){
    matrix[K,T] beta_error;
    beta_error = diag_pre_multiply(tau,L_omega)*beta_raw[r];
    beta[r] = (z[r] * gamma)' + beta_error;
  }
  
  // Regions-level
  rg_eff = diag_pre_multiply(tau_r, L_omega_r) * rg_raw;
  
  // site-level 
  for (s in 1:S) 
  st_eff[,s] = diag_pre_multiply(tau_s, L_omega_s) * st_raw[,s] + rg_eff[,rg[s]];
  
  // year-level
  yr_eff = diag_pre_multiply(tau_y, L_omega_y) * yr_raw;
}

model {
  
  // hyperpriors
  to_vector(gamma) ~ normal(0,10);
  for(r in 1:R) to_vector(beta_raw[r]) ~ std_normal();
  to_vector(rg_raw)  ~ std_normal();
  to_vector(st_raw)  ~ std_normal();
  to_vector(yr_raw)  ~ std_normal();
  L_omega ~ lkj_corr_cholesky(2);
  L_omega_r ~ lkj_corr_cholesky(2);
  L_omega_s ~ lkj_corr_cholesky(2);
  L_omega_y ~ lkj_corr_cholesky(2);
  tau ~ cauchy(0,1);
  tau_r ~ cauchy(0,1);
  tau_s ~ cauchy(0,1);
  tau_y ~ cauchy(0,1);

  
  // 1st-level priors
  sigma ~ cauchy(0,1);
  
  // Likelihood
  vector[N] mu;                       
  for (n in 1:N){
    
    int t = yr[n];
    int s = st[n];
    int r = rg[s];
    mu[n] = x[n] * (beta[r][,t] + st_eff[,s]+yr_eff[,t]);
  } 
  y ~ normal(mu,sigma);
}

generated quantities{
  
   // retrieve correlation matrix from cholesky matrix
  corr_matrix[K] cor_mat;            // year correlation matrix
  corr_matrix[K] cor_mat_s;          // site correlation matrix
  

  cor_mat = multiply_lower_tri_self_transpose(L_omega);
  cor_mat_s = multiply_lower_tri_self_transpose(L_omega_s);
  
  // Model fit and posterior checks
  array[N] real mu;                   // Predicted values
  array[N] real residuals;            // model residuals
  array[N] real y_rep;                // generated data
  vector[N] log_lik;                  // log-likelihood vector
  int<lower = 0, upper = 1> mean_gt;  // mean "pvalue""
  int<lower = 0, upper = 1> sd_gt;    // sd "pvalue"
  int<lower = 0, upper = 1> max_gt;   // max "pvalue""
  int<lower = 0, upper = 1> min_gt;   // min "pvalue"
  int neg_rep = 0;
  
  for (n in 1:N){
  
  // Predicted values
    // real mu;
    int t = yr[n];
    int s = st[n];
    int r = rg[s];
    mu[n] = x[n] * (beta[r][,t] + st_eff[,s]+yr_eff[,t]);
    
  // residuals
    residuals[n] = y[n]-mu[n];
    
  // Generated data
    y_rep[n] = normal_rng(mu[n], sigma);
    neg_rep += (y_rep[n] < 0);
    
  // Log pointwise predictive density
    log_lik[n] = normal_lpdf(y[n]| mu[n],sigma);
  }
  
  // Posteriror preditive check
  mean_gt = mean(y_rep) > mean(y);
  sd_gt = sd(y_rep) > sd(y);
  max_gt = max(y_rep) > max(y);
  min_gt = min(y_rep) > min(y);
} 
