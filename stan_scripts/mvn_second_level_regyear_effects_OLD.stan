data {
  int<lower=0> N;                     // number of data items
  int<lower=1> P;                     // random effect combined group size
  int<lower=1> L;                     // random effect group1 size
  int<lower=1> O;                     // random effect group2 size
  int<lower=1> K;                     // number of predictors
  int<lower=1> J;                     // number of 2nd-level predictors
  array[N] real y;                    // outcome array
  array[N] int<lower=1, upper=P> pp;  // combined group array
  array[N] int<lower=1, upper=L> ll;  // group1 array
  array[N] int<lower=1, upper=O> oo;  // group1 array
  array[N] row_vector[K] x;           // predictor array
  array[P] row_vector[J] z;           // 2nd-level predictor array - group1
}
parameters {
  matrix[J,K] gamma;                  // 2nd-level effect coefficents
  matrix[K,P] beta_raw;               // standardized group1-level deviations
  vector<lower=0>[K] tau;             // 2nd-level effect coefficents error
  cholesky_factor_corr[K] L_omega;    // Cholesky transformed correlation matrix
  vector[L] ll_raw;                   // standardized group1-level deviations
  vector[O] oo_raw;                   // standardized group2-level deviations
  real<lower=0> tau1;                 // group1 random effect error
  real<lower=0> tau2;                 // group2 random effect error
  real<lower=0> sigma;                // error scale
}

transformed parameters {
  matrix[K,P] beta_error;             // group1-level effect coefficents error
  array[P] vector[K] beta;            // combined group-level effect coefficents
  vector[L] ll_eff;                   // group1 random effect
  vector[O] oo_eff;                   // group2 random effect

  // Correlated Group-level error
  beta_error = diag_pre_multiply(tau,L_omega)*beta_raw;
  
  // Group-level MVN random effects - noncentered parametrization
  for (p in 1:P) beta[p] = (z[p] * gamma)' + beta_error[,p];
  ll_eff = tau1 * ll_raw;
  oo_eff = tau2 * oo_raw;
}

model {
  
  // hyperpriors
  to_vector(gamma) ~ normal(0,10);
  to_vector(beta_raw) ~ std_normal();
  ll_raw  ~ std_normal();
  oo_raw  ~ std_normal();
  L_omega ~ lkj_corr_cholesky(2);
  tau ~ cauchy(0,1);
  tau1 ~ cauchy(0,1);
  tau2 ~ cauchy(0,1);
  
  // 1st-level priors
  sigma ~ cauchy(0,1);
  
  // Likelihood
  vector[N] mu;                       // group-specific means
  
  for (n in 1:N) mu[n] = x[n] * beta[pp[n]] + ll_eff[ll[n]] + oo_eff[oo[n]];
  y ~ normal(mu,sigma);
}
generated quantities{
  
   // retrieve correlation matrix from cholesky matrix
  corr_matrix[K] cor_mat1;            // correlation matrix

  cor_mat1 = multiply_lower_tri_self_transpose(L_omega);
  
  // Model fit and posterior checks
  array[N] real mu;                   // Predicted values
  array[N] real residuals;            // model residuals
  array[N] real y_rep;                // generated data
  vector[N] log_lik;                  // log-likelihood vector
  int<lower = 0, upper = 1> mean_gt;  // mean "pvalue""
  int<lower = 0, upper = 1> sd_gt;    // sd "pvalue"
  int<lower = 0, upper = 1> max_gt;  // max "pvalue""
  int<lower = 0, upper = 1> min_gt;    // min "pvalue"
  int neg_rep = 0;
  
  for (n in 1:N){
  
  // Predicted values
    // real mu;
    mu[n] = x[n] * beta[pp[n]] + ll_eff[ll[n]] + oo_eff[oo[n]];
    
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