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
  cholesky_factor_corr[K] L_omega1;   // Cholesky transformed correlation matrix
  real<lower=0> tau2;                 // group2 effect coefficents error
  matrix[K,L] beta_raw1;              // standardized group1-level deviations
  matrix[K,O] beta_raw_s;             // standardized group1-level deviations
  // array[O] real beta_raw2;            // standardized group2-level deviations
  real<lower=0> sigma;                // error scale
}

transformed parameters {
  // matrix[K,L] beta1_error;            // group1-level effect coefficents error
  array[L] vector[K] beta1;           // group1-level effect coefficents
  // array[O] real beta2;                // group2-level effect coefficents
  array[O] vector[K] site_eff;
  
  // // Correlated Group-level error
  // beta1_error = diag_pre_multiply(tau1,L_omega1)*beta_raw1;
  // 
  // // Group-level MVN random effects - noncentered parametrization
  // for (l in 1:L) beta1[l] = (z1[l] * gamma1)' + beta1_error[,l];
  // for (o in 1:O) beta2[o] = tau2 .* beta_raw2[o];
  
  // Year-level MVN random effects - noncentered parametrization
  for(l in 1:L){
    vector[K] beta1_error = diag_pre_multiply(tau1,L_omega1)*beta_raw1[,l];
    beta1[l] = (z1[l] * gamma1)' + beta1_error;
  }
  
  // site-level MVN random effects - noncentered parametrization
  for(o in 1:O) st_eff[o] = diag_pre_multiply(tau_s, L_omega_s) * st_raw[,o];
  
}

model {
  
  // hyperpriors
  to_vector(beta_raw1) ~ std_normal();
  to_vector(beta_raw_s) ~ std_normal();
  // for (o in 1:O) beta_raw2[o] ~ std_normal();
  to_vector(gamma1) ~ normal(0,10);
  tau1 ~ cauchy(0,1);
  tau_s ~ cauchy(0,1);
  L_omega1 ~ lkj_corr_cholesky(2);
  L_omega_s ~ lkj_corr_cholesky(2);
  // tau2 ~ cauchy(0,1);
  
  // 1st-level priors
  sigma ~ cauchy(0,1);
  
  // Likelihood
  vector[N] mu;                       // group-specific means
  
  // for (n in 1:N) mu[n] = x[n] * beta1[ll[n]] + beta2[oo[n]];
  for (n in 1:N) {
    int y = ll[n];
    int s = oo[n];
    mu = x[n]*(beta1[y]+st_eff[s])
  }
  y ~ normal(mu,sigma);
}
generated quantities{
  corr_matrix[K] cor_mat1;            // correlation matrix

  // retrieve correlation matrix from cholesky matrix
  cor_mat1 = multiply_lower_tri_self_transpose(L_omega1);
  
} 
  
  