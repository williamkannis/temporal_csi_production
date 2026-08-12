data {
  int<lower=0> N;                     // number of data items
  int<lower=1> P;                     // random effect combined group size
  int<lower=1> L;                     // random effect group1 size
  int<lower=1> O;                     // random effect group2 size
  int<lower=1> K;                     // number of predictors
  int<lower=1> J;                     // number of 2nd-level predictors
  array[N] real y;                    // outcome array
  array[N] int<lower=1, upper=P> pp;  // combined group array
  array[P] int<lower=1, upper=L> ll;  // group1 array
  array[P] int<lower=1, upper=O> oo;  // group1 array
  array[N] row_vector[K] x;           // predictor array
  array[P] row_vector[J] z;           // 2nd-level predictor array - group1
}
parameters {
  matrix[J,K] gamma;                  // 2nd-level effect coefficents
  matrix[K,P] beta_raw;               // standardized group1-level deviations
  vector<lower=0>[K] tau;             // 2nd-level effect coefficents error
  cholesky_factor_corr[K] L_omega;    // Cholesky transformed correlation matrix
  matrix[K,L] l_eff;                  // group1 random effect
  matrix[K,O] o_eff;                  // group2 random effect
  vector<lower=0>[K] tau1;            // group1 random effect error
  vector<lower=0>[K] tau2;            // group2 random effect error
  real<lower=0> sigma;                // error scale
}

transformed parameters {
  matrix[K,P] beta_error;             // group1-level effect coefficents error
  array[P] vector[K] beta;            // combined group-level effect coefficents

  // Correlated Group-level error
  beta_error = diag_pre_multiply(tau,L_omega)*beta_raw;
  
  // Group-level MVN random effects - noncentered parametrization
  for (p in 1:P) 
  beta[p] = (z[p] * gamma)' + beta_error[,p] + l_eff[,ll[p]] + o_eff[,oo[p]];

}

model {
  
  // hyperpriors
  to_vector(gamma) ~ normal(0,10);
  to_vector(beta_raw) ~ std_normal();
  L_omega ~ lkj_corr_cholesky(2);
  for (k in 1:K){
    to_vector(l_eff[k,]) ~ normal(0,tau1[k]);
    to_vector(o_eff[k,]) ~ normal(0,tau2[k]);
  }
  tau ~ cauchy(0,1);
  tau1 ~ cauchy(0,1);
  tau2 ~ cauchy(0,1);
  
  // 1st-level priors
  sigma ~ cauchy(0,1);
  
  // Likelihood
  vector[N] mu;                       // group-specific means
  
  for (n in 1:N) mu[n] = x[n] * beta[pp[n]];
  y ~ normal(mu,sigma);
}
generated quantities{
  corr_matrix[K] cor_mat1;            // correlation matrix

  // retrieve correlation matrix from cholesky matrix
  cor_mat1 = multiply_lower_tri_self_transpose(L_omega);
  
} 
  
  