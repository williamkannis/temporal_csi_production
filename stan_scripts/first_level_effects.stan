data {
  int<lower=0> N;   // number of data items
  int<lower=0> K;   // number of predictors
  matrix[N, K] x;   // predictor matrix
  vector[N] y;      // outcome vector
}
parameters {
  vector[K] beta;       // coefficients for predictors
  real<lower=0> sigma;  // error scale
}
model {
  
  beta ~ normal(0,10);
  sigma ~ cauchy(0,1);
  
  
  y ~ normal(x * beta, sigma);  // likelihood
}
