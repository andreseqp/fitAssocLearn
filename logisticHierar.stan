data {
  int<lower=1> D; # # of predictor 
  int<lower=0> N; # # of outcomes
  int<lower=1> L; # # of levels
  array[N] int<lower=0, upper=1> y; # outcome
  array[N] int<lower=1, upper=L> ll; # levels
  array[N] row_vector[D] x; # predictors
}
parameters {
  array[D] real mu; # 
  array[D] real<lower=0> sigma;
  array[L] vector[D] beta;
}
model {
  for (d in 1:D) {
    mu[d] ~ normal(0, 100);
    for (l in 1:L) {
      beta[l, d] ~ normal(mu[d], sigma[d]);
    }
  }
  for (n in 1:N) {
    y[n] ~ bernoulli(inv_logit(x[n] * beta[ll[n]]));
  }
}
