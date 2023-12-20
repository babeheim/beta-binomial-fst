data {
  int<lower=1> m;
  array[m] int<lower=1> n;
  array[m] int<lower=0, upper=n> x;
}
parameters {
  real<lower=0, upper=1> mu;
  real<lower=0, upper=1> F;
}
transformed parameters {
  real theta = (1 - F) / F;
  real<lower=0> a = mu * theta;
  real<lower=0> b = (1 - mu) * theta;
}
model {
  F ~ beta(2, 2);
  mu ~ beta(2, 2);
  target += beta_binomial_lpmf(x | n, a, b);
}
