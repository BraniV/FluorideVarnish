data {
  int<lower=1> N;
  vector[N] y;
  vector<lower=0>[N] se;
  vector[N] tdec;
  vector[N] B;
}
parameters {
  real beta0;
  real beta1;
  real beta2;
  real<lower=0> tau;
}
model {
  beta0 ~ normal(0, 1);
  beta1 ~ normal(0, 0.25);
  beta2 ~ normal(0, 0.25);
  tau   ~ normal(0, 0.5);

  y ~ normal(
    beta0 + beta1 * tdec + beta2 * B,
    sqrt(square(se) + square(tau))
  );
}
