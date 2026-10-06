data {
  int<lower=1> N;
  vector[N] y;
  vector<lower=0>[N] se;
  vector[N] tdec;
}
parameters {
  real beta0;
  real beta1;
  real<lower=0> tau;
}
model {
  beta0 ~ normal(0, 1);
  beta1 ~ normal(0, 0.25);
  tau   ~ normal(0, 0.5);

  y ~ normal(
    beta0 + beta1 * tdec,
    sqrt(square(se) + square(tau))
  );
}
generated quantities {
  real PF_2000;
  real PF_2010;
  real PF_2020;
  real PF_2025;

  PF_2000 = 1 - exp(beta0);
  PF_2010 = 1 - exp(beta0 + beta1);
  PF_2020 = 1 - exp(beta0 + 2 * beta1);
  PF_2025 = 1 - exp(beta0 + 2.5 * beta1);
}
