data {
  int<lower=1> N;
  vector[N] xT;
  vector[N] xC;
  vector<lower=0>[N] vT;
  vector<lower=0>[N] vC;
  vector[N] tdec;
  vector[N] B;
  real xCbar;
  real<lower=0> s_logC;
}
parameters {
  real beta0;
  real beta1;
  real beta2;
  real beta3;
  real<lower=0> tau;
}
model {
  vector[N] C;
  vector[N] mu_T;
  vector[N] sd_T;

  beta0 ~ normal(0, 1);
  beta1 ~ normal(0, 0.25);
  beta2 ~ normal(0, 0.25);
  beta3 ~ normal(0, 0.25);
  tau   ~ normal(0, 0.5);

  C = (xC - xCbar) / s_logC;

  for (i in 1:N) {
    mu_T[i] = xC[i] + beta0 + beta1 * tdec[i] + beta2 * B[i] + beta3 * C[i];

    // Integrated measurement-error likelihood obtained with a locally flat
    // prior for the latent true log control mean.
    sd_T[i] = sqrt(
      vT[i] + square(tau)
      + square(1 + beta3 / s_logC) * vC[i]
    );
  }

  xT ~ normal(mu_T, sd_T);
}
