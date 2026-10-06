data {
  int<lower=1> N;

  vector[N] mean_C_obs;
  vector[N] mean_T_obs;

  vector<lower=0>[N] se_C;
  vector<lower=0>[N] se_T;

  // 0 = Gaussian random effects
  // 1 = standardized Student-t_4 random effects
  int<lower=0, upper=1> use_t4;

  // Index of Milsom, used only for one compact predictive-tail diagnostic.
  int<lower=1, upper=N> milsom_index;
}

parameters {
  // Study-specific baseline log means:
  // m_Ci = exp(eta_Ci) > 0.
  vector[N] eta_C;

  // Pooled log treatment-to-control mean ratio.
  real mu;

  // Between-study SD.  Prior below is Half-Normal(0,0.5^2)
  // because the parameter is constrained to be nonnegative.
  real<lower=0> tau;

  // Noncentered standardized random effects.
  vector[N] z;
}

transformed parameters {
  vector[N] theta;
  vector[N] mean_C;
  vector[N] mean_T;

  theta  = mu + tau * z;
  mean_C = exp(eta_C);
  mean_T = mean_C .* exp(theta);
}

model {
  // Priors used in the manuscript.
  mu    ~ normal(0, 1);
  tau   ~ normal(0, 0.5);
  eta_C ~ normal(0, 2.5);

  // Random-effects distribution.
  //
  // For nu=4, Student-t variance is scale^2 * nu/(nu-2).
  // Choosing scale=sqrt((nu-2)/nu)=sqrt(1/2) makes Var(z)=1,
  // so tau keeps the interpretation of a between-study SD.
  if (use_t4 == 1) {
    z ~ student_t(4, 0, sqrt(0.5));
  } else {
    z ~ std_normal();
  }

  // Aggregate arm-level likelihood.
  //
  // The reported arm SDs are treated as fixed and the sampling
  // distributions of the reported means are approximated as Gaussian.
  mean_C_obs ~ normal(mean_C, se_C);
  mean_T_obs ~ normal(mean_T, se_T);
}

generated quantities {
  real PF_pool;

  // Posterior predictive draw for a NEW TRUE study effect.
  real z_new;
  real theta_new;
  real PF_new;

  int<lower=0, upper=1> benefit_new;
  int<lower=0, upper=1> benefit20_new;

  // Study-specific true prevented fractions.
  vector[N] PF_study;

  // Posterior predictive replicated arm means.
  vector[N] mean_C_rep;
  vector[N] mean_T_rep;

  // Pointwise log likelihood for possible later LOO work.
  vector[2 * N] log_lik;

  // Posterior predictive discrepancy measures based on standardized
  // arm-level residuals.
  real T_ss_obs;
  real T_ss_rep;
  real T_max_obs;
  real T_max_rep;

  int<lower=0, upper=1> ppc_ss;
  int<lower=0, upper=1> ppc_max;

  // How unusual is Milsom's latent standardized random effect relative
  // to a fresh draw from the fitted random-effects population?
  int<lower=0, upper=1> ppc_milsom_re_extreme;

  PF_pool = 1 - exp(mu);

  if (use_t4 == 1) {
    z_new = student_t_rng(4, 0, sqrt(0.5));
  } else {
    z_new = normal_rng(0, 1);
  }

  theta_new = mu + tau * z_new;
  PF_new = 1 - exp(theta_new);

  benefit_new   = PF_new > 0;
  benefit20_new = PF_new > 0.20;

  T_ss_obs  = 0;
  T_ss_rep  = 0;
  T_max_obs = 0;
  T_max_rep = 0;

  for (i in 1:N) {
    real zc_obs;
    real zt_obs;
    real zc_rep;
    real zt_rep;

    PF_study[i] = 1 - exp(theta[i]);

    mean_C_rep[i] = normal_rng(mean_C[i], se_C[i]);
    mean_T_rep[i] = normal_rng(mean_T[i], se_T[i]);

    log_lik[i] =
      normal_lpdf(mean_C_obs[i] | mean_C[i], se_C[i]);

    log_lik[N + i] =
      normal_lpdf(mean_T_obs[i] | mean_T[i], se_T[i]);

    zc_obs = (mean_C_obs[i] - mean_C[i]) / se_C[i];
    zt_obs = (mean_T_obs[i] - mean_T[i]) / se_T[i];

    zc_rep = (mean_C_rep[i] - mean_C[i]) / se_C[i];
    zt_rep = (mean_T_rep[i] - mean_T[i]) / se_T[i];

    T_ss_obs += square(zc_obs) + square(zt_obs);
    T_ss_rep += square(zc_rep) + square(zt_rep);

    T_max_obs = fmax(T_max_obs, fabs(zc_obs));
    T_max_obs = fmax(T_max_obs, fabs(zt_obs));

    T_max_rep = fmax(T_max_rep, fabs(zc_rep));
    T_max_rep = fmax(T_max_rep, fabs(zt_rep));
  }

  ppc_ss  = T_ss_rep  >= T_ss_obs;
  ppc_max = T_max_rep >= T_max_obs;

  ppc_milsom_re_extreme =
    fabs(z_new) >= fabs(z[milsom_index]);
}
