# ==============================================================================
# 03_historical_armlevel_Stan.R
#
# Production Stan analysis for the two nonlinear Bayesian historical models:
#
#   (A) aggregate arm-level model with Gaussian random effects;
#   (B) aggregate arm-level model with standardized Student-t_4 random effects.
#
# This script:
#   * compiles historical_armlevel.stan using CmdStanR;
#   * runs 4 independent NUTS chains for each model;
#   * uses a noncentered random-effects parameterization;
#   * reports mu, tau, pooled PF, predictive PF and predictive probabilities;
#   * reports MCSE, ESS_bulk, ESS_tail and R-hat;
#   * reports divergences, maximum-treedepth hits and E-BFMI;
#   * generates posterior predictive checks at the arm-mean level;
#   * summarizes Milsom's study-specific posterior effect;
#   * generates trace plots and posterior-predictive figures;
#   * writes LaTeX-ready text for the main manuscript and Supplement.
#
# MODEL MATCH TO CURRENT MANUSCRIPT
# -------------------------
# The model is exactly the aggregate-data model described in manuscript:
#
#   Xbar_Ci | m_i ~ N(m_i, s_Ci^2/n_Ci)
#   Xbar_Ti | m_i,theta_i ~ N(m_i exp(theta_i), s_Ti^2/n_Ti)
#   log(m_i)=eta_i,  eta_i ~ N(0,2.5^2)
#   mu ~ N(0,1^2)
#   tau ~ Half-Normal(0,0.5^2)
#
# Gaussian RE:
#   theta_i = mu + tau z_i, z_i ~ N(0,1)
#
# Robust RE:
#   theta_i = mu + tau z_i,
#   z_i ~ t_4(0, sqrt((4-2)/4)),
# so Var(z_i)=1 and tau remains the between-study SD.
#
# IMPORTANT
# ---------
# The observed arm SDs are treated as fixed.  This model is an aggregate-data
# approximation and does not reconstruct zero inflation, individual-level
# correlation, or unreported cluster design effects.
#
# SOFTWARE
# --------
# Required R packages:
#   cmdstanr
#   posterior
#
# CmdStan itself must also be installed.
#
# Typical one-time installation:
#
# install.packages(
#   "cmdstanr",
#   repos = c("https://stan-dev.r-universe.dev", getOption("repos"))
# )
# cmdstanr::install_cmdstan()
#
# Then run this script from the directory containing:
#   historical9_stan_armlevel_production.R
#   historical_armlevel.stan
#
# ==============================================================================


# ------------------------------------------------------------------------------
# 0. Packages and CmdStan installation check
# ------------------------------------------------------------------------------

if (!requireNamespace("cmdstanr", quietly = TRUE)) {
  stop(
    "Package 'cmdstanr' is not installed.\n",
    "Install it first, then rerun this script."
  )
}

if (!requireNamespace("posterior", quietly = TRUE)) {
  stop("Package 'posterior' is not installed. Install it first.")
}

library(cmdstanr)
library(posterior)

cmdstan_ver <- tryCatch(
  cmdstanr::cmdstan_version(error_on_NA = FALSE),
  error = function(e) NA
)

if (length(cmdstan_ver) == 1 && is.na(cmdstan_ver)) {
  stop(
    "CmdStan is not installed or CmdStanR cannot find it.\n",
    "Run cmdstanr::install_cmdstan() once, then rerun this script."
  )
}

cat("CmdStan version:", as.character(cmdstan_ver), "\n")
cat("CmdStanR version:", as.character(packageVersion("cmdstanr")), "\n")
cat("posterior version:", as.character(packageVersion("posterior")), "\n")


# ------------------------------------------------------------------------------
# 1. Locate the script directory
# ------------------------------------------------------------------------------

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", args, value = TRUE)

  if (length(hit) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", hit[1]))))
  }

  normalizePath(getwd())
}

script_dir <- get_script_dir()

stan_file <- file.path(script_dir, "historical_armlevel.stan")

if (!file.exists(stan_file)) {
  stop("Cannot find Stan file: ", stan_file)
}

outdir <- file.path(script_dir, "output_03_historical_Stan")

if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}


# ------------------------------------------------------------------------------
# 2. Frozen nine-study arm-level data
# ------------------------------------------------------------------------------

dat <- data.frame(
  study = c(
    "Koch",
    "Modeer",
    "Clark",
    "Tewari",
    "Bravo",
    "Skold",
    "Arruda",
    "Tagliaferro",
    "Milsom"
  ),

  mean_T = c(0.90, 2.50, 2.43, 0.55, 1.48, 0.79, 4.61, 0.33, 0.66),
  sd_T   = c(3.80, 3.10, 3.09, 4.48, 1.53, 1.67, 4.39, 1.04, 0.73),
  n_T    = c(60,   87,   246,  311,  98,   190,  113,  91,   94),

  mean_C = c(4.00, 3.70, 3.11, 2.16, 2.58, 1.85, 7.72, 0.57, 0.63),
  sd_C   = c(3.75, 3.90, 3.54, 4.12, 1.89, 2.89, 5.84, 1.39, 0.66),
  n_C    = c(61,   107,  234,  307,  116,  181,  97,   86,   95)
)

stopifnot(
  all(dat$mean_T > 0),
  all(dat$mean_C > 0),
  all(dat$sd_T > 0),
  all(dat$sd_C > 0),
  all(dat$n_T > 1),
  all(dat$n_C > 1)
)

dat$se_T <- dat$sd_T / sqrt(dat$n_T)
dat$se_C <- dat$sd_C / sqrt(dat$n_C)

write.csv(
  dat,
  file.path(outdir, "stan_armlevel_input_data.csv"),
  row.names = FALSE
)

milsom_index <- which(dat$study == "Milsom")
stopifnot(length(milsom_index) == 1)


# ------------------------------------------------------------------------------
# 3. Compile the Stan model
# ------------------------------------------------------------------------------

cat("\nCompiling Stan model...\n")

mod <- cmdstanr::cmdstan_model(
  stan_file,
  force_recompile = FALSE
)


# ------------------------------------------------------------------------------
# 4. Sampling settings
# ------------------------------------------------------------------------------

# Four chains are used because rank-normalized split R-hat and ESS diagnostics
# are meaningful only with multiple chains.
#
# This is a small model, so we use deliberately conservative HMC settings.
# If all diagnostics are excellent, these settings are more than adequate.
#
# Total retained draws/model = 4 * 4000 = 16000.

chains          <- 4
iter_warmup     <- 2000
iter_sampling   <- 4000
adapt_delta     <- 0.995
max_treedepth   <- 15

available_cores <- parallel::detectCores(logical = TRUE)

parallel_chains <- if (is.na(available_cores)) {
  4
} else {
  min(chains, max(1, available_cores))
}

cat("\nSampling settings:\n")
cat("  chains          =", chains, "\n")
cat("  iter_warmup     =", iter_warmup, "\n")
cat("  iter_sampling   =", iter_sampling, "\n")
cat("  adapt_delta     =", adapt_delta, "\n")
cat("  max_treedepth   =", max_treedepth, "\n")
cat("  parallel_chains =", parallel_chains, "\n")


# ------------------------------------------------------------------------------
# 5. Stan data list
# ------------------------------------------------------------------------------

base_stan_data <- list(
  N = nrow(dat),

  mean_C_obs = dat$mean_C,
  mean_T_obs = dat$mean_T,

  se_C = dat$se_C,
  se_T = dat$se_T,

  milsom_index = milsom_index
)


# ------------------------------------------------------------------------------
# 6. Run Gaussian random-effects arm-level model
# ------------------------------------------------------------------------------

cat("\n============================================================\n")
cat("RUNNING GAUSSIAN RANDOM-EFFECTS ARM-LEVEL MODEL\n")
cat("============================================================\n\n")

normal_dir <- file.path(outdir, "cmdstan_normal")

if (!dir.exists(normal_dir)) {
  dir.create(normal_dir, recursive = TRUE)
}

data_normal <- c(
  base_stan_data,
  list(use_t4 = 0L)
)

fit_normal <- mod$sample(
  data = data_normal,

  seed = 20260927,

  chains = chains,
  parallel_chains = parallel_chains,

  iter_warmup = iter_warmup,
  iter_sampling = iter_sampling,

  adapt_delta = adapt_delta,
  max_treedepth = max_treedepth,

  refresh = 500,

  output_dir = normal_dir,
  output_basename = "armlevel_normal"
)


# ------------------------------------------------------------------------------
# 7. Run standardized t4 random-effects arm-level model
# ------------------------------------------------------------------------------

cat("\n============================================================\n")
cat("RUNNING STANDARDIZED t4 RANDOM-EFFECTS ARM-LEVEL MODEL\n")
cat("============================================================\n\n")

t4_dir <- file.path(outdir, "cmdstan_t4")

if (!dir.exists(t4_dir)) {
  dir.create(t4_dir, recursive = TRUE)
}

data_t4 <- c(
  base_stan_data,
  list(use_t4 = 1L)
)

fit_t4 <- mod$sample(
  data = data_t4,

  seed = 20260928,

  chains = chains,
  parallel_chains = parallel_chains,

  iter_warmup = iter_warmup,
  iter_sampling = iter_sampling,

  adapt_delta = adapt_delta,
  max_treedepth = max_treedepth,

  refresh = 500,

  output_dir = t4_dir,
  output_basename = "armlevel_t4"
)


# ------------------------------------------------------------------------------
# 8. Helper functions
# ------------------------------------------------------------------------------

q3 <- function(x) {
  unname(
    quantile(
      x,
      probs = c(0.025, 0.50, 0.975),
      names = FALSE
    )
  )
}

scalar_draws <- function(fit, variable) {
  x <- fit$draws(
    variables = variable,
    format = "matrix"
  )

  stopifnot(ncol(x) == 1)

  as.numeric(x[, 1])
}

scalar_summary <- function(fit, variable) {

  x <- scalar_draws(fit, variable)

  qs <- q3(x)

  data.frame(
    variable = variable,
    mean = mean(x),
    sd = sd(x),
    median = qs[2],
    q025 = qs[1],
    q975 = qs[3],
    row.names = NULL
  )
}

# Full Stan convergence diagnostics for selected parameters.
#
# MCSE is reported for the posterior mean.  The paper-level inferential
# summaries use medians and quantiles; the MCSE_mean is included because
# the requested MCMC reporting calls for a direct measure of simulation error.

mcmc_parameter_summary <- function(fit, variables) {

  fit$summary(
    variables = variables,

    "mean",
    "median",
    "sd",

    q025 = ~quantile(., 0.025),
    q975 = ~quantile(., 0.975),

    mcse_mean = posterior::mcse_mean,

    rhat = posterior::rhat,
    ess_bulk = posterior::ess_bulk,
    ess_tail = posterior::ess_tail
  )
}


model_results <- function(fit, model_name) {

  s_mu   <- scalar_summary(fit, "mu")
  s_tau  <- scalar_summary(fit, "tau")
  s_pool <- scalar_summary(fit, "PF_pool")
  s_pred <- scalar_summary(fit, "PF_new")

  b_any <- scalar_draws(fit, "benefit_new")
  b_20  <- scalar_draws(fit, "benefit20_new")

  ppc_ss   <- scalar_draws(fit, "ppc_ss")
  ppc_max  <- scalar_draws(fit, "ppc_max")
  ppc_mils <- scalar_draws(fit, "ppc_milsom_re_extreme")

  data.frame(
    model = model_name,

    mu_median = s_mu$median,
    mu_lb = s_mu$q025,
    mu_ub = s_mu$q975,

    PF_pool_median = s_pool$median,
    PF_pool_lb = s_pool$q025,
    PF_pool_ub = s_pool$q975,

    tau_median = s_tau$median,
    tau_lb = s_tau$q025,
    tau_ub = s_tau$q975,

    PF_pred_median = s_pred$median,
    PF_pred_lb = s_pred$q025,
    PF_pred_ub = s_pred$q975,

    Pr_PFnew_gt_0 = mean(b_any),
    Pr_PFnew_gt_020 = mean(b_20),

    PPC_p_SS = mean(ppc_ss),
    PPC_p_MAX = mean(ppc_max),

    Pr_new_RE_more_extreme_than_Milsom = mean(ppc_mils),

    row.names = NULL
  )
}


sampler_summary <- function(fit, model_name) {

  # HMC diagnostics by chain.
  d <- fit$diagnostic_summary(quiet = TRUE)

  # Convergence / efficiency diagnostics over ALL actual model parameters.
  # eta_C and z are vectors, so requesting their base names includes all
  # elements.
  ps <- mcmc_parameter_summary(
    fit,
    variables = c("mu", "tau", "eta_C", "z")
  )

  # MCSE for clinically reported quantities.
  rs <- mcmc_parameter_summary(
    fit,
    variables = c(
      "mu",
      "tau",
      "PF_pool",
      "PF_new",
      "benefit_new",
      "benefit20_new"
    )
  )

  data.frame(
    model = model_name,

    max_Rhat = max(ps$rhat, na.rm = TRUE),
    min_ESS_bulk = min(ps$ess_bulk, na.rm = TRUE),
    min_ESS_tail = min(ps$ess_tail, na.rm = TRUE),

    max_MCSE_mean_reported =
      max(rs$mcse_mean, na.rm = TRUE),

    divergences = sum(d$num_divergent),
    max_treedepth_hits = sum(d$num_max_treedepth),
    min_EBFMI = min(d$ebfmi),

    row.names = NULL
  )
}


study_effect_summary <- function(fit, model_name) {

  x <- fit$draws(
    variables = "PF_study",
    format = "matrix"
  )

  out <- do.call(
    rbind,
    lapply(
      seq_len(nrow(dat)),
      function(i) {

        qs <- q3(x[, paste0("PF_study[", i, "]")])

        data.frame(
          model = model_name,
          study = dat$study[i],
          PF_median = qs[2],
          PF_lb = qs[1],
          PF_ub = qs[3],
          row.names = NULL
        )
      }
    )
  )

  out
}


# ------------------------------------------------------------------------------
# 9. Production summaries
# ------------------------------------------------------------------------------

results <- rbind(
  model_results(
    fit_normal,
    "Arm-level Gaussian RE"
  ),
  model_results(
    fit_t4,
    "Arm-level standardized t4 RE"
  )
)

sampler <- rbind(
  sampler_summary(
    fit_normal,
    "Arm-level Gaussian RE"
  ),
  sampler_summary(
    fit_t4,
    "Arm-level standardized t4 RE"
  )
)

study_effects <- rbind(
  study_effect_summary(
    fit_normal,
    "Arm-level Gaussian RE"
  ),
  study_effect_summary(
    fit_t4,
    "Arm-level standardized t4 RE"
  )
)


cat("\n============================================================\n")
cat("POSTERIOR RESULTS\n")
cat("============================================================\n\n")

print(results, digits = 6, row.names = FALSE)

cat("\n============================================================\n")
cat("SAMPLER / CONVERGENCE DIAGNOSTICS\n")
cat("============================================================\n\n")

print(sampler, digits = 6, row.names = FALSE)


# ------------------------------------------------------------------------------
# 10. Detailed parameter summaries
# ------------------------------------------------------------------------------

summary_normal_reported <- mcmc_parameter_summary(
  fit_normal,
  variables = c(
    "mu",
    "tau",
    "PF_pool",
    "PF_new",
    "benefit_new",
    "benefit20_new",
    "theta",
    "PF_study"
  )
)

summary_t4_reported <- mcmc_parameter_summary(
  fit_t4,
  variables = c(
    "mu",
    "tau",
    "PF_pool",
    "PF_new",
    "benefit_new",
    "benefit20_new",
    "theta",
    "PF_study"
  )
)

write.csv(
  results,
  file.path(outdir, "stan_model_results.csv"),
  row.names = FALSE
)

write.csv(
  sampler,
  file.path(outdir, "stan_sampler_diagnostics.csv"),
  row.names = FALSE
)

write.csv(
  study_effects,
  file.path(outdir, "stan_study_specific_PF.csv"),
  row.names = FALSE
)

write.csv(
  summary_normal_reported,
  file.path(outdir, "stan_normal_full_summary.csv"),
  row.names = FALSE
)

write.csv(
  summary_t4_reported,
  file.path(outdir, "stan_t4_full_summary.csv"),
  row.names = FALSE
)

capture.output(
  fit_normal$diagnostic_summary(quiet = TRUE),
  file = file.path(outdir, "stan_normal_HMC_diagnostics.txt")
)

capture.output(
  fit_t4$diagnostic_summary(quiet = TRUE),
  file = file.path(outdir, "stan_t4_HMC_diagnostics.txt")
)

capture.output(
  fit_normal$cmdstan_diagnose(),
  file = file.path(outdir, "stan_normal_cmdstan_diagnose.txt")
)

capture.output(
  fit_t4$cmdstan_diagnose(),
  file = file.path(outdir, "stan_t4_cmdstan_diagnose.txt")
)

capture.output(
  sessionInfo(),
  file = file.path(outdir, "R_sessionInfo.txt")
)


# ------------------------------------------------------------------------------
# 11. Production-quality diagnostic warning
# ------------------------------------------------------------------------------

diagnostics_ok <- with(
  sampler,
  all(
    max_Rhat < 1.01 &
    min_ESS_bulk >= 400 &
    min_ESS_tail >= 400 &
    divergences == 0 &
    max_treedepth_hits == 0 &
    min_EBFMI >= 0.30
  )
)

if (!diagnostics_ok) {
  warning(
    "\nAt least one production diagnostic threshold was not met.\n",
    "Inspect stan_sampler_diagnostics.csv and the CmdStan diagnose files\n",
    "before using the generated manuscript text as final."
  )
} else {
  cat(
    "\nAll prespecified production diagnostic thresholds were satisfied.\n"
  )
}


# ------------------------------------------------------------------------------
# 12. Figures: trace plots
# ------------------------------------------------------------------------------

make_trace_pdf <- function(fit, filename, title_prefix) {

  d <- fit$draws(
    variables = c("mu", "tau", "PF_pool"),
    format = "df"
  )

  pdf(filename, width = 9, height = 8)

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)

  par(mfrow = c(3, 1), mar = c(3.5, 4.2, 2.4, 1.0))

  vars <- c("mu", "tau", "PF_pool")

  for (v in vars) {

    yr <- range(d[[v]])

    plot(
      NA,
      xlim = range(d$.iteration),
      ylim = yr,
      xlab = "Post-warmup iteration",
      ylab = v,
      main = paste(title_prefix, "-", v)
    )

    for (ch in sort(unique(d$.chain))) {
      ii <- d$.chain == ch

      lines(
        d$.iteration[ii],
        d[[v]][ii],
        lwd = 0.55,
        lty = ch
      )
    }
  }

  dev.off()
}

make_trace_pdf(
  fit_normal,
  file.path(outdir, "stan_normal_trace.pdf"),
  "Gaussian RE"
)

make_trace_pdf(
  fit_t4,
  file.path(outdir, "stan_t4_trace.pdf"),
  "standardized t4 RE"
)


# ------------------------------------------------------------------------------
# 13. Figures: arm-level posterior predictive checks
# ------------------------------------------------------------------------------

rep_summary <- function(fit) {

  C <- fit$draws(
    variables = "mean_C_rep",
    format = "matrix"
  )

  T <- fit$draws(
    variables = "mean_T_rep",
    format = "matrix"
  )

  out <- data.frame(
    label = c(
      paste0(dat$study, " C"),
      paste0(dat$study, " T")
    ),

    observed = c(
      dat$mean_C,
      dat$mean_T
    ),

    median = NA_real_,
    lb = NA_real_,
    ub = NA_real_
  )

  for (i in seq_len(nrow(dat))) {

    qc <- q3(C[, paste0("mean_C_rep[", i, "]")])
    qt <- q3(T[, paste0("mean_T_rep[", i, "]")])

    out$lb[i]     <- qc[1]
    out$median[i] <- qc[2]
    out$ub[i]     <- qc[3]

    j <- nrow(dat) + i

    out$lb[j]     <- qt[1]
    out$median[j] <- qt[2]
    out$ub[j]     <- qt[3]
  }

  out
}

plot_ppc_arms <- function(fit, filename, main) {

  x <- rep_summary(fit)

  pdf(filename, width = 9, height = 8)

  y <- rev(seq_len(nrow(x)))

  xr <- range(
    c(
      x$lb,
      x$ub,
      x$observed
    )
  )

  plot(
    NA,
    xlim = xr,
    ylim = c(0.5, nrow(x) + 0.5),
    yaxt = "n",
    xlab = "Reported / replicated arm mean",
    ylab = "",
    bty = "n",
    main = main
  )

  axis(
    2,
    at = y,
    labels = x$label,
    las = 1,
    cex.axis = 0.72,
    tick = FALSE
  )

  segments(
    x$lb,
    y,
    x$ub,
    y,
    lwd = 1.2
  )

  points(
    x$median,
    y,
    pch = 16
  )

  points(
    x$observed,
    y,
    pch = 1,
    cex = 1.05
  )

  legend(
    "bottomright",
    legend = c(
      "Posterior predictive median",
      "Observed mean"
    ),
    pch = c(16, 1),
    bty = "n"
  )

  dev.off()
}

plot_ppc_arms(
  fit_normal,
  file.path(outdir, "stan_normal_arm_PPC.pdf"),
  "Arm-level posterior predictive check: Gaussian RE"
)

plot_ppc_arms(
  fit_t4,
  file.path(outdir, "stan_t4_arm_PPC.pdf"),
  "Arm-level posterior predictive check: standardized t4 RE"
)


# ------------------------------------------------------------------------------
# 14. Figure: pooled and predictive PF distributions
# ------------------------------------------------------------------------------

pool_n <- scalar_draws(fit_normal, "PF_pool")
pred_n <- scalar_draws(fit_normal, "PF_new")

pool_t <- scalar_draws(fit_t4, "PF_pool")
pred_t <- scalar_draws(fit_t4, "PF_new")

pdf(
  file.path(outdir, "stan_PF_model_comparison.pdf"),
  width = 9,
  height = 7
)

oldpar <- par(no.readonly = TRUE)
par(mfrow = c(2, 1), mar = c(4, 4, 2.5, 1))

dn <- density(pool_n)
dt <- density(pool_t)

plot(
  dn,
  lwd = 2,
  xlab = "Pooled prevented fraction",
  ylab = "Posterior density",
  main = "Pooled historical effect"
)

lines(dt, lwd = 2, lty = 2)
abline(v = 0, lty = 3)

legend(
  "topleft",
  legend = c("Gaussian RE", "standardized t4 RE"),
  lty = c(1, 2),
  lwd = 2,
  bty = "n"
)

dn <- density(pred_n)
dt <- density(pred_t)

xr <- range(dn$x, dt$x)

plot(
  dn,
  lwd = 2,
  xlim = xr,
  xlab = "Prevented fraction in a new true study",
  ylab = "Posterior predictive density",
  main = "Prediction"
)

lines(dt, lwd = 2, lty = 2)
abline(v = 0, lty = 3)
abline(v = 0.20, lty = 3)

legend(
  "topleft",
  legend = c("Gaussian RE", "standardized t4 RE"),
  lty = c(1, 2),
  lwd = 2,
  bty = "n"
)

par(oldpar)
dev.off()


# ------------------------------------------------------------------------------
# 15. Figure: study-specific effects under both models
# ------------------------------------------------------------------------------

se_n <- subset(
  study_effects,
  model == "Arm-level Gaussian RE"
)

se_t <- subset(
  study_effects,
  model == "Arm-level standardized t4 RE"
)

raw_pf <- 1 - dat$mean_T / dat$mean_C

pdf(
  file.path(outdir, "stan_study_specific_PF.pdf"),
  width = 9,
  height = 6.5
)

y <- rev(seq_len(nrow(dat)))

xr <- range(
  c(
    se_n$PF_lb,
    se_n$PF_ub,
    se_t$PF_lb,
    se_t$PF_ub,
    raw_pf
  )
)

plot(
  NA,
  xlim = xr,
  ylim = c(0.5, nrow(dat) + 0.5),
  yaxt = "n",
  xlab = "Positive prevented fraction",
  ylab = "",
  bty = "n",
  main = "Study-specific posterior effects"
)

axis(
  2,
  at = y,
  labels = dat$study,
  las = 1,
  tick = FALSE
)

abline(v = 0, lty = 3)

off <- 0.13

segments(
  se_n$PF_lb, y + off,
  se_n$PF_ub, y + off,
  lwd = 1.1
)

points(
  se_n$PF_median,
  y + off,
  pch = 16
)

segments(
  se_t$PF_lb, y - off,
  se_t$PF_ub, y - off,
  lwd = 1.1,
  lty = 2
)

points(
  se_t$PF_median,
  y - off,
  pch = 17
)

points(
  raw_pf,
  y,
  pch = 1
)

legend(
  "bottomright",
  legend = c(
    "Gaussian RE posterior",
    "t4 RE posterior",
    "Raw arm-ratio PF"
  ),
  pch = c(16, 17, 1),
  lty = c(1, 2, NA),
  bty = "n"
)

dev.off()


# ------------------------------------------------------------------------------
# 16. Formatting helpers for LaTeX
# ------------------------------------------------------------------------------

fmt <- function(x, d = 3) {
  sprintf(paste0("%.", d, "f"), x)
}

fmt_int <- function(med, lb, ub, d = 3) {
  paste0(
    fmt(med, d),
    " [",
    fmt(lb, d),
    ",",
    fmt(ub, d),
    "]"
  )
}

tex_escape <- function(x) {
  gsub("_", "\\\\_", x, fixed = TRUE)
}


# ------------------------------------------------------------------------------
# 17. Write replacement rows for Table historical-models in manuscript
# ------------------------------------------------------------------------------

rn <- results[1, ]
rt <- results[2, ]

table_rows <- c(
  "Bayesian arm-level,",
  "normal random effects",
  paste0(
    "& ",
    fmt_int(
      rn$PF_pool_median,
      rn$PF_pool_lb,
      rn$PF_pool_ub
    )
  ),
  paste0("& ", fmt(rn$tau_median)),
  paste0(
    "& ",
    fmt_int(
      rn$PF_pred_median,
      rn$PF_pred_lb,
      rn$PF_pred_ub
    )
  ),
  paste0("& ", fmt(rn$Pr_PFnew_gt_0), " \\\\"),
  "",
  "Bayesian arm-level,",
  "$t_4$ random effects",
  paste0(
    "& ",
    fmt_int(
      rt$PF_pool_median,
      rt$PF_pool_lb,
      rt$PF_pool_ub
    )
  ),
  paste0("& ", fmt(rt$tau_median)),
  paste0(
    "& ",
    fmt_int(
      rt$PF_pred_median,
      rt$PF_pred_lb,
      rt$PF_pred_ub
    )
  ),
  paste0("& ", fmt(rt$Pr_PFnew_gt_0), " \\\\")
)

writeLines(
  table_rows,
  file.path(outdir, "stan_table_rows_for_manuscript.tex")
)


# ------------------------------------------------------------------------------
# 18. Generate compact main-paper replacement text
# ------------------------------------------------------------------------------

sn <- sampler[1, ]
st <- sampler[2, ]

diag_ok_text <- if (diagnostics_ok) {

  paste0(
    "Across the two fits, all rank-normalized $\\\\widehat R$ values for ",
    "the sampled model parameters were below ",
    fmt(max(sampler$max_Rhat), 3),
    ", the minimum bulk and tail effective sample sizes were ",
    sprintf("%.0f", min(sampler$min_ESS_bulk)),
    " and ",
    sprintf("%.0f", min(sampler$min_ESS_tail)),
    ", respectively, no divergent transitions or maximum-treedepth hits ",
    "were observed, and the minimum chain-specific E-BFMI was ",
    fmt(min(sampler$min_EBFMI), 2),
    "."
  )

} else {

  paste0(
    "\\textbf{The production diagnostic thresholds were not all met.} ",
    "The largest $\\\\widehat R$ was ",
    fmt(max(sampler$max_Rhat), 3),
    "; total divergences were ",
    sum(sampler$divergences),
    "; maximum-treedepth hits were ",
    sum(sampler$max_treedepth_hits),
    "; and the minimum chain-specific E-BFMI was ",
    fmt(min(sampler$min_EBFMI), 2),
    ". These diagnostics should be resolved before the numerical results ",
    "are treated as final."
  )
}

main_lines <- c(
  "\\paragraph{Stan arm-level production analyses.}",
  paste0(
    "The aggregate arm-level models were fit in Stan using four NUTS chains ",
    "with a noncentered parameterization. Under Gaussian random effects, ",
    "the posterior median pooled prevented fraction was ",
    fmt(rn$PF_pool_median),
    " with 95\\% credible interval ",
    "$[",
    fmt(rn$PF_pool_lb),
    ",",
    fmt(rn$PF_pool_ub),
    "]$, and the posterior median heterogeneity was ",
    "$\\\\tau=",
    fmt(rn$tau_median),
    "$. The posterior predictive prevented fraction for a new true study ",
    "had median ",
    fmt(rn$PF_pred_median),
    " and 95\\% predictive interval ",
    "$[",
    fmt(rn$PF_pred_lb),
    ",",
    fmt(rn$PF_pred_ub),
    "]$. The posterior probabilities of any preventive effect and of at ",
    "least a 20\\% preventive fraction were ",
    fmt(rn$Pr_PFnew_gt_0),
    " and ",
    fmt(rn$Pr_PFnew_gt_020),
    ", respectively."
  ),
  "",
  paste0(
    "Replacing the Gaussian random-effects distribution by the standardized ",
    "$t_4$ distribution left the pooled center similar: the posterior median ",
    "pooled prevented fraction was ",
    fmt(rt$PF_pool_median),
    " with 95\\% credible interval ",
    "$[",
    fmt(rt$PF_pool_lb),
    ",",
    fmt(rt$PF_pool_ub),
    "]$, while the posterior median heterogeneity was ",
    "$\\\\tau=",
    fmt(rt$tau_median),
    "$. The corresponding posterior predictive median was ",
    fmt(rt$PF_pred_median),
    " with 95\\% interval ",
    "$[",
    fmt(rt$PF_pred_lb),
    ",",
    fmt(rt$PF_pred_ub),
    "]$. The posterior probabilities of any preventive effect and of at ",
    "least a 20\\% effect were ",
    fmt(rt$Pr_PFnew_gt_0),
    " and ",
    fmt(rt$Pr_PFnew_gt_020),
    "."
  ),
  "",
  diag_ok_text,
  "",
  paste0(
    "Posterior predictive checks based on the sum of squared standardized ",
    "arm-level residuals gave tail areas ",
    fmt(rn$PPC_p_SS),
    " and ",
    fmt(rt$PPC_p_SS),
    " for the Gaussian and $t_4$ models, respectively; the corresponding ",
    "maximum-residual checks gave ",
    fmt(rn$PPC_p_MAX),
    " and ",
    fmt(rt$PPC_p_MAX),
    ". Detailed MCMC diagnostics and posterior predictive assessments are ",
    "reported in the Supplementary Material."
  )
)

writeLines(
  main_lines,
  file.path(outdir, "stan_main_results_for_manuscript.tex")
)


# ------------------------------------------------------------------------------
# 19. Generate extensive Supplementary Material section
# ------------------------------------------------------------------------------

# Milsom summaries
mils_n <- subset(
  study_effects,
  model == "Arm-level Gaussian RE" &
    study == "Milsom"
)

mils_t <- subset(
  study_effects,
  model == "Arm-level standardized t4 RE" &
    study == "Milsom"
)

# Study-specific LaTeX rows.
study_rows <- character(0)

for (i in seq_len(nrow(dat))) {

  a <- se_n[i, ]
  b <- se_t[i, ]

  study_rows <- c(
    study_rows,
    paste0(
      dat$study[i], " & ",
      fmt_int(a$PF_median, a$PF_lb, a$PF_ub), " & ",
      fmt_int(b$PF_median, b$PF_lb, b$PF_ub),
      " \\\\"
    )
  )
}

supp_lines <- c(
  "% ============================================================================",
  "% GENERATED by historical9_stan_armlevel_production.R",
  "% Do not edit numerical values by hand; rerun the production script instead.",
  "% ============================================================================",
  "",
  "\\section{Stan Production Analysis for the Aggregate Arm-Level Models}",
  "\\label{sec:supp-stan-armlevel}",
  "",
  paste0(
    "The aggregate arm-level Gaussian and robust $t_4$ random-effects models ",
    "were regenerated in Stan using the frozen nine-study arm summaries. ",
    "These models are complementary to the contrast-based normal--normal ",
    "hierarchy: rather than taking a first-stage log ratio and estimated ",
    "sampling variance as data, they model the reported treatment and control ",
    "sample means directly. Stan was used because the likelihood is nonlinear ",
    "in the latent baseline means and treatment effects and because Hamiltonian ",
    "Monte Carlo provides transparent convergence and geometry diagnostics.",
    "\\cite{Carpenter2017,Vehtari2021}"
  ),
  "",
  "\\subsection{Aggregate arm-level likelihood}",
  "",
  "For study $i$, let $m_i>0$ denote the latent control mean and let",
  "$\\theta_i$ denote the study-specific log treatment-to-control mean ratio.",
  "The sampling model is",
  "$$",
  "\\bar X_{Ci}\\mid m_i",
  "\\sim",
  "N\\!\\left(m_i,\\frac{s_{Ci}^2}{n_{Ci}}\\right),",
  "$$",
  "and",
  "$$",
  "\\bar X_{Ti}\\mid m_i,\\theta_i",
  "\\sim",
  "N\\!\\left(m_i e^{\\theta_i},\\frac{s_{Ti}^2}{n_{Ti}}\\right).",
  "$$",
  "Positivity of the latent baseline mean is imposed through",
  "$$",
  "m_i=\\exp(\\eta_i),",
  "\\qquad",
  "\\eta_i\\sim N(0,2.5^2).",
  "$$",
  "The common priors are",
  "$$",
  "\\mu\\sim N(0,1^2),",
  "\\qquad",
  "\\tau\\sim\\mathrm{Half\\mbox{-}Normal}(0,0.5^2).",
  "$$",
  "",
  "The reported arm standard deviations and sample sizes are treated as fixed.",
  "Consequently, this remains an aggregate-data approximation: it does not",
  "reconstruct zero inflation, individual-level longitudinal dependence,",
  "cluster effects that are absent from the published uncertainty, or",
  "individual-level covariates.",
  "",
  "\\subsection{Noncentered Gaussian and standardized $t_4$ random effects}",
  "",
  "For the Gaussian model, the study effects are written in noncentered form as",
  "$$",
  "\\theta_i=\\mu+\\tau z_i,",
  "\\qquad",
  "z_i\\sim N(0,1).",
  "$$",
  "The robust analysis retains the same definition of $\\mu$ and $\\tau$ but",
  "replaces the standardized Gaussian effect by",
  "$$",
  "\\theta_i=\\mu+\\tau z_i,",
  "\\qquad",
  "z_i\\sim t_4\\!\\left(0,\\sqrt{\\frac{4-2}{4}}\\right).",
  "$$",
  "For a Student-$t_\\nu$ distribution with scale $s$,",
  "$\\operatorname{Var}(Z)=s^2\\nu/(\\nu-2)$ when $\\nu>2$.",
  "Thus the scale $s=\\sqrt{(4-2)/4}=1/\\sqrt{2}$ makes",
  "$\\operatorname{Var}(z_i)=1$, so $\\tau$ retains the interpretation of the",
  "between-study standard deviation in both models. The degrees of freedom are",
  "fixed rather than estimated because nine studies contain little information",
  "about a separate tail-shape parameter.",
  "",
  "\\subsection{Posterior prediction and clinically interpretable targets}",
  "",
  "For every posterior draw, a new true study effect was generated from the",
  "same fitted random-effects population. The reported pooled and predictive",
  "quantities are",
  "$$",
  "\\PF_{+,\\mathrm{pool}}=1-e^\\mu,",
  "\\qquad",
  "\\PF_{+,\\mathrm{new}}=1-e^{\\theta_{\\mathrm{new}}}.",
  "$$",
  "We additionally estimate",
  "$$",
  "\\Pr(\\PF_{+,\\mathrm{new}}>0\\mid\\mathcal D)",
  "\\quad\\hbox{and}\\quad",
  "\\Pr(\\PF_{+,\\mathrm{new}}>0.20\\mid\\mathcal D).",
  "$$",
  "",
  "\\subsection{NUTS configuration and convergence criteria}",
  "",
  paste0(
    "Each model used four independently initialized NUTS chains, ",
    iter_warmup,
    " warm-up iterations per chain, and ",
    iter_sampling,
    " retained iterations per chain, for ",
    chains * iter_sampling,
    " post-warm-up draws per model. The target acceptance probability was ",
    fmt(adapt_delta, 3),
    " and the maximum tree depth was ",
    max_treedepth,
    ". A noncentered parameterization was used for the random effects."
  ),
  "",
  "Production diagnostics were evaluated using rank-normalized split",
  "$\\widehat R$, bulk and tail effective sample sizes, Monte Carlo standard",
  "errors, divergent transitions, maximum-treedepth hits, and E-BFMI.",
  "\\cite{Vehtari2021}",
  "The prespecified numerical checks were $\\widehat R<1.01$, bulk and tail",
  "ESS at least 400 for every sampled model parameter, no divergences, no",
  "maximum-treedepth hits, and chain-specific E-BFMI at least 0.30.",
  "",
  "\\begin{table}[ht]",
  "\\centering",
  "\\caption{Posterior summaries from the Stan aggregate arm-level models. Entries are posterior medians with central 95\\% credible or predictive intervals.}",
  "\\label{tab:supp-stan-results}",
  "\\small",
  "\\begin{tabular}{p{3.8cm}cc}",
  "\\toprule",
  "Quantity & Gaussian RE & Standardized $t_4$ RE \\\\",
  "\\midrule",
  paste0(
    "Pooled $\\PF_+$ & ",
    fmt_int(rn$PF_pool_median, rn$PF_pool_lb, rn$PF_pool_ub),
    " & ",
    fmt_int(rt$PF_pool_median, rt$PF_pool_lb, rt$PF_pool_ub),
    " \\\\"
  ),
  paste0(
    "$\\mu$ & ",
    fmt_int(rn$mu_median, rn$mu_lb, rn$mu_ub),
    " & ",
    fmt_int(rt$mu_median, rt$mu_lb, rt$mu_ub),
    " \\\\"
  ),
  paste0(
    "$\\tau$ & ",
    fmt_int(rn$tau_median, rn$tau_lb, rn$tau_ub),
    " & ",
    fmt_int(rt$tau_median, rt$tau_lb, rt$tau_ub),
    " \\\\"
  ),
  paste0(
    "Predictive $\\PF_{+,\\mathrm{new}}$ & ",
    fmt_int(rn$PF_pred_median, rn$PF_pred_lb, rn$PF_pred_ub),
    " & ",
    fmt_int(rt$PF_pred_median, rt$PF_pred_lb, rt$PF_pred_ub),
    " \\\\"
  ),
  paste0(
    "$\\Pr(\\PF_{+,\\mathrm{new}}>0)$ & ",
    fmt(rn$Pr_PFnew_gt_0),
    " & ",
    fmt(rt$Pr_PFnew_gt_0),
    " \\\\"
  ),
  paste0(
    "$\\Pr(\\PF_{+,\\mathrm{new}}>0.20)$ & ",
    fmt(rn$Pr_PFnew_gt_020),
    " & ",
    fmt(rt$Pr_PFnew_gt_020),
    " \\\\"
  ),
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}",
  "",
  "\\subsection{Hamiltonian Monte Carlo diagnostics}",
  "",
  diag_ok_text,
  "",
  "\\begin{table}[ht]",
  "\\centering",
  "\\caption{Stan convergence and Hamiltonian Monte Carlo diagnostics. ESS values are minima over all sampled model parameters; $\\widehat R$ is the corresponding maximum.}",
  "\\label{tab:supp-stan-diagnostics}",
  "\\small",
  "\\begin{tabular}{lcc}",
  "\\toprule",
  "Diagnostic & Gaussian RE & Standardized $t_4$ RE \\\\",
  "\\midrule",
  paste0(
    "Maximum $\\widehat R$ & ",
    fmt(sn$max_Rhat, 3),
    " & ",
    fmt(st$max_Rhat, 3),
    " \\\\"
  ),
  paste0(
    "Minimum bulk ESS & ",
    sprintf("%.0f", sn$min_ESS_bulk),
    " & ",
    sprintf("%.0f", st$min_ESS_bulk),
    " \\\\"
  ),
  paste0(
    "Minimum tail ESS & ",
    sprintf("%.0f", sn$min_ESS_tail),
    " & ",
    sprintf("%.0f", st$min_ESS_tail),
    " \\\\"
  ),
  paste0(
    "Largest MCSE(mean), reported targets & ",
    fmt(sn$max_MCSE_mean_reported, 4),
    " & ",
    fmt(st$max_MCSE_mean_reported, 4),
    " \\\\"
  ),
  paste0(
    "Divergent transitions & ",
    sn$divergences,
    " & ",
    st$divergences,
    " \\\\"
  ),
  paste0(
    "Maximum-treedepth hits & ",
    sn$max_treedepth_hits,
    " & ",
    st$max_treedepth_hits,
    " \\\\"
  ),
  paste0(
    "Minimum E-BFMI & ",
    fmt(sn$min_EBFMI, 3),
    " & ",
    fmt(st$min_EBFMI, 3),
    " \\\\"
  ),
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}",
  "",
  "\\begin{figure}[ht]",
  "\\centering",
  "\\includegraphics[width=0.88\\textwidth]{output_03_historical_Stan/stan_normal_trace.pdf}",
  "\\caption{Post-warm-up trace plots for $\\mu$, $\\tau$, and the pooled prevented fraction under the Gaussian random-effects arm-level model.}",
  "\\label{fig:supp-stan-normal-trace}",
  "\\end{figure}",
  "",
  "\\begin{figure}[ht]",
  "\\centering",
  "\\includegraphics[width=0.88\\textwidth]{output_03_historical_Stan/stan_t4_trace.pdf}",
  "\\caption{Post-warm-up trace plots for $\\mu$, $\\tau$, and the pooled prevented fraction under the standardized $t_4$ random-effects arm-level model.}",
  "\\label{fig:supp-stan-t4-trace}",
  "\\end{figure}",
  "",
  "\\subsection{Posterior predictive checks}",
  "",
  "Replicated treatment and control sample means were generated from the",
  "fitted arm-level sampling model for every posterior draw. Two global",
  "discrepancy statistics were used. The first is the sum of squared",
  "standardized residuals over all 18 observed arm means,",
  "$$",
  "T_{\\mathrm{SS}}",
  "=",
  "\\sum_{i=1}^{9}",
  "\\left\\{",
  "\\left(\\frac{\\bar X_{Ci}-m_i}{s_{Ci}/\\sqrt{n_{Ci}}}\\right)^2",
  "+",
  "\\left(\\frac{\\bar X_{Ti}-m_i e^{\\theta_i}}{s_{Ti}/\\sqrt{n_{Ti}}}\\right)^2",
  "\\right\\}.",
  "$$",
  "The second is the maximum absolute standardized arm-level residual,",
  "$$",
  "T_{\\max}",
  "=",
  "\\max_{i,a\\in\\{C,T\\}} |r_{ia}|.",
  "$$",
  "For each discrepancy, the posterior predictive tail area is the fraction",
  "of posterior draws for which the replicated statistic exceeds its observed",
  "counterpart.",
  "",
  "\\begin{table}[ht]",
  "\\centering",
  "\\caption{Posterior predictive checks for the aggregate arm-level models. Values near 0 or 1 would indicate that the observed discrepancy lies in a tail of the replicated distribution.}",
  "\\label{tab:supp-stan-ppc}",
  "\\small",
  "\\begin{tabular}{lcc}",
  "\\toprule",
  "Check & Gaussian RE & Standardized $t_4$ RE \\\\",
  "\\midrule",
  paste0(
    "$\\Pr(T_{\\mathrm{SS}}^{\\mathrm{rep}}\\geq T_{\\mathrm{SS}}^{\\mathrm{obs}}\\mid\\mathcal D)$ & ",
    fmt(rn$PPC_p_SS),
    " & ",
    fmt(rt$PPC_p_SS),
    " \\\\"
  ),
  paste0(
    "$\\Pr(T_{\\max}^{\\mathrm{rep}}\\geq T_{\\max}^{\\mathrm{obs}}\\mid\\mathcal D)$ & ",
    fmt(rn$PPC_p_MAX),
    " & ",
    fmt(rt$PPC_p_MAX),
    " \\\\"
  ),
  paste0(
    "$\\Pr(|z_{\\mathrm{new}}|\\geq |z_{\\mathrm{Milsom}}|\\mid\\mathcal D)$ & ",
    fmt(rn$Pr_new_RE_more_extreme_than_Milsom),
    " & ",
    fmt(rt$Pr_new_RE_more_extreme_than_Milsom),
    " \\\\"
  ),
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}",
  "",
  paste0(
    "Under the Gaussian model, the study-specific posterior prevented fraction ",
    "for Milsom had median ",
    fmt(mils_n$PF_median),
    " and 95\\% credible interval $[",
    fmt(mils_n$PF_lb),
    ",",
    fmt(mils_n$PF_ub),
    "]$; under the standardized $t_4$ model the corresponding summary was ",
    fmt(mils_t$PF_median),
    " $[",
    fmt(mils_t$PF_lb),
    ",",
    fmt(mils_t$PF_ub),
    "]$. The additional tail-area calculation in Table~\\ref{tab:supp-stan-ppc} ",
    "places Milsom's latent standardized random effect relative to a fresh draw ",
    "from the fitted population rather than labeling the study mechanically as ",
    "an outlier."
  ),
  "",
  "\\begin{figure}[ht]",
  "\\centering",
  "\\includegraphics[width=0.86\\textwidth]{output_03_historical_Stan/stan_normal_arm_PPC.pdf}",
  "\\caption{Arm-level posterior predictive intervals under Gaussian random effects. Open points are the observed arm means; filled points and horizontal intervals are posterior predictive medians and central 95\\% intervals.}",
  "\\label{fig:supp-stan-normal-ppc}",
  "\\end{figure}",
  "",
  "\\begin{figure}[ht]",
  "\\centering",
  "\\includegraphics[width=0.86\\textwidth]{output_03_historical_Stan/stan_t4_arm_PPC.pdf}",
  "\\caption{Arm-level posterior predictive intervals under standardized $t_4$ random effects.}",
  "\\label{fig:supp-stan-t4-ppc}",
  "\\end{figure}",
  "",
  "\\subsection{Study-specific posterior effects}",
  "",
  "\\begin{table}[ht]",
  "\\centering",
  "\\caption{Study-specific posterior prevented fractions under the two aggregate arm-level models. Entries are posterior medians with central 95\\% credible intervals.}",
  "\\label{tab:supp-stan-study-effects}",
  "\\small",
  "\\begin{tabular}{lcc}",
  "\\toprule",
  "Study & Gaussian RE & Standardized $t_4$ RE \\\\",
  "\\midrule",
  study_rows,
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}",
  "",
  "\\begin{figure}[ht]",
  "\\centering",
  "\\includegraphics[width=0.84\\textwidth]{output_03_historical_Stan/stan_study_specific_PF.pdf}",
  "\\caption{Study-specific posterior prevented fractions from the Gaussian and standardized $t_4$ arm-level models, together with the raw prevented fractions calculated from the reported arm means.}",
  "\\label{fig:supp-stan-study-effects}",
  "\\end{figure}",
  "",
  "\\subsection{Robustness of the pooled center versus predictive tails}",
  "",
  paste0(
    "The Gaussian model gives a pooled posterior median prevented fraction of ",
    fmt(rn$PF_pool_median),
    ", whereas the standardized $t_4$ model gives ",
    fmt(rt$PF_pool_median),
    ". The corresponding posterior median heterogeneity estimates are ",
    fmt(rn$tau_median),
    " and ",
    fmt(rt$tau_median),
    ". The more consequential difference is therefore assessed through the ",
    "posterior predictive distribution rather than through the pooled center. ",
    "The Gaussian predictive interval is $[",
    fmt(rn$PF_pred_lb),
    ",",
    fmt(rn$PF_pred_ub),
    "]$, compared with $[",
    fmt(rt$PF_pred_lb),
    ",",
    fmt(rt$PF_pred_ub),
    "]$ under standardized $t_4$ random effects."
  ),
  "",
  "\\begin{figure}[ht]",
  "\\centering",
  "\\includegraphics[width=0.84\\textwidth]{output_03_historical_Stan/stan_PF_model_comparison.pdf}",
  "\\caption{Posterior distributions of the pooled prevented fraction (top) and posterior predictive distributions for a new true study effect (bottom) under Gaussian and standardized $t_4$ random effects.}",
  "\\label{fig:supp-stan-pf-comparison}",
  "\\end{figure}",
  "",
  "\\subsection{Reproducibility}",
  "",
  paste0(
    "The production run used R version ",
    R.version.string,
    ", CmdStanR version ",
    as.character(packageVersion("cmdstanr")),
    ", posterior version ",
    as.character(packageVersion("posterior")),
    ", and CmdStan version ",
    as.character(cmdstan_ver),
    ". The random-number seeds were 20260927 for the Gaussian model and ",
    "20260928 for the standardized $t_4$ model."
  ),
  "",
  "The archived production directory contains the original CmdStan CSV files,",
  "the complete parameter summaries, HMC diagnostic summaries, CmdStan",
  "diagnostic output, posterior predictive summaries, figures, the input data,",
  "and the R session information. Re-running the driver script regenerates the",
  "main-manuscript text and every numerical table in this section from the",
  "same frozen model fits."
)

writeLines(
  supp_lines,
  file.path(outdir, "stan_supplement_results.tex")
)


# ------------------------------------------------------------------------------
# 20. Version / configuration record
# ------------------------------------------------------------------------------

config <- c(
  paste0("R: ", R.version.string),
  paste0("cmdstanr: ", packageVersion("cmdstanr")),
  paste0("posterior: ", packageVersion("posterior")),
  paste0("CmdStan: ", cmdstan_ver),
  paste0("chains: ", chains),
  paste0("iter_warmup: ", iter_warmup),
  paste0("iter_sampling: ", iter_sampling),
  paste0("adapt_delta: ", adapt_delta),
  paste0("max_treedepth: ", max_treedepth),
  "seed normal: 20260927",
  "seed t4: 20260928"
)

writeLines(
  config,
  file.path(outdir, "stan_run_configuration.txt")
)


# ------------------------------------------------------------------------------
# 21. Independent check against the earlier provisional values
# ------------------------------------------------------------------------------

# These are NOT hard constraints.  They are only a gross check that the
# production fit has not accidentally changed model, data, or transformation.

benchmark <- data.frame(
  model = c(
    "Arm-level Gaussian RE",
    "Arm-level standardized t4 RE"
  ),
  PF_pool_provisional = c(0.432, 0.425),
  tau_provisional = c(0.385, 0.431),
  P_any_provisional = c(0.915, 0.909)
)

check <- data.frame(
  model = results$model,
  PF_pool_difference =
    results$PF_pool_median - benchmark$PF_pool_provisional,
  tau_difference =
    results$tau_median - benchmark$tau_provisional,
  P_any_difference =
    results$Pr_PFnew_gt_0 - benchmark$P_any_provisional
)

cat("\n============================================================\n")
cat("CHECK AGAINST EARLIER PROVISIONAL VALUES\n")
cat("============================================================\n\n")

print(check, digits = 6, row.names = FALSE)

write.csv(
  check,
  file.path(outdir, "stan_check_vs_provisional.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 22. Final console message
# ------------------------------------------------------------------------------

cat("\n============================================================\n")
cat("PRODUCTION RUN COMPLETE\n")
cat("============================================================\n\n")

cat(
  "Main manuscript text:\n  ",
  file.path(outdir, "stan_main_results_for_manuscript.tex"),
  "\n\n",
  sep = ""
)

cat(
  "Replacement table rows for manuscript:\n  ",
  file.path(outdir, "stan_table_rows_for_manuscript.tex"),
  "\n\n",
  sep = ""
)

cat(
  "Extensive Supplement section:\n  ",
  file.path(outdir, "stan_supplement_results.tex"),
  "\n\n",
  sep = ""
)

cat(
  "All output is in:\n  ",
  normalizePath(outdir),
  "\n",
  sep = ""
)

# ==============================================================================
# END
# ==============================================================================
