# ==============================================================================
# 02_historical_bayesmeta_prior_sensitivity.R
#
# Bayesian normal-normal hierarchical meta-analysis for the 9 historical
# fluoride-varnish trials.
#
# PURPOSE
# -------
# This is the production analysis for the contrast-based Bayesian NNHM:
#
#   y_i | theta_i  ~ N(theta_i, s_i^2)
#   theta_i | mu,tau ~ N(mu, tau^2)
#
# with
#
#   mu ~ N(0, 1^2)
#
# and four targeted priors for tau:
#
#   1. Primary:          Half-Normal(scale = 0.50)
#   2. More concentrated Half-Normal(scale = 0.25)
#   3. Heavier tailed:   Half-t_4(scale = 0.50)
#   4. Alternative shape: DuMouchel prior
#
# IMPORTANT COMPUTATIONAL POINT
# -----------------------------
# bayesmeta uses deterministic numerical integration for this low-dimensional
# normal-normal hierarchical model. There are NO MCMC chains, and therefore
# R-hat, ESS, MCSE, divergences, treedepth, etc. are not relevant for this
# analysis. Those diagnostics are reserved for the later Stan arm-level and
# heavy-tailed random-effects models.
#
# Primary outputs:
#   * posterior median and 95% central CrI for mu
#   * posterior median and 95% central CrI for tau
#   * pooled prevented fraction PF+ = 1-exp(mu)
#   * posterior predictive PF+ for a new true study effect
#   * Pr(PF_new > 0 | data)
#   * Pr(PF_new > 0.20 | data)
#   * primary-model study-specific shrinkage summaries
#   * prior/posterior sensitivity plots
#
# Expected values from an independent deterministic check are approximately:
#
# Primary Half-N(0,0.5^2):
#   pooled PF median 0.389, 95% CrI [0.215, 0.556]
#   tau median 0.302, 95% CrI [0.100, 0.661]
#   predictive PF median 0.386, 95% PI [-0.284, 0.732]
#   Pr(PF_new>0) 0.927; Pr(PF_new>0.20) 0.809
#
# Minor differences in the last decimal may arise from package numerical
# tolerances. The CSV files produced below are the production source.
# ==============================================================================


# ------------------------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------------------------

# Install once if needed:
# install.packages("bayesmeta")

library(bayesmeta)


# ------------------------------------------------------------------------------
# 1. Nine historical trials: arm-level summary data
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
  all(dat$sd_T >= 0),
  all(dat$sd_C >= 0),
  all(dat$n_T > 1),
  all(dat$n_C > 1)
)


# ------------------------------------------------------------------------------
# 2. First-stage log ratios and first-order delta-method standard errors
# ------------------------------------------------------------------------------

# Primary contrast:
#
#   y_i = log(mean_T / mean_C)
#
# First-order variance:
#
#   s_i^2 =
#     sd_T^2/(n_T mean_T^2) +
#     sd_C^2/(n_C mean_C^2)
#
# The second-order response-ratio correction is intentionally NOT used here.
# It will be analyzed separately as a sensitivity analysis.

dat$yi <- log(dat$mean_T / dat$mean_C)

dat$vi <- (dat$sd_T^2 / (dat$n_T * dat$mean_T^2)) +
          (dat$sd_C^2 / (dat$n_C * dat$mean_C^2))

dat$sei <- sqrt(dat$vi)

cat("\n============================================================\n")
cat("NINE-STUDY LOG-RATIO INPUT DATA\n")
cat("============================================================\n\n")

print(
  dat[, c("study", "yi", "sei")],
  digits = 7,
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 3. Priors
# ------------------------------------------------------------------------------

# Common prior for pooled log ratio:
#
#   mu ~ N(0, 1^2)
#
# Heterogeneity-prior sensitivity:
#
# Primary:
tau_HN_050 <- function(tau) {
  bayesmeta::dhalfnormal(tau, scale = 0.50)
}

# More concentrated:
tau_HN_025 <- function(tau) {
  bayesmeta::dhalfnormal(tau, scale = 0.25)
}

# Heavier tailed:
tau_HT4_050 <- function(tau) {
  bayesmeta::dhalft(tau, df = 4, scale = 0.50)
}

# DuMouchel is supplied directly by bayesmeta via tau.prior="DuMouchel".
#
# Its data-dependent scale is
#
#   s0^2 = k / sum(1/s_i^2)
#
# and
#
#   p(tau) = s0 / (s0 + tau)^2, tau >= 0.

k <- nrow(dat)
s0_DuMouchel <- sqrt(k / sum(1 / dat$vi))

cat("\nDuMouchel scale s0 =", sprintf("%.7f", s0_DuMouchel), "\n")


# ------------------------------------------------------------------------------
# 4. Fit the four Bayesian NNHMs
# ------------------------------------------------------------------------------

# interval.type="central" ensures equal-tailed 95% intervals, matching the
# reporting used in the manuscript.
#
# delta and epsilon are set somewhat tighter than the defaults for a production
# analysis. With only nine studies, the computational burden remains small.

fit_model <- function(tau_prior) {

  bayesmeta(
    y = dat$yi,
    sigma = dat$sei,
    labels = dat$study,

    mu.prior.mean = 0,
    mu.prior.sd   = 1,

    tau.prior = tau_prior,

    interval.type = "central",

    delta = 0.005,
    epsilon = 1e-6
  )
}

cat("\nFitting primary Half-Normal(0,0.5^2)...\n")
fit_HN050 <- fit_model(tau_HN_050)

cat("Fitting concentrated Half-Normal(0,0.25^2)...\n")
fit_HN025 <- fit_model(tau_HN_025)

cat("Fitting heavier-tailed Half-t4(scale=0.5)...\n")
fit_HT4 <- fit_model(tau_HT4_050)

cat("Fitting DuMouchel prior...\n")
fit_DM <- fit_model("DuMouchel")


fits <- list(
  "Half-N(0,0.5^2)"  = fit_HN050,
  "Half-N(0,0.25^2)" = fit_HN025,
  "Half-t4(0,0.5)"    = fit_HT4,
  "DuMouchel"         = fit_DM
)


# ------------------------------------------------------------------------------
# 5. Helper functions for transformed posterior summaries
# ------------------------------------------------------------------------------

pf_from_log <- function(x) {
  1 - exp(x)
}

# q is assumed to contain the .025, .50, .975 quantiles of a log-scale
# parameter. Since PF=1-exp(x) is decreasing, the transformed endpoints reverse.

transform_log_quantiles_to_PF <- function(q) {

  stopifnot(length(q) == 3)

  c(
    PF_lb     = pf_from_log(q[3]),
    PF_median = pf_from_log(q[2]),
    PF_ub     = pf_from_log(q[1])
  )
}


summarize_bayesmeta <- function(fit, model_name) {

  probs <- c(0.025, 0.50, 0.975)

  q_mu <- unname(
    fit$qposterior(mu.p = probs)
  )

  q_tau <- unname(
    fit$qposterior(tau.p = probs)
  )

  # Posterior predictive distribution for a NEW TRUE study effect theta_new.
  # This is not a future noisy observed y; it is the distribution of the true
  # underlying study effect.
  q_pred_theta <- unname(
    fit$qposterior(
      theta.p = probs,
      predict = TRUE
    )
  )

  PF_pool <- transform_log_quantiles_to_PF(q_mu)
  PF_pred <- transform_log_quantiles_to_PF(q_pred_theta)

  # PF_new > 0  <=> theta_new < 0
  p_any <- unname(
    fit$pposterior(
      theta = 0,
      predict = TRUE
    )
  )

  # PF_new > 0.20 <=> 1-exp(theta_new)>0.20
  #                  <=> theta_new < log(0.80)
  p_20 <- unname(
    fit$pposterior(
      theta = log(0.80),
      predict = TRUE
    )
  )

  data.frame(
    model = model_name,

    mu_lb = q_mu[1],
    mu_median = q_mu[2],
    mu_ub = q_mu[3],

    PF_pool_lb = PF_pool["PF_lb"],
    PF_pool_median = PF_pool["PF_median"],
    PF_pool_ub = PF_pool["PF_ub"],

    tau_lb = q_tau[1],
    tau_median = q_tau[2],
    tau_ub = q_tau[3],

    theta_pred_lb = q_pred_theta[1],
    theta_pred_median = q_pred_theta[2],
    theta_pred_ub = q_pred_theta[3],

    PF_pred_lb = PF_pred["PF_lb"],
    PF_pred_median = PF_pred["PF_median"],
    PF_pred_ub = PF_pred["PF_ub"],

    Pr_PFnew_gt_0 = p_any,
    Pr_PFnew_gt_020 = p_20,

    row.names = NULL
  )
}


results <- do.call(
  rbind,
  Map(summarize_bayesmeta, fits, names(fits))
)

cat("\n============================================================\n")
cat("BAYESMETA PRIOR-SENSITIVITY RESULTS\n")
cat("============================================================\n\n")

print(results, digits = 6, row.names = FALSE)


# ------------------------------------------------------------------------------
# 6. Compact manuscript table
# ------------------------------------------------------------------------------

paper_table <- results[, c(
  "model",
  "PF_pool_median", "PF_pool_lb", "PF_pool_ub",
  "tau_median", "tau_lb", "tau_ub",
  "PF_pred_median", "PF_pred_lb", "PF_pred_ub",
  "Pr_PFnew_gt_0", "Pr_PFnew_gt_020"
)]

cat("\n============================================================\n")
cat("COMPACT PAPER TABLE\n")
cat("============================================================\n\n")

print(paper_table, digits = 4, row.names = FALSE)


# ------------------------------------------------------------------------------
# 7. Study-specific posterior shrinkage under the PRIMARY prior
# ------------------------------------------------------------------------------

probs <- c(0.025, 0.50, 0.975)

shrinkage <- do.call(
  rbind,
  lapply(seq_len(nrow(dat)), function(i) {

    q_theta <- unname(
      fit_HN050$qposterior(
        theta.p = probs,
        individual = i
      )
    )

    q_pf <- transform_log_quantiles_to_PF(q_theta)

    data.frame(
      study = dat$study[i],
      raw_log_ratio = dat$yi[i],
      raw_PF = pf_from_log(dat$yi[i]),
      theta_lb = q_theta[1],
      theta_median = q_theta[2],
      theta_ub = q_theta[3],
      PF_lb = q_pf["PF_lb"],
      PF_median = q_pf["PF_median"],
      PF_ub = q_pf["PF_ub"],
      row.names = NULL
    )
  })
)

cat("\n============================================================\n")
cat("PRIMARY-PRIOR STUDY-SPECIFIC SHRINKAGE\n")
cat("============================================================\n\n")

print(shrinkage, digits = 5, row.names = FALSE)


# ------------------------------------------------------------------------------
# 8. Prior calibration summaries
# ------------------------------------------------------------------------------

# For the proper fixed-scale priors, report prior median and 95th percentile.
# For DuMouchel:
#   F(tau) = tau/(s0+tau)
# so q_p = s0*p/(1-p).

prior_summary <- data.frame(
  prior = c(
    "Half-N(0,0.5^2)",
    "Half-N(0,0.25^2)",
    "Half-t4(0,0.5)",
    "DuMouchel"
  ),
  median = c(
    bayesmeta::qhalfnormal(0.50, scale = 0.50),
    bayesmeta::qhalfnormal(0.50, scale = 0.25),
    bayesmeta::qhalft(0.50, df = 4, scale = 0.50),
    s0_DuMouchel
  ),
  q95 = c(
    bayesmeta::qhalfnormal(0.95, scale = 0.50),
    bayesmeta::qhalfnormal(0.95, scale = 0.25),
    bayesmeta::qhalft(0.95, df = 4, scale = 0.50),
    s0_DuMouchel * 0.95 / 0.05
  )
)

cat("\n============================================================\n")
cat("HETEROGENEITY PRIOR CALIBRATION\n")
cat("============================================================\n\n")

print(prior_summary, digits = 5, row.names = FALSE)


# ------------------------------------------------------------------------------
# 9. Production output directory
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

outdir <- file.path(script_dir, "output_02_historical_bayesmeta")

if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

write.csv(
  dat,
  file = file.path(outdir, "bayesmeta_input_data.csv"),
  row.names = FALSE
)

write.csv(
  results,
  file = file.path(outdir, "bayesmeta_prior_sensitivity_full.csv"),
  row.names = FALSE
)

write.csv(
  paper_table,
  file = file.path(outdir, "bayesmeta_prior_sensitivity_paper.csv"),
  row.names = FALSE
)

write.csv(
  shrinkage,
  file = file.path(outdir, "bayesmeta_primary_shrinkage.csv"),
  row.names = FALSE
)

write.csv(
  prior_summary,
  file = file.path(outdir, "bayesmeta_tau_prior_calibration.csv"),
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = file.path(outdir, "R_sessionInfo.txt")
)


# ------------------------------------------------------------------------------
# 10. Plot: prior densities for tau
# ------------------------------------------------------------------------------

tau_max <- max(
  1.5,
  fit_HN050$qposterior(tau.p = 0.995),
  fit_HN025$qposterior(tau.p = 0.995),
  fit_HT4$qposterior(tau.p = 0.995),
  fit_DM$qposterior(tau.p = 0.995)
)

tau_grid <- seq(0, tau_max, length.out = 700)

pdf(
  file.path(outdir, "bayesmeta_tau_prior_densities.pdf"),
  width = 8,
  height = 6
)

plot(
  tau_grid,
  fit_HN050$dprior(tau = tau_grid),
  type = "l",
  lwd = 2,
  xlab = expression(tau),
  ylab = "Prior density",
  main = "Heterogeneity priors used in the sensitivity analysis"
)

lines(
  tau_grid,
  fit_HN025$dprior(tau = tau_grid),
  lwd = 2,
  lty = 2
)

lines(
  tau_grid,
  fit_HT4$dprior(tau = tau_grid),
  lwd = 2,
  lty = 3
)

lines(
  tau_grid,
  fit_DM$dprior(tau = tau_grid),
  lwd = 2,
  lty = 4
)

legend(
  "topright",
  legend = names(fits),
  lty = 1:4,
  lwd = 2,
  bty = "n"
)

dev.off()


# ------------------------------------------------------------------------------
# 11. Plot: posterior densities for tau
# ------------------------------------------------------------------------------

pdf(
  file.path(outdir, "bayesmeta_tau_posterior_densities.pdf"),
  width = 8,
  height = 6
)

plot(
  tau_grid,
  fit_HN050$dposterior(tau = tau_grid),
  type = "l",
  lwd = 2,
  xlab = expression(tau),
  ylab = "Posterior density",
  main = "Posterior heterogeneity under alternative priors"
)

lines(
  tau_grid,
  fit_HN025$dposterior(tau = tau_grid),
  lwd = 2,
  lty = 2
)

lines(
  tau_grid,
  fit_HT4$dposterior(tau = tau_grid),
  lwd = 2,
  lty = 3
)

lines(
  tau_grid,
  fit_DM$dposterior(tau = tau_grid),
  lwd = 2,
  lty = 4
)

legend(
  "topright",
  legend = names(fits),
  lty = 1:4,
  lwd = 2,
  bty = "n"
)

dev.off()


# ------------------------------------------------------------------------------
# 12. Plot: posterior densities for the pooled effect mu
# ------------------------------------------------------------------------------

mu_lo <- min(sapply(fits, function(f) f$qposterior(mu.p = 0.001)))
mu_hi <- max(sapply(fits, function(f) f$qposterior(mu.p = 0.999)))

mu_grid <- seq(mu_lo, mu_hi, length.out = 700)

pdf(
  file.path(outdir, "bayesmeta_mu_posterior_densities.pdf"),
  width = 8,
  height = 6
)

plot(
  mu_grid,
  fit_HN050$dposterior(mu = mu_grid),
  type = "l",
  lwd = 2,
  xlab = expression(mu),
  ylab = "Posterior density",
  main = "Posterior pooled log ratio under alternative heterogeneity priors"
)

lines(
  mu_grid,
  fit_HN025$dposterior(mu = mu_grid),
  lwd = 2,
  lty = 2
)

lines(
  mu_grid,
  fit_HT4$dposterior(mu = mu_grid),
  lwd = 2,
  lty = 3
)

lines(
  mu_grid,
  fit_DM$dposterior(mu = mu_grid),
  lwd = 2,
  lty = 4
)

abline(v = 0, lty = 3)

legend(
  "topleft",
  legend = names(fits),
  lty = 1:4,
  lwd = 2,
  bty = "n"
)

dev.off()


# ------------------------------------------------------------------------------
# 13. Plot: posterior predictive densities on the PF scale
# ------------------------------------------------------------------------------

# If p = PF = 1-exp(theta), then theta=log(1-p) and
#
#   f_PF(p) = f_theta(log(1-p)) / (1-p),   p < 1.
#
# The plotted range is chosen to show the scientifically relevant mass.

pf_grid <- seq(-0.75, 0.90, length.out = 800)
theta_from_pf <- log(1 - pf_grid)

pred_pf_density <- function(fit) {
  fit$dposterior(
    theta = theta_from_pf,
    predict = TRUE
  ) / (1 - pf_grid)
}

dens_HN050 <- pred_pf_density(fit_HN050)
dens_HN025 <- pred_pf_density(fit_HN025)
dens_HT4   <- pred_pf_density(fit_HT4)
dens_DM    <- pred_pf_density(fit_DM)

pdf(
  file.path(outdir, "bayesmeta_predictive_PF_densities.pdf"),
  width = 8,
  height = 6
)

plot(
  pf_grid,
  dens_HN050,
  type = "l",
  lwd = 2,
  xlab = expression(PF["+,"*"new"]),
  ylab = "Posterior predictive density",
  main = "Posterior prediction under alternative heterogeneity priors"
)

lines(pf_grid, dens_HN025, lwd = 2, lty = 2)
lines(pf_grid, dens_HT4,   lwd = 2, lty = 3)
lines(pf_grid, dens_DM,    lwd = 2, lty = 4)

abline(v = 0, lty = 3)
abline(v = 0.20, lty = 3)

legend(
  "topleft",
  legend = names(fits),
  lty = 1:4,
  lwd = 2,
  bty = "n"
)

dev.off()


# ------------------------------------------------------------------------------
# 14. Plot: primary-prior raw and shrinkage estimates on PF scale
# ------------------------------------------------------------------------------

pdf(
  file.path(outdir, "bayesmeta_primary_shrinkage_PF.pdf"),
  width = 8.5,
  height = 6.5
)

ypos <- rev(seq_len(nrow(shrinkage)))

xlim <- range(
  c(
    shrinkage$PF_lb,
    shrinkage$PF_ub,
    shrinkage$raw_PF
  )
)

plot(
  NA,
  xlim = xlim,
  ylim = c(0.5, nrow(shrinkage) + 0.5),
  yaxt = "n",
  ylab = "",
  xlab = "Positive prevented fraction",
  bty = "n"
)

axis(
  2,
  at = ypos,
  labels = shrinkage$study,
  las = 1,
  tick = FALSE
)

abline(v = 0, lty = 3)

segments(
  shrinkage$PF_lb,
  ypos,
  shrinkage$PF_ub,
  ypos,
  lwd = 1.4
)

points(
  shrinkage$PF_median,
  ypos,
  pch = 16
)

points(
  shrinkage$raw_PF,
  ypos,
  pch = 1,
  cex = 1.1
)

legend(
  "bottomright",
  legend = c("Posterior median", "Raw study estimate"),
  pch = c(16, 1),
  bty = "n"
)

title("Study-specific shrinkage under the primary Half-Normal(0,0.5^2) prior")

dev.off()


# ------------------------------------------------------------------------------
# 15. Write a LaTeX table from the production output
# ------------------------------------------------------------------------------

fmt_interval <- function(med, lb, ub, digits = 3) {
  paste0(
    sprintf(paste0("%.", digits, "f"), med),
    " [",
    sprintf(paste0("%.", digits, "f"), lb),
    ",",
    sprintf(paste0("%.", digits, "f"), ub),
    "]"
  )
}

latex_lines <- c(
  "\\begin{table}[ht]",
  "\\centering",
  "\\caption{Prior sensitivity of the Bayesian normal--normal hierarchical model for the nine historical trials. Entries are posterior medians; intervals are central 95\\% credible or predictive intervals.}",
  "\\label{tab:bayesmeta-prior-sensitivity}",
  "\\small",
  "\\begin{tabular}{p{3.2cm}ccccc}",
  "\\toprule",
  "Prior for $\\tau$ & Pooled $\\PF_+$ & $\\tau$ & Predictive $\\PF_{+,\\mathrm{new}}$ & $\\Pr(\\PF_{+,\\mathrm{new}}>0)$ & $\\Pr(\\PF_{+,\\mathrm{new}}>0.20)$ \\\\",
  "\\midrule"
)

for (i in seq_len(nrow(results))) {

  prior_tex <- c(
    "Half-N$(0,0.5^2)$",
    "Half-N$(0,0.25^2)$",
    "Half-$t_4(0,0.5)$",
    "DuMouchel"
  )[i]

  latex_lines <- c(
    latex_lines,
    paste0(
      prior_tex, " & ",
      fmt_interval(
        results$PF_pool_median[i],
        results$PF_pool_lb[i],
        results$PF_pool_ub[i]
      ), " & ",
      fmt_interval(
        results$tau_median[i],
        results$tau_lb[i],
        results$tau_ub[i]
      ), " & ",
      fmt_interval(
        results$PF_pred_median[i],
        results$PF_pred_lb[i],
        results$PF_pred_ub[i]
      ), " & ",
      sprintf("%.3f", results$Pr_PFnew_gt_0[i]), " & ",
      sprintf("%.3f", results$Pr_PFnew_gt_020[i]),
      " \\\\"
    )
  )
}

latex_lines <- c(
  latex_lines,
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}"
)

writeLines(
  latex_lines,
  con = file.path(outdir, "bayesmeta_prior_sensitivity_table.tex")
)


# ------------------------------------------------------------------------------
# 16. Independent benchmark checks against expected production values
# ------------------------------------------------------------------------------

# These checks are deliberately loose enough to accommodate tiny numerical
# differences across bayesmeta versions, while still catching a wrong prior,
# wrong effect scale, or reversed PF transformation.

expected <- data.frame(
  model = c(
    "Half-N(0,0.5^2)",
    "Half-N(0,0.25^2)",
    "Half-t4(0,0.5)",
    "DuMouchel"
  ),
  PF_pool = c(0.389, 0.382, 0.389, 0.380),
  tau = c(0.302, 0.253, 0.300, 0.249),
  PF_pred = c(0.386, 0.379, 0.386, 0.375),
  p_any = c(0.927, 0.951, 0.927, 0.947),
  p_20 = c(0.809, 0.838, 0.809, 0.839)
)

check <- data.frame(
  model = results$model,
  PF_pool = results$PF_pool_median,
  tau = results$tau_median,
  PF_pred = results$PF_pred_median,
  p_any = results$Pr_PFnew_gt_0,
  p_20 = results$Pr_PFnew_gt_020
)

cat("\n============================================================\n")
cat("CHECK AGAINST INDEPENDENT DETERMINISTIC BENCHMARK\n")
cat("============================================================\n\n")

check_print <- data.frame(
  model = check$model,
  PF_pool_diff = check$PF_pool - expected$PF_pool,
  tau_diff = check$tau - expected$tau,
  PF_pred_diff = check$PF_pred - expected$PF_pred,
  p_any_diff = check$p_any - expected$p_any,
  p_20_diff = check$p_20 - expected$p_20
)

print(check_print, digits = 6, row.names = FALSE)

stopifnot(
  max(abs(check_print$PF_pool_diff)) < 0.01,
  max(abs(check_print$tau_diff)) < 0.01,
  max(abs(check_print$PF_pred_diff)) < 0.01,
  max(abs(check_print$p_any_diff)) < 0.01,
  max(abs(check_print$p_20_diff)) < 0.01
)


# ------------------------------------------------------------------------------
# 17. Final concise production statement
# ------------------------------------------------------------------------------

cat("\n============================================================\n")
cat("PRIMARY BAYESMETA RESULT\n")
cat("============================================================\n\n")

r <- results[1, ]

cat(
  sprintf(
    paste0(
      "Primary Half-N(0,0.5^2): pooled PF median %.3f ",
      "(95%% CrI %.3f to %.3f); tau median %.3f ",
      "(95%% CrI %.3f to %.3f); predictive PF median %.3f ",
      "(95%% PI %.3f to %.3f); Pr(PF_new>0)=%.3f; ",
      "Pr(PF_new>0.20)=%.3f.\n"
    ),
    r$PF_pool_median,
    r$PF_pool_lb,
    r$PF_pool_ub,
    r$tau_median,
    r$tau_lb,
    r$tau_ub,
    r$PF_pred_median,
    r$PF_pred_lb,
    r$PF_pred_ub,
    r$Pr_PFnew_gt_0,
    r$Pr_PFnew_gt_020
  )
)

cat("\nOutput written to: ", normalizePath(outdir), "\n", sep = "")

# ==============================================================================
# END
# ==============================================================================
