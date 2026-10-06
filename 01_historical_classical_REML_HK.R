# ==============================================================================
# 01_historical_classical_REML_HK.R
#
# Frequentist benchmark for the 9 historical fluoride-varnish trials
# Effect measure: log ratio of means (log-RoM / log response ratio)
# Estimator: REML random-effects meta-analysis
# Inference: Hartung-Knapp (Knapp-Hartung) small-sample adjustment
# Prediction interval: Riley-type t_{k-2} interval, to match the manuscript
#
# Primary software:
#   R + metafor
#
# This script is intended to be the reproducible frequentist benchmark against
# which the Bayesian hierarchical analyses are compared.
#
# IMPORTANT DATA NOTE
# -------------------
# The arm-level values below are the defined 9-study analysis dataset used in
# the current manuscript:
#     treatment/control means, SDs, and sample sizes.
#
# The script derives the log ratio and its first-order delta-method variance
# directly from those arm summaries. This is preferable to entering rounded
# LRRs and SEs from a manuscript table.
#
# Before final submission, any study-level design corrections (especially
# cluster randomization) should be incorporated at the within-study variance
# stage if a valid design-adjusted SE is available. The present code reproduces
# the currently defined benchmark.
#
# Expected benchmark from these exact arm summaries:
#   REML mu (log-RoM)            about -0.5050
#   REML tau                     about  0.2835
#   pooled PF+                   about  0.3965
#   HK 95% CI for pooled PF+     about [0.180, 0.556]
#   Riley 95% PI for new PF+     about [-0.265, 0.712]
#
# ==============================================================================


# ------------------------------------------------------------------------------
# 0. Package
# ------------------------------------------------------------------------------

# Install once if needed:
# install.packages("metafor")

library(metafor)


# ------------------------------------------------------------------------------
# 1. Enter the 9-study arm-level summary data
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

  # Fluoride-varnish arm
  mean_T = c(0.90, 2.50, 2.43, 0.55, 1.48, 0.79, 4.61, 0.33, 0.66),
  sd_T   = c(3.80, 3.10, 3.09, 4.48, 1.53, 1.67, 4.39, 1.04, 0.73),
  n_T    = c(60,   87,   246,  311,  98,   190,  113,  91,   94),

  # Control arm
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
# 2. Construct the log ratio of means and sampling variance
# ------------------------------------------------------------------------------

# Effect:
#
#       y_i = log(mean_T / mean_C)
#
# First-order delta-method variance:
#
#       v_i = sd_T^2/(n_T mean_T^2) + sd_C^2/(n_C mean_C^2)
#
# Positive prevented fraction used only for interpretation:
#
#       PF+ = 1 - exp(y)
#
# Negative y favors fluoride varnish; positive PF+ denotes prevention.

dat$yi <- log(dat$mean_T / dat$mean_C)

dat$vi <- (dat$sd_T^2 / (dat$n_T * dat$mean_T^2)) +
          (dat$sd_C^2 / (dat$n_C * dat$mean_C^2))

dat$sei <- sqrt(dat$vi)
dat$PF  <- 1 - exp(dat$yi)


# Independent metafor construction of the same effect measure.
# measure="ROM" is the log ratio of means.
#
# IMPORTANT (metafor >= 5.1-16):
# metafor now applies a second-order bias correction to measure="ROM" by
# default.  Our primary manuscript estimand is the ordinary first-order
# log ratio log(mean_T/mean_C), with the usual large-sample (LS) delta
# variance.  Therefore we must explicitly set correct=FALSE and vtype="LS"
# to reproduce the analysis defined in the manuscript.
#
# This serves as a check that the manually stated formula and metafor agree.

dat_es <- escalc(
  measure = "ROM",
  m1i = mean_T, sd1i = sd_T, n1i = n_T,
  m2i = mean_C, sd2i = sd_C, n2i = n_C,
  data = dat,
  vtype = "LS",
  correct = FALSE
)

stopifnot(
  max(abs(dat$yi - dat_es$yi)) < 1e-10,
  max(abs(dat$vi - dat_es$vi)) < 1e-10
)


# Descriptive, study-specific 95% normal intervals.
# These are NOT the Hartung-Knapp pooled interval; they are simply useful
# study-level summaries on the log-RoM scale.

z975 <- qnorm(0.975)

dat$yi_ci_lb <- dat$yi - z975 * dat$sei
dat$yi_ci_ub <- dat$yi + z975 * dat$sei

# Because PF+(y)=1-exp(y) is decreasing, interval endpoints reverse.
dat$PF_ci_lb <- 1 - exp(dat$yi_ci_ub)
dat$PF_ci_ub <- 1 - exp(dat$yi_ci_lb)

cat("\n============================================================\n")
cat("STUDY-LEVEL DATA AND DERIVED EFFECTS\n")
cat("============================================================\n\n")

print(
  dat[, c(
    "study", "mean_T", "sd_T", "n_T",
    "mean_C", "sd_C", "n_C",
    "yi", "sei", "PF", "PF_ci_lb", "PF_ci_ub"
  )],
  digits = 4,
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 3. Primary frequentist benchmark: REML + Hartung-Knapp
# ------------------------------------------------------------------------------

# method="REML" estimates tau^2 by restricted maximum likelihood.
#
# test="knha" applies the Knapp-Hartung / Hartung-Knapp adjustment to
# the uncertainty of the pooled mean. With k=9 and an intercept-only model,
# the pooled-effect t reference distribution has k-p = 8 df.

fit_hk <- rma.uni(
  yi = yi,
  vi = vi,
  data = dat,
  method = "REML",
  test = "knha",
  level = 95
)

cat("\n============================================================\n")
cat("PRIMARY MODEL: REML + HARTUNG-KNAPP\n")
cat("============================================================\n\n")

print(fit_hk, digits = 5)


# ------------------------------------------------------------------------------
# 4. Check the Hartung-Knapp scale factor
# ------------------------------------------------------------------------------

# For an intercept-only model, the usual HK variance estimator can be written
#
#       q_HK / sum(w_i),
#
# where
#
#       q_HK = [ sum w_i (y_i - mu_hat)^2 ] / (k-1)
#       w_i  = 1 / (v_i + tau_hat^2).
#
# If q_HK < 1, the ordinary HK adjustment can yield a smaller SE than the
# conventional REML/Wald SE. The "adhoc" modification in metafor prevents that
# by using the larger SE.
#
# Here q_HK > 1, so ordinary HK and the modified safeguard should coincide.

k <- nrow(dat)
p <- 1L

tau2_hat <- as.numeric(fit_hk$tau2)
tau_hat  <- sqrt(tau2_hat)
mu_hat   <- as.numeric(fit_hk$b[1])

w_re <- 1 / (dat$vi + tau2_hat)

q_HK <- sum(w_re * (dat$yi - mu_hat)^2) / (k - p)
se_HK_manual <- sqrt(q_HK / sum(w_re))

fit_adhoc <- rma.uni(
  yi = yi,
  vi = vi,
  data = dat,
  method = "REML",
  test = "adhoc",
  level = 95
)

cat("\nHartung-Knapp scale factor q_HK =", format(q_HK, digits = 6), "\n")
cat("Manual HK SE                  =", format(se_HK_manual, digits = 6), "\n")
cat("metafor HK SE                 =", format(as.numeric(fit_hk$se), digits = 6), "\n")
cat("metafor modified-HK SE        =", format(as.numeric(fit_adhoc$se), digits = 6), "\n")

if (q_HK > 1) {
  cat("q_HK > 1: ordinary HK and modified HK should coincide.\n")
} else {
  cat("q_HK < 1: modified HK may be wider than ordinary HK.\n")
}


# ------------------------------------------------------------------------------
# 5. Transform pooled estimate and HK confidence interval to PF+
# ------------------------------------------------------------------------------

# PF+(y) = 1-exp(y) is monotonically decreasing.
# Therefore:
#   lower PF bound = 1-exp(upper log bound)
#   upper PF bound = 1-exp(lower log bound)

pf_from_log <- function(x) {
  1 - exp(x)
}

reverse_transform_interval <- function(lb, ub) {
  c(
    lb = 1 - exp(unname(ub)),
    ub = 1 - exp(unname(lb))
  )
}

pooled_pf <- pf_from_log(mu_hat)

ci_log <- c(
  lb = as.numeric(fit_hk$ci.lb),
  ub = as.numeric(fit_hk$ci.ub)
)

ci_pf <- reverse_transform_interval(ci_log["lb"], ci_log["ub"])


# ------------------------------------------------------------------------------
# 6. Prediction interval for a new true study effect
# ------------------------------------------------------------------------------

# IMPORTANT:
# The manuscript's prediction interval uses the Riley-style small-sample
# interval. For a random-effects model with one pooled-effect parameter, this
# uses t_{k-2}, here t_7, and combines tau_hat^2 with the HK uncertainty of
# mu_hat.
#
# metafor's default prediction interval with test="knha" uses t_{k-p}=t_8.
# The current manuscript values correspond to predtype="Riley", so we use that
# as the principal benchmark and print the metafor default only for comparison.

pred_riley <- predict(
  fit_hk,
  level = 95,
  predtype = "Riley"
)

pred_default <- predict(
  fit_hk,
  level = 95
)

pi_log <- c(
  lb = as.numeric(pred_riley$pi.lb),
  ub = as.numeric(pred_riley$pi.ub)
)

pi_pf <- reverse_transform_interval(pi_log["lb"], pi_log["ub"])


# Manual Riley calculation as an independent check:
#
#   mu_hat +/- t_{k-2, .975} sqrt(tau_hat^2 + SE_HK^2)

crit_riley <- qt(0.975, df = k - 2)

pi_log_manual <- c(
  lb = mu_hat - crit_riley * sqrt(tau2_hat + se_HK_manual^2),
  ub = mu_hat + crit_riley * sqrt(tau2_hat + se_HK_manual^2)
)

stopifnot(
  max(abs(pi_log - pi_log_manual)) < 1e-8
)


# ------------------------------------------------------------------------------
# 7. Heterogeneity summaries
# ------------------------------------------------------------------------------

Q  <- as.numeric(fit_hk$QE)
Qp <- as.numeric(fit_hk$QEp)
I2 <- as.numeric(fit_hk$I2)
H2 <- as.numeric(fit_hk$H2)

cat("\n============================================================\n")
cat("PRIMARY BENCHMARK SUMMARY\n")
cat("============================================================\n\n")

cat("k                         =", k, "\n")
cat("Pooled log-RoM mu_hat     =", sprintf("%.6f", mu_hat), "\n")
cat("REML tau^2                =", sprintf("%.6f", tau2_hat), "\n")
cat("REML tau                  =", sprintf("%.6f", tau_hat), "\n")
cat("HK pooled SE              =", sprintf("%.6f", as.numeric(fit_hk$se)), "\n")
cat("Pooled PF+                =", sprintf("%.6f", pooled_pf), "\n")
cat("HK 95% CI for log-RoM     = [",
    sprintf("%.6f", ci_log["lb"]), ", ",
    sprintf("%.6f", ci_log["ub"]), "]\n", sep = "")
cat("HK 95% CI for PF+         = [",
    sprintf("%.6f", ci_pf["lb"]), ", ",
    sprintf("%.6f", ci_pf["ub"]), "]\n", sep = "")
cat("Riley 95% PI for log-RoM  = [",
    sprintf("%.6f", pi_log["lb"]), ", ",
    sprintf("%.6f", pi_log["ub"]), "]\n", sep = "")
cat("Riley 95% PI for PF+      = [",
    sprintf("%.6f", pi_pf["lb"]), ", ",
    sprintf("%.6f", pi_pf["ub"]), "]\n", sep = "")
cat("Cochran Q                 =", sprintf("%.4f", Q),
    "  p =", sprintf("%.5f", Qp), "\n")
cat("I^2 (%)                   =", sprintf("%.2f", I2), "\n")
cat("H^2                       =", sprintf("%.3f", H2), "\n")


# ------------------------------------------------------------------------------
# 8. Numerical reproduction check against the current manuscript
# ------------------------------------------------------------------------------

# These are checks, not inputs to the analysis.
# Tolerances allow for ordinary floating-point / package-level differences.

expected <- c(
  mu       = -0.5050,
  tau      =  0.2835,
  pooledPF =  0.3965,
  ciPF_lb  =  0.180,
  ciPF_ub  =  0.556,
  piPF_lb  = -0.265,
  piPF_ub  =  0.712
)

obtained <- c(
  mu       = mu_hat,
  tau      = tau_hat,
  pooledPF = pooled_pf,
  ciPF_lb  = ci_pf["lb"],
  ciPF_ub  = ci_pf["ub"],
  piPF_lb  = pi_pf["lb"],
  piPF_ub  = pi_pf["ub"]
)

cat("\n============================================================\n")
cat("CHECK AGAINST MANUSCRIPT BENCHMARK\n")
cat("============================================================\n\n")

check_table <- data.frame(
  quantity = names(expected),
  manuscript = as.numeric(expected),
  reproduced = as.numeric(obtained),
  difference = as.numeric(obtained - expected)
)

print(check_table, digits = 6, row.names = FALSE)


# ------------------------------------------------------------------------------
# 9. Leave-one-study-out REML-HK influence analysis
# ------------------------------------------------------------------------------

# We refit the complete model after deleting each study. This is intentionally
# done in a loop so that the same Riley prediction interval is obtained for
# every deletion and transformed consistently to PF+.

loo_one <- function(j, data) {

  d <- data[-j, , drop = FALSE]

  f <- rma.uni(
    yi = yi,
    vi = vi,
    data = d,
    method = "REML",
    test = "knha",
    level = 95
  )

  pr <- predict(
    f,
    level = 95,
    predtype = "Riley"
  )

  mu_j  <- as.numeric(f$b[1])
  tau_j <- sqrt(as.numeric(f$tau2))

  ci_j <- reverse_transform_interval(
    as.numeric(f$ci.lb),
    as.numeric(f$ci.ub)
  )

  pi_j <- reverse_transform_interval(
    as.numeric(pr$pi.lb),
    as.numeric(pr$pi.ub)
  )

  data.frame(
    omitted = data$study[j],
    k = nrow(d),
    mu_log = mu_j,
    pooled_PF = pf_from_log(mu_j),
    CI_PF_lb = ci_j["lb"],
    CI_PF_ub = ci_j["ub"],
    tau = tau_j,
    PI_PF_lb = pi_j["lb"],
    PI_PF_ub = pi_j["ub"],
    row.names = NULL
  )
}

loo <- do.call(
  rbind,
  lapply(seq_len(nrow(dat)), loo_one, data = dat)
)

cat("\n============================================================\n")
cat("LEAVE-ONE-STUDY-OUT REML-HK ANALYSIS\n")
cat("============================================================\n\n")

print(loo, digits = 4, row.names = FALSE)


# ------------------------------------------------------------------------------
# 10. Profile-likelihood confidence interval for heterogeneity
# ------------------------------------------------------------------------------

# This is supplementary. It is useful to show how uncertain tau is with k=9.
# confint() reports intervals for tau^2 and tau (depending on metafor version /
# print method). We print the result rather than hard-code extraction.

cat("\n============================================================\n")
cat("PROFILE-LIKELIHOOD INTERVAL FOR HETEROGENEITY\n")
cat("============================================================\n\n")

tau_ci <- confint(fit_hk, type = "PL")
print(tau_ci, digits = 5)


# ------------------------------------------------------------------------------
# 11. Save machine-readable results
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

outdir <- file.path(script_dir, "output_01_historical_classical")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

write.csv(
  dat,
  file = file.path(outdir, "historical9_study_effects.csv"),
  row.names = FALSE
)

summary_out <- data.frame(
  k = k,
  mu_log = mu_hat,
  tau2_REML = tau2_hat,
  tau_REML = tau_hat,
  HK_SE = as.numeric(fit_hk$se),
  pooled_PF = pooled_pf,
  pooled_log_CI_lb = ci_log["lb"],
  pooled_log_CI_ub = ci_log["ub"],
  pooled_PF_CI_lb = ci_pf["lb"],
  pooled_PF_CI_ub = ci_pf["ub"],
  Riley_log_PI_lb = pi_log["lb"],
  Riley_log_PI_ub = pi_log["ub"],
  Riley_PF_PI_lb = pi_pf["lb"],
  Riley_PF_PI_ub = pi_pf["ub"],
  Q = Q,
  Q_p = Qp,
  I2_percent = I2,
  H2 = H2,
  HK_scale_q = q_HK
)

write.csv(
  summary_out,
  file = file.path(outdir, "historical9_REML_HK_summary.csv"),
  row.names = FALSE
)

write.csv(
  loo,
  file = file.path(outdir, "historical9_leave_one_out.csv"),
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = file.path(outdir, "R_sessionInfo.txt")
)


# ------------------------------------------------------------------------------
# 12. Publication-style forest plot on the positive prevented-fraction scale
# ------------------------------------------------------------------------------

# Study intervals are descriptive first-order normal intervals.
# The pooled interval is the HK 95% confidence interval.
# The pooled diamond is plotted at the transformed pooled effect.

pdf(
  file.path(outdir, "historical9_REML_HK_forest_PF.pdf"),
  width = 8.5,
  height = 6.5
)

op <- par(
  mar = c(5, 9, 2, 2),
  las = 1,
  xaxs = "i"
)

ypos <- rev(seq_len(nrow(dat))) + 1
xlim <- c(-0.50, 1.00)

plot(
  NA,
  xlim = xlim,
  ylim = c(0.5, nrow(dat) + 2.0),
  xlab = "Prevented fraction  PF+ = 1 - treatment/control mean ratio",
  ylab = "",
  yaxt = "n",
  bty = "n"
)

abline(v = 0, lty = 2)

axis(
  side = 2,
  at = ypos,
  labels = dat$study,
  tick = FALSE,
  las = 1
)

# Study-specific intervals and estimates
segments(
  x0 = dat$PF_ci_lb,
  y0 = ypos,
  x1 = dat$PF_ci_ub,
  y1 = ypos,
  lwd = 1.2
)

points(
  dat$PF,
  ypos,
  pch = 15,
  cex = 0.9
)

# Pooled HK row
ypool <- 0.9

axis(
  side = 2,
  at = ypool,
  labels = "REML-HK",
  tick = FALSE,
  las = 1
)

segments(
  x0 = ci_pf["lb"],
  y0 = ypool,
  x1 = ci_pf["ub"],
  y1 = ypool,
  lwd = 2
)

# Diamond for pooled estimate and HK CI
diamond_y <- c(ypool, ypool + 0.16, ypool, ypool - 0.16)
diamond_x <- c(ci_pf["lb"], pooled_pf, ci_pf["ub"], pooled_pf)

polygon(
  diamond_x,
  diamond_y
)

title(
  main = "Nine historical fluoride-varnish trials: REML-Hartung-Knapp"
)

par(op)
dev.off()


# ------------------------------------------------------------------------------
# 13. Optional diagnostic plots
# ------------------------------------------------------------------------------

# Standard metafor influence diagnostics. These are useful for the Supplement
# but do not need to appear in the main paper.

pdf(
  file.path(outdir, "historical9_influence_diagnostics.pdf"),
  width = 9,
  height = 7
)

infl <- influence(fit_hk)
plot(infl)

dev.off()


# ------------------------------------------------------------------------------
# 14. Final concise console statement
# ------------------------------------------------------------------------------

cat("\n============================================================\n")
cat("FINAL MANUSCRIPT-SCALE BENCHMARK\n")
cat("============================================================\n\n")

cat(
  sprintf(
    paste0(
      "REML-HK: mu = %.4f, tau = %.4f, pooled PF+ = %.4f ",
      "(95%% HK CI %.3f to %.3f); Riley 95%% prediction interval ",
      "%.3f to %.3f.\n"
    ),
    mu_hat,
    tau_hat,
    pooled_pf,
    ci_pf["lb"],
    ci_pf["ub"],
    pi_pf["lb"],
    pi_pf["ub"]
  )
)

cat("\nOutput files written to: ", normalizePath(outdir), "\n", sep = "")

# ==============================================================================
# END
# ==============================================================================
