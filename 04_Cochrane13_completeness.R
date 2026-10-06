# ==============================================================================
# 04_Cochrane13_completeness.R
#
# Historical completeness sensitivity analysis:
#   9-study historical subset versus the 13-study permanent-surface
#   fluoride-varnish evidence set in Marinho et al. (Cochrane 2013).
#
# IMPORTANT:
# This is NOT the primary log-ratio-of-means analysis.
#
# The purpose of this script is narrower: place the original 9-study subset
# and the complete Cochrane-13 permanent-surface set on the SAME historical
# prevented-fraction (PF) scale and ask whether adding the four omitted
# historical trials materially changes the pooled historical center.
#
# The study-specific PF estimates and SEs below are transcribed from the
# Cochrane Analysis 1.1 forest plot:
#   "D(M)FS increment (prevented fraction - nearest to 3 years (13 trials))".
#
# Positive PF favors fluoride varnish.
#
# The four studies not in the 9-study core are:
#   Holm 1984, Borutta 1991, Hardman 2007, Gugwad 2011.
#
# Software:
#   R + metafor
#
# Expected REML point estimates from these Cochrane-scale data:
#   9-study subset:  mu ~= 0.4394, tau ~= 0.2169
#   13-study set:    mu ~= 0.4346, tau ~= 0.1974
#
# These reproduce the manuscript statement that adding the four studies
# changes the pooled historical center by less than one percentage point.
# ==============================================================================


# ------------------------------------------------------------------------------
# 0. Package
# ------------------------------------------------------------------------------

# Install once if needed:
# install.packages("metafor")

library(metafor)


# ------------------------------------------------------------------------------
# 1. Cochrane permanent-surface prevented-fraction data
# ------------------------------------------------------------------------------

# These values are the study-level PF estimates and standard errors shown in
# the Cochrane permanent-tooth-surface forest plot. They are intentionally
# entered directly rather than reconstructed from our current arm-level
# summaries, because this sensitivity analysis is designed to compare 9 and 13
# studies on exactly the same historical Cochrane scale.

dat <- data.frame(
  study = c(
    "Koch 1975",
    "Modeer 1984",
    "Holm 1984",
    "Clark 1985",
    "Tewari 1990",
    "Borutta 1991",
    "Bravo 1997",
    "Skold 2005",
    "Hardman 2007",
    "Tagliaferro 2011",
    "Gugwad 2011",
    "Milsom 2011",
    "Arruda 2012"
  ),

  PF = c(
     0.78,
     0.30,
     0.55,
     0.20,
     0.75,
     0.29,
     0.43,
     0.60,
    -0.18,
     0.42,
     1.88,
     0.03,
     0.40
  ),

  SE = c(
    0.126,
    0.154,
    0.078,
    0.077,
    0.124,
    0.106,
    0.141,
    0.060,
    0.451,
    0.245,
    0.832,
    0.144,
    0.091
  ),

  # TRUE = one of the original 9 historical studies.
  core9 = c(
    TRUE,  # Koch
    TRUE,  # Modeer
    FALSE, # Holm
    TRUE,  # Clark
    TRUE,  # Tewari
    FALSE, # Borutta
    TRUE,  # Bravo
    TRUE,  # Skold
    FALSE, # Hardman
    TRUE,  # Tagliaferro
    FALSE, # Gugwad
    TRUE,  # Milsom
    TRUE   # Arruda
  )
)

dat$vi <- dat$SE^2

stopifnot(
  nrow(dat) == 13,
  sum(dat$core9) == 9,
  all(dat$SE > 0)
)

dat9  <- dat[dat$core9, ]
dat13 <- dat


# ------------------------------------------------------------------------------
# 2. Fit REML random-effects models
# ------------------------------------------------------------------------------

# We fit the same normal-normal PF-scale model to both evidence sets:
#
#   z_i | phi_i ~ N(phi_i, s_i^2)
#   phi_i | mu,tau ~ N(mu, tau^2)
#
# This PF-scale Gaussian model is NOT the preferred primary model in the paper.
# It is used only for this historical completeness comparison.

fit_reml <- function(d, test = "knha") {
  rma.uni(
    yi = PF,
    sei = SE,
    data = d,
    method = "REML",
    test = test,
    level = 95
  )
}

fit9_hk  <- fit_reml(dat9,  test = "knha")
fit13_hk <- fit_reml(dat13, test = "knha")

# Modified Hartung-Knapp safeguard ("adhoc") is also computed.
# For the 9-study PF-scale fit, the ordinary HK scale factor is slightly below
# one, so the modified interval is marginally wider. For the 13-study fit, the
# ordinary and modified intervals coincide.

fit9_mhk  <- fit_reml(dat9,  test = "adhoc")
fit13_mhk <- fit_reml(dat13, test = "adhoc")


# ------------------------------------------------------------------------------
# 3. Riley-type small-sample prediction intervals
# ------------------------------------------------------------------------------

pred9 <- predict(
  fit9_mhk,
  level = 95,
  predtype = "Riley"
)

pred13 <- predict(
  fit13_mhk,
  level = 95,
  predtype = "Riley"
)


# ------------------------------------------------------------------------------
# 4. Hartung-Knapp scale factor
# ------------------------------------------------------------------------------

hk_scale <- function(fit, d) {

  tau2 <- as.numeric(fit$tau2)
  mu   <- as.numeric(fit$b[1])

  w <- 1 / (d$vi + tau2)

  sum(w * (d$PF - mu)^2) / (nrow(d) - 1)
}

q9  <- hk_scale(fit9_hk, dat9)
q13 <- hk_scale(fit13_hk, dat13)


# ------------------------------------------------------------------------------
# 5. Concise numerical comparison
# ------------------------------------------------------------------------------

get_summary <- function(fit_hk, fit_mhk, pred, d, qhk) {

  data.frame(
    studies = nrow(d),
    mu_REML = as.numeric(fit_hk$b[1]),
    tau2_REML = as.numeric(fit_hk$tau2),
    tau_REML = sqrt(as.numeric(fit_hk$tau2)),
    HK_scale = qhk,

    # Ordinary Hartung-Knapp
    HK_CI_lb = as.numeric(fit_hk$ci.lb),
    HK_CI_ub = as.numeric(fit_hk$ci.ub),

    # Modified Hartung-Knapp safeguard
    mHK_CI_lb = as.numeric(fit_mhk$ci.lb),
    mHK_CI_ub = as.numeric(fit_mhk$ci.ub),

    # Small-sample prediction interval based on the modified-HK fit
    Riley_PI_lb = as.numeric(pred$pi.lb),
    Riley_PI_ub = as.numeric(pred$pi.ub),

    Q = as.numeric(fit_hk$QE),
    Q_p = as.numeric(fit_hk$QEp),
    I2_percent = as.numeric(fit_hk$I2),
    H2 = as.numeric(fit_hk$H2)
  )
}

sum9  <- get_summary(fit9_hk, fit9_mhk, pred9, dat9, q9)
sum13 <- get_summary(fit13_hk, fit13_mhk, pred13, dat13, q13)

comparison <- rbind(
  cbind(set = "Nine-study Cochrane-scale subset", sum9),
  cbind(set = "Cochrane-13", sum13)
)

cat("\n============================================================\n")
cat("COCHRANE HISTORICAL COMPLETENESS SENSITIVITY ANALYSIS\n")
cat("============================================================\n\n")

print(comparison, digits = 5, row.names = FALSE)

cat("\nChange in pooled PF after adding four studies:\n")
cat(
  sprintf(
    "  %.6f - %.6f = %.6f (%.2f percentage points)\n",
    sum13$mu_REML,
    sum9$mu_REML,
    sum13$mu_REML - sum9$mu_REML,
    100 * (sum13$mu_REML - sum9$mu_REML)
  )
)

cat("\nChange in REML tau:\n")
cat(
  sprintf(
    "  %.6f - %.6f = %.6f\n",
    sum13$tau_REML,
    sum9$tau_REML,
    sum13$tau_REML - sum9$tau_REML
  )
)


# ------------------------------------------------------------------------------
# 6. Reproduction checks
# ------------------------------------------------------------------------------

# These are checks against the rounded manuscript values, not model inputs.

expected <- data.frame(
  set = c("Nine-study Cochrane-scale subset", "Cochrane-13"),
  mu  = c(0.439, 0.435),
  tau = c(0.217, 0.197)
)

obtained <- data.frame(
  set = c("Nine-study Cochrane-scale subset", "Cochrane-13"),
  mu = c(sum9$mu_REML, sum13$mu_REML),
  tau = c(sum9$tau_REML, sum13$tau_REML)
)

check <- merge(expected, obtained, by = "set", suffixes = c("_expected", "_obtained"))
check$mu_difference  <- check$mu_obtained  - check$mu_expected
check$tau_difference <- check$tau_obtained - check$tau_expected

cat("\n============================================================\n")
cat("CHECK AGAINST MANUSCRIPT VALUES\n")
cat("============================================================\n\n")

print(check, digits = 6, row.names = FALSE)

stopifnot(
  max(abs(check$mu_difference)) < 0.001,
  max(abs(check$tau_difference)) < 0.001
)


# ------------------------------------------------------------------------------
# 7. Forest plot for the 13-study historical set
# ------------------------------------------------------------------------------

# An asterisk marks the four studies that are not part of the original
# 9-study historical core.

study_label <- ifelse(
  dat13$core9,
  dat13$study,
  paste0(dat13$study, " *")
)


get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", args, value = TRUE)
  if (length(hit) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", hit[1]))))
  }
  normalizePath(getwd())
}
script_dir <- get_script_dir()

outdir <- file.path(script_dir, "output_04_Cochrane13")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

pdf(
  file.path(outdir, "Cochrane13_PF_forest.pdf"),
  width = 8.5,
  height = 8.0
)

forest(
  fit13_mhk,
  slab = study_label,
  xlab = "Prevented fraction (positive = prevention)",
  mlab = "REML--modified HK, 13 studies",
  alim = c(-1.0, 2.0),
  at = seq(-1.0, 2.0, by = 0.5),
  header = c("Study", "PF [95% CI]")
)

abline(v = 0, lty = 2)

mtext(
  "* Added in the Cochrane-13 completeness analysis",
  side = 1,
  line = 4.0,
  adj = 0
)

dev.off()


# ------------------------------------------------------------------------------
# 8. Direct visual comparison of the pooled 9- and 13-study estimates
# ------------------------------------------------------------------------------

pdf(
  file.path(outdir, "Cochrane9_vs_13_summary.pdf"),
  width = 7.5,
  height = 4.5
)

xlim <- range(
  c(
    sum9$mHK_CI_lb, sum9$mHK_CI_ub,
    sum13$mHK_CI_lb, sum13$mHK_CI_ub
  )
)

xlim <- xlim + c(-0.08, 0.08)

plot(
  NA,
  xlim = xlim,
  ylim = c(0.5, 2.5),
  yaxt = "n",
  ylab = "",
  xlab = "Prevented fraction",
  bty = "n"
)

abline(v = 0, lty = 2)

axis(
  2,
  at = c(2, 1),
  labels = c("Nine-study subset", "Cochrane-13"),
  las = 1,
  tick = FALSE
)

segments(sum9$mHK_CI_lb, 2, sum9$mHK_CI_ub, 2, lwd = 2)
points(sum9$mu_REML, 2, pch = 15, cex = 1.1)

segments(sum13$mHK_CI_lb, 1, sum13$mHK_CI_ub, 1, lwd = 2)
points(sum13$mu_REML, 1, pch = 15, cex = 1.1)

title("Historical completeness sensitivity: 9 versus 13 studies")

dev.off()


# ------------------------------------------------------------------------------
# 9. Save machine-readable output
# ------------------------------------------------------------------------------

write.csv(
  dat13,
  file = file.path(outdir, "Cochrane13_PF_data.csv"),
  row.names = FALSE
)

write.csv(
  comparison,
  file = file.path(outdir, "Cochrane9_vs_13_REML_HK_summary.csv"),
  row.names = FALSE
)

capture.output(
  fit9_hk,
  file = file.path(outdir, "fit9_HK.txt")
)

capture.output(
  fit13_hk,
  file = file.path(outdir, "fit13_HK.txt")
)

capture.output(
  sessionInfo(),
  file = file.path(outdir, "R_sessionInfo.txt")
)



# ------------------------------------------------------------------------------
# 10. Optional Bayesian completeness check
# ------------------------------------------------------------------------------

# The Supplement also reports a Bayesian 9-versus-13 comparison on this same
# legacy PF scale.  This is intentionally separate from the primary log-ratio
# Bayesian model.  If bayesmeta is installed, reproduce it here using a diffuse
# prior for the pooled PF and Half-Normal(0,0.5^2) for tau.

if (requireNamespace("bayesmeta", quietly = TRUE)) {

  tau_HN_050 <- function(tau) {
    bayesmeta::dhalfnormal(tau, scale = 0.50)
  }

  fit_bayes_pf <- function(d) {
    bayesmeta::bayesmeta(
      y = d$PF,
      sigma = d$SE,
      labels = d$study,
      mu.prior.mean = 0,
      mu.prior.sd = 10,
      tau.prior = tau_HN_050,
      interval.type = "central",
      delta = 0.005,
      epsilon = 1e-6
    )
  }

  summarize_bayes_pf <- function(fit, set_name) {
    probs <- c(0.025, 0.50, 0.975)
    q_mu <- unname(fit$qposterior(mu.p = probs))
    q_tau <- unname(fit$qposterior(tau.p = probs))
    q_pred <- unname(
      fit$qposterior(theta.p = probs, predict = TRUE)
    )

    data.frame(
      set = set_name,
      mu_lb = q_mu[1],
      mu_median = q_mu[2],
      mu_ub = q_mu[3],
      tau_lb = q_tau[1],
      tau_median = q_tau[2],
      tau_ub = q_tau[3],
      pred_lb = q_pred[1],
      pred_median = q_pred[2],
      pred_ub = q_pred[3],
      stringsAsFactors = FALSE
    )
  }

  fit9_bayes <- fit_bayes_pf(dat9)
  fit13_bayes <- fit_bayes_pf(dat13)

  bayes_comparison <- rbind(
    summarize_bayes_pf(fit9_bayes, "Nine-study Cochrane-scale subset"),
    summarize_bayes_pf(fit13_bayes, "Cochrane-13")
  )

  write.csv(
    bayes_comparison,
    file = file.path(outdir, "Cochrane9_vs_13_Bayesian_summary.csv"),
    row.names = FALSE
  )

  cat("\nBayesian Cochrane-scale completeness check:\n")
  print(bayes_comparison, digits = 4, row.names = FALSE)

} else {
  message(
    "Package 'bayesmeta' not installed: skipping optional Bayesian ",
    "Cochrane-13 completeness check."
  )
}


# ------------------------------------------------------------------------------
# 11. Final concise statement

# ------------------------------------------------------------------------------

cat("\n============================================================\n")
cat("FINAL HISTORICAL-COMPLETENESS RESULT\n")
cat("============================================================\n\n")

cat(
  sprintf(
    paste0(
      "Nine-study Cochrane-scale subset: PF = %.3f, tau = %.3f.\n",
      "Cochrane-13:                     PF = %.3f, tau = %.3f.\n",
      "Difference in pooled PF = %.3f (%.1f percentage points).\n"
    ),
    sum9$mu_REML,
    sum9$tau_REML,
    sum13$mu_REML,
    sum13$tau_REML,
    sum13$mu_REML - sum9$mu_REML,
    100 * (sum13$mu_REML - sum9$mu_REML)
  )
)

cat("\nOutput files written to: ", normalizePath(outdir), "\n", sep = "")

# ==============================================================================
# END
# ==============================================================================
