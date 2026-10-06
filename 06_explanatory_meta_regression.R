# ==============================================================================
# 06_explanatory_meta_regression.R
#
# Final 14-study Bayesian explanatory meta-regression.
#
# Models:
#   A. Time only (uses calendar_time_meta_regression.stan)
#   B. Time + organized background prevention
#   C. Time + organized background prevention + error-aware control burden
#
# The primary dataset is the same 14-study continuous-outcome set used in the
# final calendar-time model: historical 9 + Munoz-Millan + Latifi-Xhemajli +
# McMahon + Zeng + approximately cluster-adjusted Wang.
#
# Expected manuscript-scale posterior medians are approximately:
#   time only:                         beta1 = 0.091
#   time + organized prevention:      beta1 = 0.057
#   time + prevention + burden:       beta1 = 0.053
#
# These are explanatory study-level moderator analyses, not causal mediation.
# ==============================================================================

if (!requireNamespace("cmdstanr", quietly = TRUE)) {
  stop("Install 'cmdstanr' and CmdStan before running this script.")
}
if (!requireNamespace("posterior", quietly = TRUE)) {
  stop("Install 'posterior' before running this script.")
}

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", args, value = TRUE)
  if (length(hit) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", hit[1]))))
  }
  normalizePath(getwd())
}
script_dir <- get_script_dir()
outdir <- file.path(script_dir, "output_06_explanatory_regression")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

lrr_components <- function(mean_T, sd_T, n_T, mean_C, sd_C, n_C,
                           design_effect = 1) {
  xT <- log(mean_T)
  xC <- log(mean_C)

  vT <- sd_T^2 / (n_T * mean_T^2)
  vC <- sd_C^2 / (n_C * mean_C^2)

  # For Wang, both arm-level sampling variances are inflated by the same
  # publication-based cluster design effect.
  vT <- vT * design_effect
  vC <- vC * design_effect

  data.frame(
    xT = xT,
    xC = xC,
    vT = vT,
    vC = vC,
    yi = xT - xC,
    vi = vT + vC,
    sei = sqrt(vT + vC)
  )
}

dat <- data.frame(
  study = c(
    "Koch", "Modeer", "Clark", "Tewari", "Bravo", "Skold", "Arruda",
    "Tagliaferro", "Milsom", "Munoz-Millan", "Latifi-Xhemajli",
    "McMahon", "Zeng", "Wang"
  ),
  year = c(
    1975, 1984, 1985, 1990, 1997, 2005, 2012,
    2011, 2011, 2018, 2019, 2020, 2025, 2022
  ),
  mean_T = c(
    0.90, 2.50, 2.43, 0.55, 1.48, 0.79, 4.61,
    0.33, 0.66, 1.60, 5.20, 3.50, 2.05, 0.38
  ),
  sd_T = c(
    3.80, 3.10, 3.09, 4.48, 1.53, 1.67, 4.39,
    1.04, 0.73, 2.40, 10.50, 5.90, 2.64, 1.21
  ),
  n_T = c(
    60, 87, 246, 311, 98, 190, 113,
    91, 94, 131, 218, 577, 79, 2385
  ),
  mean_C = c(
    4.00, 3.70, 3.11, 2.16, 2.58, 1.85, 7.72,
    0.57, 0.63, 2.10, 10.10, 3.50, 1.96, 0.61
  ),
  sd_C = c(
    3.75, 3.90, 3.54, 4.12, 1.89, 2.89, 5.84,
    1.39, 0.66, 2.50, 12.90, 4.90, 2.56, 1.60
  ),
  n_C = c(
    61, 107, 234, 307, 116, 181, 97,
    86, 95, 144, 209, 573, 80, 2620
  ),
  B = c(
    0, 0, 0, 0, 0, 0, 0,
    1, 0, 1, 0, 1, 1, 0
  ),
  stringsAsFactors = FALSE
)

# Wang: ICC=.20 and planning cluster size m=40 -> DE=8.8.
dat$design_effect <- 1
dat$design_effect[dat$study == "Wang"] <- 1 + (40 - 1) * 0.20

parts <- do.call(
  rbind,
  lapply(seq_len(nrow(dat)), function(i) {
    lrr_components(
      dat$mean_T[i], dat$sd_T[i], dat$n_T[i],
      dat$mean_C[i], dat$sd_C[i], dat$n_C[i],
      dat$design_effect[i]
    )
  })
)

dat <- cbind(dat, parts)
dat$tdec <- (dat$year - 2000) / 10

write.csv(
  dat,
  file.path(outdir, "explanatory_14study_input.csv"),
  row.names = FALSE
)

summarize_fit <- function(fit, model_name, vars) {
  d <- posterior::as_draws_df(fit$draws(variables = vars))

  rows <- lapply(vars, function(v) {
    q <- unname(quantile(d[[v]], c(0.025, 0.50, 0.975)))
    data.frame(
      model = model_name,
      parameter = v,
      q025 = q[1],
      median = q[2],
      q975 = q[3],
      Pr_gt_0 = mean(d[[v]] > 0),
      Pr_lt_0 = mean(d[[v]] < 0),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

sample_model <- function(stan_file, data_list, seed) {
  mod <- cmdstanr::cmdstan_model(stan_file)
  mod$sample(
    data = data_list,
    seed = seed,
    chains = 4,
    parallel_chains = min(4, parallel::detectCores()),
    iter_warmup = 2000,
    iter_sampling = 4000,
    adapt_delta = 0.995,
    max_treedepth = 15,
    refresh = 250
  )
}

# ------------------------------------------------------------------------------
# A. Time only
# ------------------------------------------------------------------------------
time_file <- file.path(script_dir, "calendar_time_meta_regression.stan")
if (!file.exists(time_file)) stop("Missing ", time_file)

fit_time <- sample_model(
  time_file,
  list(
    N = nrow(dat),
    y = dat$yi,
    se = dat$sei,
    tdec = dat$tdec
  ),
  seed = 20260931
)

# ------------------------------------------------------------------------------
# B. Time + organized background prevention
# ------------------------------------------------------------------------------
prevention_file <- file.path(script_dir, "time_background_prevention.stan")
if (!file.exists(prevention_file)) stop("Missing ", prevention_file)

fit_prevention <- sample_model(
  prevention_file,
  list(
    N = nrow(dat),
    y = dat$yi,
    se = dat$sei,
    tdec = dat$tdec,
    B = dat$B
  ),
  seed = 20260932
)

# ------------------------------------------------------------------------------
# C. Time + prevention + error-aware control burden
# ------------------------------------------------------------------------------
burden_file <- file.path(script_dir, "time_background_burden.stan")
if (!file.exists(burden_file)) stop("Missing ", burden_file)

xCbar <- mean(dat$xC)
s_logC <- sd(dat$xC)

fit_burden <- sample_model(
  burden_file,
  list(
    N = nrow(dat),
    xT = dat$xT,
    xC = dat$xC,
    vT = dat$vT,
    vC = dat$vC,
    tdec = dat$tdec,
    B = dat$B,
    xCbar = xCbar,
    s_logC = s_logC
  ),
  seed = 20260933
)

results <- rbind(
  summarize_fit(fit_time, "Time only", c("beta0", "beta1", "tau")),
  summarize_fit(
    fit_prevention,
    "Time + organized prevention",
    c("beta0", "beta1", "beta2", "tau")
  ),
  summarize_fit(
    fit_burden,
    "Time + organized prevention + burden",
    c("beta0", "beta1", "beta2", "beta3", "tau")
  )
)

write.csv(
  results,
  file.path(outdir, "explanatory_meta_regression_posterior_summary.csv"),
  row.names = FALSE
)

capture.output(
  fit_time$diagnostic_summary(),
  file = file.path(outdir, "time_only_diagnostics.txt")
)
capture.output(
  fit_prevention$diagnostic_summary(),
  file = file.path(outdir, "time_prevention_diagnostics.txt")
)
capture.output(
  fit_burden$diagnostic_summary(),
  file = file.path(outdir, "time_prevention_burden_diagnostics.txt")
)
capture.output(
  sessionInfo(),
  file = file.path(outdir, "R_sessionInfo.txt")
)

cat("\n============================================================\n")
cat("FINAL 14-STUDY EXPLANATORY META-REGRESSION\n")
cat("============================================================\n\n")
print(results, row.names = FALSE, digits = 4)

cat(
  "\nExpected manuscript-scale beta1 medians: ",
  "0.091 (time), 0.057 (time+prevention), ",
  "0.053 (time+prevention+burden).\n",
  sep = ""
)
cat("\nOutputs written to: ", normalizePath(outdir), "\n", sep = "")
