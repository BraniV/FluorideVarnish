# ==============================================================================
# 05_temporal_meta_regression.R
#
# Fluoride-varnish continuous-outcome production / sensitivity workflow
#
# PURPOSE
# -------
# This is the final continuous-outcome production/sensitivity workflow used
# for the historical-to-contemporary calendar-time analysis.
#
# It:
#   1. rebuilds the 14-study primary historical + contemporary dataset;
#   2. uses primary-source arm summaries where available;
#   3. adds Wang (2022) with a publication-based cluster design-effect
#      approximation and a second cluster-size sensitivity;
#   4. treats Jayasinghe (2022) as sensitivity-only because the paper reports
#      means but no SD/SE for the 2-year new-caries count;
#   5. treats Ghasemi (2025) as sensitivity-only because only four schools
#      were randomized and no outcome ICC or school-level summary is reported;
#   6. runs REML--Hartung--Knapp calendar-time meta-regression;
#   7. runs historical-vs-contemporary meta-regression;
#   8. pools the contemporary studies separately;
#   9. compares time functional forms by ML AICc;
#  10. runs leave-one-contemporary-study-out analyses;
#  11. optionally runs the Bayesian calendar-time model with CmdStanR;
#  12. generates temporal figures and CSV/LaTeX output.
#
# IMPORTANT INTERPRETATION
# ------------------------
# The strongest publication-based expanded primary analysis is "core14_wang":
# the nine historical studies + the four already usable modern studies + Wang.
#
# Jayasinghe and Ghasemi cannot be assigned an exact defensible LRR variance
# from the published aggregate summaries alone.  They are therefore included
# only in explicit variance-sensitivity analyses.
#
# ==============================================================================
#
# Required:
#   install.packages("metafor")
#
# Optional Bayesian production model:
#   cmdstanr, posterior, and a working CmdStan installation.
#
# ==============================================================================

if (!requireNamespace("metafor", quietly = TRUE)) {
  stop("Install 'metafor' first: install.packages('metafor')")
}

library(metafor)

RUN_STAN <- TRUE

if (RUN_STAN) {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop("RUN_STAN=TRUE but package 'cmdstanr' is not installed.")
  }
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("RUN_STAN=TRUE but package 'posterior' is not installed.")
  }
}

# ------------------------------------------------------------------------------
# 0. Output directory
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
outdir <- file.path(script_dir, "output_05_temporal_regression")

if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# 1. Helper functions
# ------------------------------------------------------------------------------

lrr_from_arms <- function(mean_T, sd_T, n_T, mean_C, sd_C, n_C) {

  stopifnot(
    mean_T > 0, mean_C > 0,
    sd_T >= 0, sd_C >= 0,
    n_T > 1, n_C > 1
  )

  yi <- log(mean_T / mean_C)

  vi <-
    sd_T^2 / (n_T * mean_T^2) +
    sd_C^2 / (n_C * mean_C^2)

  c(yi = yi, vi = vi, sei = sqrt(vi))
}


make_row_from_arms <- function(study, year, era,
                               mean_T, sd_T, n_T,
                               mean_C, sd_C, n_C,
                               endpoint,
                               variance_source,
                               variance_status = "direct delta-method") {

  z <- lrr_from_arms(
    mean_T, sd_T, n_T,
    mean_C, sd_C, n_C
  )

  data.frame(
    study = study,
    year = year,
    era = era,

    mean_T = mean_T,
    sd_T = sd_T,
    n_T = n_T,

    mean_C = mean_C,
    sd_C = sd_C,
    n_C = n_C,

    yi = unname(z["yi"]),
    vi = unname(z["vi"]),
    sei = unname(z["sei"]),

    endpoint = endpoint,
    variance_source = variance_source,
    variance_status = variance_status,

    stringsAsFactors = FALSE
  )
}


fit_time_reml_hk <- function(dat) {

  d <- dat
  d$tdec <- (d$year - 2000) / 10

  fit <- metafor::rma.uni(
    yi = yi,
    vi = vi,
    mods = ~ tdec,
    data = d,
    method = "REML",
    test = "knha"
  )

  fit
}


fit_era_reml_hk <- function(dat) {

  d <- dat
  d$modern <- as.integer(d$era == "contemporary")

  metafor::rma.uni(
    yi = yi,
    vi = vi,
    mods = ~ modern,
    data = d,
    method = "REML",
    test = "knha"
  )
}


fit_modern_pool <- function(dat) {

  d <- subset(dat, era == "contemporary")

  metafor::rma.uni(
    yi = yi,
    vi = vi,
    data = d,
    method = "REML",
    test = "knha"
  )
}


extract_time_result <- function(fit, scenario, k) {

  data.frame(
    scenario = scenario,
    k = k,

    beta0 = unname(coef(fit)[1]),
    beta1 = unname(coef(fit)[2]),

    beta1_se_HK = fit$se[2],
    beta1_ci_lb = fit$ci.lb[2],
    beta1_ci_ub = fit$ci.ub[2],

    tau = sqrt(fit$tau2),

    QE = fit$QE,
    QEp = fit$QEp,

    stringsAsFactors = FALSE
  )
}


extract_era_result <- function(fit, scenario, k) {

  data.frame(
    scenario = scenario,
    k = k,

    beta_hist = unname(coef(fit)[1]),
    beta_modern_minus_hist = unname(coef(fit)[2]),

    contrast_se_HK = fit$se[2],
    contrast_ci_lb = fit$ci.lb[2],
    contrast_ci_ub = fit$ci.ub[2],

    tau = sqrt(fit$tau2),

    stringsAsFactors = FALSE
  )
}


extract_modern_pool <- function(fit, scenario, k_modern) {

  mu <- as.numeric(coef(fit)[1])

  # PF = 1-exp(mu) is decreasing in mu; interval endpoints reverse.
  pf <- 1 - exp(mu)
  pf_lb <- 1 - exp(fit$ci.ub)
  pf_ub <- 1 - exp(fit$ci.lb)

  data.frame(
    scenario = scenario,
    k_modern = k_modern,

    mu_modern = mu,
    mu_ci_lb = fit$ci.lb,
    mu_ci_ub = fit$ci.ub,

    PF_modern = pf,
    PF_ci_lb = pf_lb,
    PF_ci_ub = pf_ub,

    tau_modern = sqrt(fit$tau2),

    stringsAsFactors = FALSE
  )
}


aicc_from_ml_fit <- function(fit) {

  ll <- as.numeric(logLik(fit))

  # Fixed-effect coefficients + tau^2.
  p <- length(coef(fit)) + 1
  k <- fit$k

  aic <- -2 * ll + 2 * p

  if (k <= p + 1) {
    return(NA_real_)
  }

  aic + 2 * p * (p + 1) / (k - p - 1)
}


model_comparison_ml <- function(dat, scenario) {

  d <- dat
  d$tdec <- (d$year - 2000) / 10
  d$modern <- as.integer(d$era == "contemporary")

  d$era3 <- cut(
    d$year,
    breaks = c(-Inf, 1999, 2017, Inf),
    labels = c("pre2000", "2000_2017", "2018plus"),
    right = TRUE
  )

  fits <- list(
    no_calendar = metafor::rma.uni(
      yi, vi, data = d,
      method = "ML"
    ),

    linear_time = metafor::rma.uni(
      yi, vi, mods = ~ tdec, data = d,
      method = "ML"
    ),

    historical_vs_contemporary = metafor::rma.uni(
      yi, vi, mods = ~ modern, data = d,
      method = "ML"
    ),

    quadratic_time = metafor::rma.uni(
      yi, vi, mods = ~ tdec + I(tdec^2), data = d,
      method = "ML"
    ),

    three_eras = metafor::rma.uni(
      yi, vi, mods = ~ era3, data = d,
      method = "ML"
    )
  )

  out <- data.frame(
    scenario = scenario,
    model = names(fits),
    AICc = vapply(fits, aicc_from_ml_fit, numeric(1)),
    stringsAsFactors = FALSE
  )

  out$delta_AICc <- out$AICc - min(out$AICc, na.rm = TRUE)

  w <- exp(-0.5 * out$delta_AICc)
  out$Akaike_weight <- w / sum(w)

  out[order(out$AICc), ]
}


leave_one_modern_out <- function(dat, scenario) {

  modern_studies <- dat$study[dat$era == "contemporary"]

  ans <- do.call(
    rbind,
    lapply(
      modern_studies,
      function(s) {

        d <- dat[dat$study != s, , drop = FALSE]

        f <- fit_time_reml_hk(d)

        data.frame(
          scenario = scenario,
          omitted = s,
          k = nrow(d),

          beta1 = unname(coef(f)[2]),
          beta1_se_HK = f$se[2],
          beta1_ci_lb = f$ci.lb[2],
          beta1_ci_ub = f$ci.ub[2],

          tau = sqrt(f$tau2),

          stringsAsFactors = FALSE
        )
      }
    )
  )

  ans
}


# ------------------------------------------------------------------------------
# 2. Historical nine-study dataset
# ------------------------------------------------------------------------------

historical <- data.frame(
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

  year = c(
    1975, 1984, 1985, 1990, 1997,
    2005, 2012, 2011, 2011
  ),

  mean_T = c(
    0.90, 2.50, 2.43, 0.55, 1.48,
    0.79, 4.61, 0.33, 0.66
  ),

  sd_T = c(
    3.80, 3.10, 3.09, 4.48, 1.53,
    1.67, 4.39, 1.04, 0.73
  ),

  n_T = c(
    60, 87, 246, 311, 98,
    190, 113, 91, 94
  ),

  mean_C = c(
    4.00, 3.70, 3.11, 2.16, 2.58,
    1.85, 7.72, 0.57, 0.63
  ),

  sd_C = c(
    3.75, 3.90, 3.54, 4.12, 1.89,
    2.89, 5.84, 1.39, 0.66
  ),

  n_C = c(
    61, 107, 234, 307, 116,
    181, 97, 86, 95
  ),

  stringsAsFactors = FALSE
)

historical$era <- "historical"
historical$endpoint <- "Historical DMF-type continuous outcome"
historical$variance_source <- "Arm means/SD/n; first-order LRR delta method"
historical$variance_status <- "historical core"

tmp <- mapply(
  FUN = function(mt, st, nt, mc, sc, nc) {
    lrr_from_arms(mt, st, nt, mc, sc, nc)
  },
  historical$mean_T,
  historical$sd_T,
  historical$n_T,
  historical$mean_C,
  historical$sd_C,
  historical$n_C
)

historical$yi <- tmp["yi", ]
historical$vi <- tmp["vi", ]
historical$sei <- tmp["sei", ]


# ------------------------------------------------------------------------------
# 3. Four contemporary studies already suitable for the weighted analysis
# ------------------------------------------------------------------------------

# IMPORTANT UPDATE FOR MUNOZ-MILLAN:
# The primary article reports SDs 2.4 and 2.5 for the 24-month dmft means
# 1.6 and 2.1.  Some secondary extractions report 2.0 and 2.6.
# This production-audit script uses the PRIMARY-ARTICLE values 2.4 and 2.5.

munoz <- make_row_from_arms(
  study = "Munoz-Millan",
  year = 2018,
  era = "contemporary",

  mean_T = 1.6,
  sd_T = 2.4,
  n_T = 131,

  mean_C = 2.1,
  sd_C = 2.5,
  n_C = 144,

  endpoint = "24-month dmft endpoint burden",
  variance_source = "Primary article: mean(SD) 1.6(2.4) vs 2.1(2.5)"
)


latifi <- make_row_from_arms(
  study = "Latifi-Xhemajli",
  year = 2019,
  era = "contemporary",

  mean_T = 5.2,
  sd_T = 10.5,
  n_T = 218,

  mean_C = 10.1,
  sd_C = 12.9,
  n_C = 209,

  endpoint = "24-month dmfs endpoint burden",
  variance_source = "Primary article: endpoint mean(SD) 5.2(10.5) vs 10.1(12.9)"
)


mcmahon <- make_row_from_arms(
  study = "McMahon",
  year = 2020,
  era = "contemporary",

  mean_T = 3.5,
  sd_T = 5.9,
  n_T = 577,

  mean_C = 3.5,
  sd_C = 4.9,
  n_C = 573,

  endpoint = "24-month d3mfs endpoint burden",
  variance_source = paste(
    "Primary trial reports final means 3.5 vs 3.5;",
    "SDs 5.9 and 4.9 from USPSTF evidence extraction"
  )
)


# Zeng Table 5: TFV+OHE versus OHE at 24 months.
# Table 5 gives caries incidence counts 31/79 and 29/80, which identify
# the corresponding analyzed group denominators.

zeng <- make_row_from_arms(
  study = "Zeng",
  year = 2025,
  era = "contemporary",

  mean_T = 2.05,
  sd_T = 2.64,
  n_T = 79,

  mean_C = 1.96,
  sd_C = 2.56,
  n_C = 80,

  endpoint = "24-month full-mouth Delta DMFT/dmft",
  variance_source = "Primary article Table 5; TFV+OHE versus OHE"
)


core_modern4 <- rbind(
  munoz,
  latifi,
  mcmahon,
  zeng
)

core13 <- rbind(
  historical,
  core_modern4
)


# ------------------------------------------------------------------------------
# 4. Wang 2022: publication-based cluster adjustment
# ------------------------------------------------------------------------------

# Published 24-month first-permanent-molar DFS increment:
#     test    0.38 (SD 1.21), n=2385
#     control 0.61 (SD 1.60), n=2620
#
# 107 classes: 51 test, 56 control.
#
# The article states an ICC of 0.20 and expected m=40 children/class in its
# sample-size discussion.  It does NOT provide a directly cluster-robust
# standard error for the log ratio.
#
# We therefore:
#   (a) calculate the ordinary individual-level LRR variance;
#   (b) multiply it by DE = 1 + (m-1)*ICC;
#   (c) use m=40, ICC=.20 as the publication-based approximation;
#   (d) repeat with actual average completed cluster size 5005/107 as a
#       cluster-size sensitivity.
#
# This should be described as "approximately cluster-adjusted", not as an
# exact variance from a published cluster model.

wang_raw <- make_row_from_arms(
  study = "Wang",
  year = 2022,
  era = "contemporary",

  mean_T = 0.38,
  sd_T = 1.21,
  n_T = 2385,

  mean_C = 0.61,
  sd_C = 1.60,
  n_C = 2620,

  endpoint = "24-month DFS increment, first permanent molars",
  variance_source = "Primary article arm mean(SD)/n",
  variance_status = "naive individual-level; requires cluster inflation"
)

wang_ICC <- 0.20
wang_m_planning <- 40
wang_DE_planning <- 1 + (wang_m_planning - 1) * wang_ICC

wang_primary <- wang_raw
wang_primary$vi <- wang_raw$vi * wang_DE_planning
wang_primary$sei <- sqrt(wang_primary$vi)
wang_primary$variance_source <- paste0(
  "Primary article LRR variance x DE; ICC=0.20, m=40, DE=",
  sprintf("%.3f", wang_DE_planning)
)
wang_primary$variance_status <- "publication-based approximate cluster adjustment"

wang_m_actual <- 5005 / 107
wang_DE_actual <- 1 + (wang_m_actual - 1) * wang_ICC

wang_actual_m <- wang_raw
wang_actual_m$vi <- wang_raw$vi * wang_DE_actual
wang_actual_m$sei <- sqrt(wang_actual_m$vi)
wang_actual_m$variance_source <- paste0(
  "Sensitivity: ICC=0.20, actual mean completed cluster size=5005/107, DE=",
  sprintf("%.3f", wang_DE_actual)
)
wang_actual_m$variance_status <- "cluster-size sensitivity"


# ------------------------------------------------------------------------------
# 5. Jayasinghe 2022: variance cannot be recovered exactly
# ------------------------------------------------------------------------------

# Published 2-year mean new caries:
#     intervention = 1.50
#     control      = 1.97
#
# Table 3 reports the means and a Mann-Whitney p-value but NO SD or SE.
# Therefore an ordinary LRR sampling variance cannot be uniquely recovered.
#
# For sensitivity only, treat the count variance as:
#     Var(Y) = phi * E(Y)
#
# which gives the approximate log-mean variance:
#     Var[log(mean)] ~= phi / (n * mean).
#
# phi=1 is equidispersion; larger phi permits overdispersion.
#
# The published table labels Control N=161 and Intervention N=162.
# We use those table values here.

jay_mean_T <- 1.50
jay_n_T <- 162

jay_mean_C <- 1.97
jay_n_C <- 161

jay_yi <- log(jay_mean_T / jay_mean_C)

jay_make <- function(phi) {

  vi <-
    phi / (jay_n_T * jay_mean_T) +
    phi / (jay_n_C * jay_mean_C)

  data.frame(
    study = "Jayasinghe",
    year = 2022,
    era = "contemporary",

    mean_T = jay_mean_T,
    sd_T = NA_real_,
    n_T = jay_n_T,

    mean_C = jay_mean_C,
    sd_C = NA_real_,
    n_C = jay_n_C,

    yi = jay_yi,
    vi = vi,
    sei = sqrt(vi),

    endpoint = "2-year mean number of new caries",
    variance_source = paste0(
      "SENSITIVITY ONLY: quasi-Poisson count variance, phi=", phi
    ),
    variance_status = "assumption-based sensitivity; published SD unavailable",

    stringsAsFactors = FALSE
  )
}

jay_phi_grid <- c(1, 1.5, 2, 3, 4)

jay_rows <- lapply(
  jay_phi_grid,
  jay_make
)

names(jay_rows) <- paste0("phi_", jay_phi_grid)


# ------------------------------------------------------------------------------
# 6. Ghasemi 2025: only four randomized schools; no outcome ICC reported
# ------------------------------------------------------------------------------

# Published 6-month DMFS increment:
#     intervention: 0.2 (SD .77), n=173
#     control:      0.8 (SD 1.59), n=179
#
# Four schools were randomized: TWO schools per arm.
#
# The paper reports individual-level t tests for the continuous outcomes and
# does not report an outcome ICC or school-specific DMFS increments.  Therefore
# an exact cluster-aware LRR variance cannot be reconstructed from the paper.
#
# We first calculate the naive individual-level variance, then use an ICC
# sensitivity grid:
#
#     DE = 1 + (m-1)*ICC
#
# with m approximately 88 analyzed children per school.
#
# With only two clusters/arm, this DE sensitivity should NOT be mistaken for
# a full cluster-randomized-trial reanalysis.  It is only an influence /
# uncertainty sensitivity analysis.

ghasemi_raw <- make_row_from_arms(
  study = "Ghasemi",
  year = 2025,
  era = "contemporary",

  mean_T = 0.2,
  sd_T = 0.77,
  n_T = 173,

  mean_C = 0.8,
  sd_C = 1.59,
  n_C = 179,

  endpoint = "6-month DMFS increment",
  variance_source = "Primary article Table 3 arm mean(SD)/n",
  variance_status = "naive individual-level; only 2 randomized schools per arm"
)

ghasemi_m <- (173 + 179) / 4
ghasemi_rho_grid <- c(0.005, 0.01, 0.02, 0.05, 0.10, 0.20)

ghasemi_make <- function(rho) {

  DE <- 1 + (ghasemi_m - 1) * rho

  z <- ghasemi_raw
  z$vi <- ghasemi_raw$vi * DE
  z$sei <- sqrt(z$vi)

  z$variance_source <- paste0(
    "SENSITIVITY ONLY: DE with rho=", rho,
    ", m=", sprintf("%.1f", ghasemi_m),
    ", DE=", sprintf("%.3f", DE)
  )

  z$variance_status <-
    "assumption-based cluster sensitivity; 2 schools/arm; no published ICC"

  z
}

ghasemi_rows <- lapply(
  ghasemi_rho_grid,
  ghasemi_make
)

names(ghasemi_rows) <- paste0("rho_", ghasemi_rho_grid)


# ------------------------------------------------------------------------------
# 7. Freeze analysis scenarios
# ------------------------------------------------------------------------------

# Primary publication-based datasets:
#
# core13:
#   historical 9 + Munoz + Latifi + McMahon + Zeng.
#
# core14_wang:
#   core13 + Wang with the article's ICC=.20 and planning m=40.
#
# Jayasinghe and Ghasemi are deliberately NOT in the primary weighted dataset
# because their exact sampling variances are not recoverable from the paper.

scenarios <- list(

  core13 = core13,

  core14_wang = rbind(
    core13,
    wang_primary
  ),

  core14_wang_actual_cluster_size = rbind(
    core13,
    wang_actual_m
  )
)

# Representative Jayasinghe sensitivity, phi=2.
scenarios$core15_wang_jay_phi2 <- rbind(
  core13,
  wang_primary,
  jay_rows[["phi_2"]]
)

# Full Ghasemi ICC grid with Jayasinghe phi=2.
for (rho in ghasemi_rho_grid) {

  nm <- paste0(
    "all16_jayphi2_ghasrho",
    format(rho, trim = TRUE, scientific = FALSE)
  )

  scenarios[[nm]] <- rbind(
    core13,
    wang_primary,
    jay_rows[["phi_2"]],
    ghasemi_make(rho)
  )
}


# ------------------------------------------------------------------------------
# 8. Save frozen study datasets
# ------------------------------------------------------------------------------

for (nm in names(scenarios)) {

  write.csv(
    scenarios[[nm]],
    file.path(outdir, paste0("data_", nm, ".csv")),
    row.names = FALSE
  )
}


# ------------------------------------------------------------------------------
# 9. REML--HK calendar-time, era contrast, and modern pooling
# ------------------------------------------------------------------------------

time_results <- list()
era_results <- list()
modern_results <- list()

for (nm in names(scenarios)) {

  d <- scenarios[[nm]]

  f_time <- fit_time_reml_hk(d)
  f_era <- fit_era_reml_hk(d)
  f_mod <- fit_modern_pool(d)

  time_results[[nm]] <-
    extract_time_result(
      f_time,
      nm,
      nrow(d)
    )

  era_results[[nm]] <-
    extract_era_result(
      f_era,
      nm,
      nrow(d)
    )

  modern_results[[nm]] <-
    extract_modern_pool(
      f_mod,
      nm,
      sum(d$era == "contemporary")
    )
}

time_results <- do.call(rbind, time_results)
era_results <- do.call(rbind, era_results)
modern_results <- do.call(rbind, modern_results)

row.names(time_results) <- NULL
row.names(era_results) <- NULL
row.names(modern_results) <- NULL

write.csv(
  time_results,
  file.path(outdir, "REML_HK_calendar_time_scenarios.csv"),
  row.names = FALSE
)

write.csv(
  era_results,
  file.path(outdir, "REML_HK_historical_vs_contemporary_scenarios.csv"),
  row.names = FALSE
)

write.csv(
  modern_results,
  file.path(outdir, "REML_HK_contemporary_pool_scenarios.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 10. Full Jayasinghe x Ghasemi variance-sensitivity grid
# ------------------------------------------------------------------------------

grid_results <- list()

counter <- 1

for (phi in jay_phi_grid) {

  for (rho in ghasemi_rho_grid) {

    d <- rbind(
      core13,
      wang_primary,
      jay_make(phi),
      ghasemi_make(rho)
    )

    ft <- fit_time_reml_hk(d)
    fm <- fit_modern_pool(d)

    grid_results[[counter]] <- data.frame(
      jay_phi = phi,
      ghasemi_rho = rho,

      beta1 = unname(coef(ft)[2]),
      beta1_se_HK = ft$se[2],
      beta1_ci_lb = ft$ci.lb[2],
      beta1_ci_ub = ft$ci.ub[2],
      tau_time = sqrt(ft$tau2),

      mu_modern = unname(coef(fm)[1]),
      PF_modern = 1 - exp(unname(coef(fm)[1])),
      tau_modern = sqrt(fm$tau2),

      stringsAsFactors = FALSE
    )

    counter <- counter + 1
  }
}

grid_results <- do.call(rbind, grid_results)

write.csv(
  grid_results,
  file.path(outdir, "Jayasinghe_Ghasemi_variance_sensitivity_grid.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 11. ML model comparison by AICc
# ------------------------------------------------------------------------------

# Run model comparison for the original core13 and the strongest
# publication-based expanded analysis core14_wang.

aicc <- rbind(
  model_comparison_ml(
    scenarios$core13,
    "core13"
  ),

  model_comparison_ml(
    scenarios$core14_wang,
    "core14_wang"
  )
)

write.csv(
  aicc,
  file.path(outdir, "ML_AICc_model_comparison.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 12. Leave-one-contemporary-study-out analysis
# ------------------------------------------------------------------------------

loo <- rbind(

  leave_one_modern_out(
    scenarios$core13,
    "core13"
  ),

  leave_one_modern_out(
    scenarios$core14_wang,
    "core14_wang"
  ),

  leave_one_modern_out(
    scenarios$all16_jayphi2_ghasrho0.05,
    "all16_jayphi2_ghasrho0.05"
  )
)

write.csv(
  loo,
  file.path(outdir, "leave_one_contemporary_out.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 13. Audit table
# ------------------------------------------------------------------------------

audit <- data.frame(
  study = c(
    "Munoz-Millan",
    "Latifi-Xhemajli",
    "McMahon",
    "Zeng",
    "Wang",
    "Jayasinghe",
    "Ghasemi"
  ),

  primary_weighted_status = c(
    "YES",
    "YES",
    "YES, endpoint-burden sensitivity",
    "YES",
    "YES, approximately cluster-adjusted",
    "NO: sensitivity only",
    "NO: sensitivity only"
  ),

  key_variance_issue = c(
    "Primary article SDs are 2.4 and 2.5; secondary USPSTF table differs",
    "Arm SDs and analyzed n are available",
    "Continuous final d3mfs SDs taken from USPSTF evidence extraction",
    "Arm SDs are available; analyzed denominators obtained from primary table",
    "Class randomized; article gives ICC=.20 and planning m=40 but no LRR cluster SE",
    "2-year means reported without SD/SE; Mann-Whitney test used",
    "Only four randomized schools (2/arm); no outcome ICC or school-level increment means"
  ),

  action_taken = c(
    "Use primary article 1.6(2.4), n=131 vs 2.1(2.5), n=144",
    "Direct first-order LRR delta variance",
    "Direct first-order LRR delta variance using extracted SDs",
    "Direct first-order LRR delta variance",
    "Inflate naive LRR variance by DE=1+(40-1)*.20=8.8; actual-m sensitivity also run",
    "Quasi-Poisson overdispersion grid phi=1,1.5,2,3,4",
    "Design-effect ICC grid rho=.005,.01,.02,.05,.10,.20 with m=88"
  ),

  recommendation = c(
    "Retain",
    "Retain, label endpoint rather than pure increment",
    "Retain, label endpoint rather than pure increment",
    "Retain",
    "Retain in expanded analysis, describe cluster variance as approximate",
    "Do not use in primary weighted model without additional variance information",
    "Do not use in primary weighted model without cluster-level information"
  ),

  stringsAsFactors = FALSE
)

write.csv(
  audit,
  file.path(outdir, "contemporary_variance_audit.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 14. Temporal figure
# ------------------------------------------------------------------------------

make_temporal_figure <- function(dat, scenario, filename) {

  d <- dat
  d$tdec <- (d$year - 2000) / 10

  fit <- fit_time_reml_hk(d)

  year_grid <- seq(
    min(d$year) - 1,
    max(d$year) + 1,
    length.out = 250
  )

  t_grid <- (year_grid - 2000) / 10

  pr <- predict(
    fit,
    newmods = matrix(t_grid, ncol = 1)
  )

  yrange <- range(
    c(
      d$yi - 1.96 * d$sei,
      d$yi + 1.96 * d$sei,
      pr$ci.lb,
      pr$ci.ub
    ),
    finite = TRUE
  )

  pdf(
    file.path(outdir, filename),
    width = 9.0,
    height = 6.5
  )

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)

  par(mar = c(4.5, 4.8, 3.0, 5.0))

  plot(
    d$year,
    d$yi,
    type = "n",
    ylim = yrange,
    xlab = "Publication year",
    ylab = "Log ratio of means",
    main = paste0("Calendar-time meta-regression: ", scenario)
  )

  polygon(
    c(year_grid, rev(year_grid)),
    c(pr$ci.lb, rev(pr$ci.ub)),
    border = NA,
    density = 18,
    angle = 45
  )

  lines(
    year_grid,
    pr$pred,
    lwd = 2
  )

  abline(
    h = 0,
    lty = 3
  )

  # Study-level 95% sampling intervals.
  segments(
    d$year,
    d$yi - 1.96 * d$sei,
    d$year,
    d$yi + 1.96 * d$sei,
    lwd = 0.8
  )

  is_mod <- d$era == "contemporary"

  points(
    d$year[!is_mod],
    d$yi[!is_mod],
    pch = 1,
    cex = 1.05
  )

  points(
    d$year[is_mod],
    d$yi[is_mod],
    pch = 17,
    cex = 1.05
  )

  # Label modern studies only; historical labels would overcrowd the figure.
  text(
    d$year[is_mod],
    d$yi[is_mod],
    labels = d$study[is_mod],
    pos = 4,
    cex = 0.68,
    offset = 0.35
  )

  legend(
    "bottomleft",
    legend = c(
      "Historical",
      "Contemporary",
      "REML-HK fitted mean"
    ),
    pch = c(1, 17, NA),
    lty = c(NA, NA, 1),
    lwd = c(NA, NA, 2),
    bty = "n"
  )

  # Secondary right axis: selected PF values mapped to the log-ratio scale.
  # PF = 1-exp(y), so y=log(1-PF).
  pf_ticks <- c(-0.50, 0, 0.20, 0.40, 0.60, 0.80)
  y_ticks <- log(1 - pf_ticks)

  keep <- y_ticks >= yrange[1] & y_ticks <= yrange[2]

  axis(
    4,
    at = y_ticks[keep],
    labels = paste0(
      sprintf("%.0f", 100 * pf_ticks[keep]),
      "%"
    )
  )

  mtext(
    expression(paste("Prevented fraction  ", PF["+"] == 1-exp(y))),
    side = 4,
    line = 3
  )

  dev.off()
}


make_temporal_figure(
  scenarios$core13,
  "core13",
  "temporal_core13.pdf"
)

make_temporal_figure(
  scenarios$core14_wang,
  "core14_wang",
  "temporal_core14_wang.pdf"
)

make_temporal_figure(
  scenarios$all16_jayphi2_ghasrho0.05,
  "all16 sensitivity: Jay phi=2, Ghasemi rho=.05",
  "temporal_all16_sensitivity.pdf"
)


# ------------------------------------------------------------------------------
# 15. Bayesian calendar-time model with CmdStanR
# ------------------------------------------------------------------------------

stan_code <- '
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
'

stan_file <- file.path(
  script_dir,
  "calendar_time_meta_regression.stan"
)

writeLines(
  stan_code,
  stan_file
)


run_bayes_time <- function(dat, scenario, model) {

  d <- dat
  d$tdec <- (d$year - 2000) / 10

  fit <- model$sample(
    data = list(
      N = nrow(d),
      y = d$yi,
      se = d$sei,
      tdec = d$tdec
    ),

    seed = 20260927,

    chains = 4,
    parallel_chains = min(
      4,
      max(
        1,
        parallel::detectCores(logical = TRUE)
      )
    ),

    iter_warmup = 2000,
    iter_sampling = 4000,

    adapt_delta = 0.995,
    max_treedepth = 15,

    refresh = 500,

    output_dir = file.path(
      outdir,
      paste0("cmdstan_", scenario)
    )
  )

  s <- fit$summary(
    variables = c(
      "beta0",
      "beta1",
      "tau",
      "PF_2000",
      "PF_2010",
      "PF_2020",
      "PF_2025"
    ),

    "median",

    q025 = ~quantile(., 0.025),
    q975 = ~quantile(., 0.975),

    mcse_mean = posterior::mcse_mean,
    rhat = posterior::rhat,
    ess_bulk = posterior::ess_bulk,
    ess_tail = posterior::ess_tail
  )

  b1 <- fit$draws(
    variables = "beta1",
    format = "matrix"
  )[, 1]

  data.frame(
    scenario = scenario,

    beta1_median =
      s$median[s$variable == "beta1"],

    beta1_q025 =
      s$q025[s$variable == "beta1"],

    beta1_q975 =
      s$q975[s$variable == "beta1"],

    Pr_beta1_gt_0 =
      mean(b1 > 0),

    tau_median =
      s$median[s$variable == "tau"],

    PF_2000_median =
      s$median[s$variable == "PF_2000"],

    PF_2010_median =
      s$median[s$variable == "PF_2010"],

    PF_2020_median =
      s$median[s$variable == "PF_2020"],

    PF_2025_median =
      s$median[s$variable == "PF_2025"],

    max_Rhat =
      max(s$rhat, na.rm = TRUE),

    min_ESS_bulk =
      min(s$ess_bulk, na.rm = TRUE),

    min_ESS_tail =
      min(s$ess_tail, na.rm = TRUE),

    stringsAsFactors = FALSE
  )
}


bayes_results <- NULL

if (RUN_STAN) {

  library(cmdstanr)
  library(posterior)

  cat("\nCompiling Bayesian calendar-time Stan model...\n")

  mod <- cmdstanr::cmdstan_model(
    stan_file,
    force_recompile = FALSE
  )

  # Primary Bayesian comparison:
  #   existing core13 versus publication-based core14+Wang.
  #
  # One representative all-16 sensitivity fit is also run.
  bayes_scenarios <- c(
    "core13",
    "core14_wang",
    "all16_jayphi2_ghasrho0.05"
  )

  bayes_results <- do.call(
    rbind,
    lapply(
      bayes_scenarios,
      function(nm) {
        run_bayes_time(
          scenarios[[nm]],
          nm,
          mod
        )
      }
    )
  )

  write.csv(
    bayes_results,
    file.path(outdir, "Bayesian_calendar_time_results.csv"),
    row.names = FALSE
  )
}


# ------------------------------------------------------------------------------
# 16. Main text and Supplement-ready automatically generated summaries
# ------------------------------------------------------------------------------

fmt <- function(x, d = 3) {
  sprintf(paste0("%.", d, "f"), x)
}

primary_time <- subset(
  time_results,
  scenario == "core14_wang"
)

primary_era <- subset(
  era_results,
  scenario == "core14_wang"
)

primary_modern <- subset(
  modern_results,
  scenario == "core14_wang"
)


main_lines <- c(
  "% AUTO-GENERATED: strongest publication-based expanded analysis",
  "% Primary expanded set = core13 + Wang with article ICC=.20 and m=40.",
  "% Jayasinghe and Ghasemi remain sensitivity-only because exact variances",
  "% cannot be recovered from their published aggregate reports.",
  "",
  "\\paragraph{Expanded continuous-outcome production audit.}",
  paste0(
    "The strongest publication-based expansion adds Wang et al. to the ",
    "13-study dataset using an approximate cluster design correction based on ",
    "the trial's stated ICC of 0.20 and planning cluster size of 40. ",
    "The resulting REML--Hartung--Knapp calendar-time coefficient is ",
    "$\\\\widehat\\\\beta_1=",
    fmt(primary_time$beta1),
    "$ per decade (95\\\\% CI $[",
    fmt(primary_time$beta1_ci_lb),
    ",",
    fmt(primary_time$beta1_ci_ub),
    "]$), with residual heterogeneity $\\\\widehat\\\\tau=",
    fmt(primary_time$tau),
    "$."
  ),
  "",
  paste0(
    "The corresponding historical-versus-contemporary coefficient is ",
    fmt(primary_era$beta_modern_minus_hist),
    " (95\\\\% CI $[",
    fmt(primary_era$contrast_ci_lb),
    ",",
    fmt(primary_era$contrast_ci_ub),
    "]$). ",
    "Within the five contemporary comparisons in this expanded primary set, ",
    "the REML pooled log ratio is ",
    fmt(primary_modern$mu_modern),
    ", corresponding to a prevented fraction of ",
    fmt(primary_modern$PF_modern),
    "."
  ),
  "",
  "Jayasinghe and Ghasemi are not assigned primary inverse-variance weights:",
  "Jayasinghe reports the two-year means without an SD or SE, whereas Ghasemi",
  "randomized only four schools and reports neither an outcome ICC nor",
  "school-specific DMFS increments. They are therefore examined through",
  "explicit variance-sensitivity analyses in the Supplement rather than",
  "through a falsely precise single reconstructed variance."
)

if (!is.null(bayes_results)) {

  bp <- subset(
    bayes_results,
    scenario == "core14_wang"
  )

  main_lines <- c(
    main_lines,
    "",
    paste0(
      "Under the prespecified Bayesian calendar-time model, the posterior ",
      "median slope is $",
      fmt(bp$beta1_median),
      "$ per decade (95\\\\% CrI $[",
      fmt(bp$beta1_q025),
      ",",
      fmt(bp$beta1_q975),
      "]$), with ",
      "$\\\\Pr(\\\\beta_1>0\\\\mid\\\\mathcal D)=",
      fmt(bp$Pr_beta1_gt_0),
      "$."
    )
  )
}

writeLines(
  main_lines,
  file.path(
    outdir,
    "main_expanded_time_results.tex"
  )
)


# Supplement table for the scenario analysis.
supp_header <- c(
  "% AUTO-GENERATED sensitivity table",
  "\\begin{table}[ht]",
  "\\centering",
  "\\caption{Sensitivity of the REML--Hartung--Knapp calendar-time coefficient to the treatment of the three contemporary studies with unresolved or cluster-dependent sampling variances.}",
  "\\label{tab:supp-time-variance-sensitivity}",
  "\\small",
  "\\begin{tabular}{lcccc}",
  "\\toprule",
  "Scenario & $k$ & $\\widehat\\beta_1$ & 95\\% CI & $\\widehat\\tau$ \\\\",
  "\\midrule"
)

supp_rows <- character(0)

for (i in seq_len(nrow(time_results))) {

  z <- time_results[i, ]

  label <- gsub("_", "\\\\_", z$scenario, fixed = TRUE)

  supp_rows <- c(
    supp_rows,
    paste0(
      label, " & ",
      z$k, " & ",
      fmt(z$beta1), " & $[",
      fmt(z$beta1_ci_lb), ",",
      fmt(z$beta1_ci_ub), "]$ & ",
      fmt(z$tau),
      " \\\\"
    )
  )
}

supp_footer <- c(
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}"
)

writeLines(
  c(
    supp_header,
    supp_rows,
    supp_footer
  ),
  file.path(
    outdir,
    "supp_time_variance_sensitivity_table.tex"
  )
)


# ------------------------------------------------------------------------------
# 17. Human-readable production recommendation
# ------------------------------------------------------------------------------

recommendation <- c(
  "CONTEMPORARY CONTINUOUS-OUTCOME PRODUCTION AUDIT",
  "================================================",
  "",
  "Primary weighted analysis that can be defended from published information:",
  "  nine historical studies + Munoz-Millan + Latifi-Xhemajli + McMahon + Zeng + Wang.",
  "",
  "Wang:",
  paste0(
    "  naive LRR SE = ", fmt(wang_raw$sei, 4),
    "; publication-based DE=8.8 gives SE = ", fmt(wang_primary$sei, 4),
    "; using actual mean completed cluster size gives SE = ",
    fmt(wang_actual_m$sei, 4), "."
  ),
  "  Recommendation: include, but call the variance approximately cluster-adjusted.",
  "",
  "Jayasinghe:",
  "  The 2-year means 1.50 vs 1.97 are published, but no corresponding SD/SE is reported.",
  "  Recommendation: do not assign a single primary inverse-variance weight.",
  "  Use the quasi-Poisson overdispersion grid as sensitivity only.",
  "",
  "Ghasemi:",
  "  Four schools were randomized (2 per arm); no outcome ICC or school-level DMFS increments are reported.",
  "  Recommendation: do not assign a single primary inverse-variance weight.",
  "  Use the ICC/design-effect grid as sensitivity only.",
  "",
  "The script generates:",
  "  REML_HK_calendar_time_scenarios.csv",
  "  REML_HK_historical_vs_contemporary_scenarios.csv",
  "  REML_HK_contemporary_pool_scenarios.csv",
  "  Jayasinghe_Ghasemi_variance_sensitivity_grid.csv",
  "  ML_AICc_model_comparison.csv",
  "  leave_one_contemporary_out.csv",
  "  contemporary_variance_audit.csv",
  "  temporal_core13.pdf",
  "  temporal_core14_wang.pdf",
  "  temporal_all16_sensitivity.pdf",
  "  Bayesian_calendar_time_results.csv (if RUN_STAN=TRUE)",
  "  main_expanded_time_results.tex",
  "  supp_time_variance_sensitivity_table.tex"
)

writeLines(
  recommendation,
  file.path(
    outdir,
    "PRODUCTION_RECOMMENDATION.txt"
  )
)


# ------------------------------------------------------------------------------
# 18. Independent benchmark checks
# ------------------------------------------------------------------------------

# These checks are deliberately loose.  They detect accidental changes in the
# data/model but do not force exact equality across metafor versions.

b13 <- subset(
  time_results,
  scenario == "core13"
)$beta1

b14 <- subset(
  time_results,
  scenario == "core14_wang"
)$beta1

if (abs(b13 - 0.1116) > 0.01) {
  warning("core13 time slope differs materially from the independent benchmark.")
}

if (abs(b14 - 0.0986) > 0.01) {
  warning("core14+Wang time slope differs materially from the independent benchmark.")
}


# ------------------------------------------------------------------------------
# 19. Console summary
# ------------------------------------------------------------------------------

cat("\n============================================================\n")
cat("CONTINUOUS-OUTCOME PRODUCTION / SENSITIVITY ANALYSIS COMPLETE\n")
cat("============================================================\n\n")

cat("Primary publication-based expanded scenario: core14_wang\n\n")

print(
  subset(
    time_results,
    scenario %in% c(
      "core13",
      "core14_wang",
      "core14_wang_actual_cluster_size",
      "core15_wang_jay_phi2",
      "all16_jayphi2_ghasrho0.05"
    )
  ),
  digits = 5,
  row.names = FALSE
)

cat("\nOutput directory:\n")
cat(normalizePath(outdir), "\n")

# ==============================================================================
# END
# ==============================================================================
