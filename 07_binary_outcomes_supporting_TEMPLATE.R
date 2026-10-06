# ==============================================================================
# 07_binary_outcomes_supporting_TEMPLATE.R
#
# Supporting contemporary binary-outcome synthesis.
#
# IMPORTANT:
# The uploaded R_Code.zip did not contain the frozen child-level event-count
# extraction (or design-adjusted log-risk-ratio standard errors) used for the
# manuscript's supportive binary analysis.  This script is therefore a CLEAN
# TEMPLATE, not a claimed reproduction of those results.
#
# Fill data/binary_event_data_TEMPLATE.csv with verified counts (or directly
# with a design-adjusted logRR and SE), then run this script.
#
# Manuscript-scale reported crude RRs were approximately:
#   Munoz-Millan 0.81, McMahon 0.85, Wang 0.72, Zeng 1.08, He 0.81.
#
# Do NOT reconstruct cluster-randomized variances from nominal child counts if
# a design-adjusted effect/SE is available.
# ==============================================================================

if (!requireNamespace("metafor", quietly = TRUE)) {
  stop("Install 'metafor' first: install.packages('metafor')")
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
input_file <- file.path(script_dir, "data", "binary_event_data_TEMPLATE.csv")
outdir <- file.path(script_dir, "output_07_binary_supporting")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

d <- read.csv(input_file, stringsAsFactors = FALSE)

required <- c(
  "study", "year", "event_T", "n_T", "event_C", "n_C",
  "logRR_adjusted", "SE_logRR_adjusted"
)
stopifnot(all(required %in% names(d)))

# Prefer a verified design-adjusted logRR/SE when supplied.
use_adjusted <- is.finite(d$logRR_adjusted) & is.finite(d$SE_logRR_adjusted)

# Otherwise calculate a crude individual-level log risk ratio from event counts.
# This is appropriate only when such a crude variance is scientifically valid.
can_crude <- with(
  d,
  is.finite(event_T) & is.finite(n_T) &
  is.finite(event_C) & is.finite(n_C) &
  event_T > 0 & event_C > 0 &
  event_T < n_T & event_C < n_C
)

d$yi <- NA_real_
d$sei <- NA_real_

d$yi[use_adjusted] <- d$logRR_adjusted[use_adjusted]
d$sei[use_adjusted] <- d$SE_logRR_adjusted[use_adjusted]

idx <- !use_adjusted & can_crude
d$yi[idx] <- log(
  (d$event_T[idx] / d$n_T[idx]) /
  (d$event_C[idx] / d$n_C[idx])
)
d$sei[idx] <- sqrt(
  1 / d$event_T[idx] - 1 / d$n_T[idx] +
  1 / d$event_C[idx] - 1 / d$n_C[idx]
)

if (any(!is.finite(d$yi) | !is.finite(d$sei))) {
  stop(
    "Binary input is incomplete. Fill verified event counts or adjusted ",
    "logRR/SE values in data/binary_event_data_TEMPLATE.csv."
  )
}

d$vi <- d$sei^2

fit <- metafor::rma.uni(
  yi = yi,
  vi = vi,
  data = d,
  method = "REML",
  test = "knha"
)

pred <- predict(fit, transf = exp)

write.csv(
  d,
  file.path(outdir, "binary_effects_used.csv"),
  row.names = FALSE
)
capture.output(
  fit,
  file = file.path(outdir, "binary_REML_HK_fit.txt")
)
capture.output(
  sessionInfo(),
  file = file.path(outdir, "R_sessionInfo.txt")
)

cat("\nSupporting binary-outcome REML-HK synthesis\n")
print(fit, digits = 4)
cat("\nPooled RR =", exp(as.numeric(coef(fit)[1])), "\n")
cat("Output written to: ", normalizePath(outdir), "\n", sep = "")
