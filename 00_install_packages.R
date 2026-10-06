# ==============================================================================
# 00_install_packages.R
# Packages used by the fluoride-varnish production analyses.
# ==============================================================================

cran_packages <- c("metafor", "bayesmeta", "posterior")

missing <- cran_packages[
  !vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing) > 0) {
  install.packages(missing)
}

# CmdStanR is commonly installed from the Stan repository rather than CRAN.
if (!requireNamespace("cmdstanr", quietly = TRUE)) {
  install.packages(
    "cmdstanr",
    repos = c(
      "https://stan-dev.r-universe.dev",
      getOption("repos")
    )
  )
}

# Install CmdStan once if needed.
if (requireNamespace("cmdstanr", quietly = TRUE)) {
  ok <- tryCatch({
    cmdstanr::cmdstan_version()
    TRUE
  }, error = function(e) FALSE)

  if (!ok) {
    message("CmdStan is not installed. Installing the current CmdStan release...")
    cmdstanr::install_cmdstan()
  }
}

message("Package setup complete.")
