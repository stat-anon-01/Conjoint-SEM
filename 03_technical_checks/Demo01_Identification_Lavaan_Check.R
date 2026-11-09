# DEMONSTRATION 1: Numerical identification and scale diagnostics
# Run from the repository's 03_technical_checks directory:
#   source("Demo01_Identification_Lavaan_Check.R")
#
# This script reuses DEMO 1's original DGP and fitted model syntax.
# It writes Demo01_Identification_Lavaan_Results.csv in this directory.
# Numerical full rank at an estimated solution supports local, not global,
# identification. It does not constitute an analytical identification proof.

if (!requireNamespace("lavaan", quietly = TRUE))
  stop('Install lavaan with install.packages("lavaan")')
library(lavaan)

# Locate the repository independently of the R working directory.
# When sourced, R records the executing file in a source frame ($ofile).
source_files <- lapply(sys.frames(), function(frame) frame$ofile)
source_files <- Filter(function(x) is.character(x) && length(x) == 1L &&
                       grepl("Demo01_Identification_Lavaan_Check\\.R$", x),
                       source_files)
script_dir <- if (length(source_files)) {
  dirname(normalizePath(tail(source_files, 1)[[1]], mustWork = TRUE))
} else {
  getwd()
}

load_demo01 <- function() {
  candidates <- unique(normalizePath(c(
    script_dir,
    getwd(),
    file.path(getwd(), "03_technical_checks"),
    file.path(getwd(), "03_check")
  ), mustWork = FALSE))
  demo_dirs <- unique(unlist(lapply(candidates, function(d) c(
    file.path(d, "../02_demonstrations/demo01"),
    file.path(d, "02_demonstrations/demo01"),
    file.path(d, "demo01")
  ))))
  required <- c("Demo01_00_DGP.R", "Demo01_01_ANALYSIS.R")
  found <- Filter(function(d) all(file.exists(file.path(d, required))), demo_dirs)
  if (!length(found)) {
    stop("Cannot locate DEMO 1 scripts. Keep 02_demonstrations/demo01 and ",
         "03_technical_checks (or 03_check) in the same GITHUB folder. ",
         "Current working directory: ", getwd())
  }
  e <- new.env(parent = globalenv())
  # chdir=TRUE resolves relative paths inside the sourced analysis scripts.
  source(file.path(found[[1]], required[[1]]), local = e, chdir = TRUE)
  source(file.path(found[[1]], required[[2]]), local = e, chdir = TRUE)
  e
}

demo <- load_demo01()
fits <- list(M00 = demo$fit_M00, M10 = demo$fit_M10,
             M01 = demo$fit_M01, M11 = demo$fit_M11)
# DGP quantities used for scale normalization below:
Phi <- demo$Phi
mu <- demo$mu
tau_B <- demo$tau_B
sd_B <- demo$sd_B

# 1. Basic fit and identification diagnostics.
# lavInspect(fit, "delta") is the derivative of model-implied sample
# statistics with respect to model parameters, evaluated at the fit.
# Numerical rank depends on the tolerance and fitted parameter point.
check_one <- function(fit) {
  delta <- lavInspect(fit, "delta")
  info <- lavInspect(fit, "information")
  nfree <- lavInspect(fit, "npar")
  pe <- parameterEstimates(fit)
  pt <- parTable(fit)
  residual_variances <- subset(pe, op == "~~" & lhs == rhs &
                                     lhs %in% paste0("B", 1:8))
  data.frame(
    converged = lavInspect(fit, "converged"),
    free_parameters = nfree,
    df = fitMeasures(fit, "df"),
    delta_rank = qr(delta, tol = 1e-7)$rank,
    information_rank = qr(info, tol = 1e-7)$rank,
    smallest_delta_singular_value = min(svd(delta)$d),
    negative_residual_variances = sum(residual_variances$est < 0),
    missing_standard_errors = sum(is.na(pt$se[pt$free > 0])),
    row.names = NULL
  )
}

audit <- do.call(rbind, lapply(fits, check_one))
audit <- cbind(model = names(fits), audit)
print(audit, row.names = FALSE)

# Rank checks are empirical diagnostics at the estimated parameter point.
# A full-rank derivative does not establish global identification.
write.csv(audit, file.path(script_dir, "Demo01_Identification_Lavaan_Results.csv"),
          row.names = FALSE)

# 2. DGP-to-model scale normalization.
# Var(F1) is fixed at 1 in the fitted binary models.
# The common DGP threshold (.50) is absorbed into the F0 mean.
scale_multiplier <- sqrt(1 / Phi[2, 2])
reference_means <- scale_multiplier * mu
reference_means["F0"] <- scale_multiplier * (mu["F0"] - tau_B["B1"])
reference_covariance <- scale_multiplier^2 * Phi
reference_residual_variances <- scale_multiplier^2 * sd_B^2
reference_B8_threshold <- scale_multiplier * (tau_B["B8"] - tau_B["B1"])

cat("\nScale multiplier:", scale_multiplier, "\n")
cat("Reference factor means:\n")
print(reference_means)
cat("Reference factor covariance matrix:\n")
print(reference_covariance)
cat("Reference residual variances:\n")
print(reference_residual_variances)
cat("Reference B8 threshold deviation:", reference_B8_threshold, "\n")

# 3. Inspect fitted estimates and WLSMV nested-model tests.
print(parameterEstimates(fits$M11))
print(lavTestLRT(fits$M00, fits$M10))
print(lavTestLRT(fits$M01, fits$M11))
print(lavTestLRT(fits$M00, fits$M01))
print(lavTestLRT(fits$M10, fits$M11))
