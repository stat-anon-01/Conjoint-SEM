# DEMONSTRATION 3: Cross-format preference comparability
# 01_ANALYSIS.R — Fit the models used in the article
#
# Run in order from this folder:
#   1. source("01_ANALYSIS.R")   # generates data and fits models
#   2. source("02_TABLE.R")      # exports tables without refitting
#
# Dependencies: MASS (data generation), lavaan (estimation).
# Models and their fixed loadings/restrictions are written explicitly so
# readers can inspect the assumptions; do not substitute generated syntax.
# M0: baseline; M1: equal means; M2: equal covariance matrices; M3: common individual factors.
# Running from a fresh R session is recommended.
# ------------------------------------------------------------

if (!requireNamespace("lavaan", quietly=TRUE)) stop("Install lavaan")
source("00_DGP.R")
sim <- make_study2_data(N=1000L, formats=c("R","B"), rho=.70, seed=12345)
dat <- sim$data
ordered_items <- paste0("B", 1:8)

model_M0 <- '
  # M0: explicitly specified fixed conjoint measurement model
  R_F0 =~ 1*R1 + 1*R2 + 1*R3 + 1*R4 + 1*R5 + 1*R6 + 1*R7 + 1*R8
  R_F1 =~ -1*R1 + 1*R2 + -1*R3 + 1*R4 + -1*R5 + 1*R6 + -1*R7 + 1*R8
  R_F2 =~ -1*R1 + -1*R2 + 1*R3 + 1*R4 + -1*R5 + -1*R6 + 1*R7 + 1*R8
  R_F3 =~ -1*R1 + -1*R2 + -1*R3 + -1*R4 + 1*R5 + 1*R6 + 1*R7 + 1*R8
  B_F0 =~ 1*B1 + 1*B2 + 1*B3 + 1*B4 + 1*B5 + 1*B6 + 1*B7 + 1*B8
  B_F1 =~ -1*B1 + 1*B2 + -1*B3 + 1*B4 + -1*B5 + 1*B6 + -1*B7 + 1*B8
  B_F2 =~ -1*B1 + -1*B2 + 1*B3 + 1*B4 + -1*B5 + -1*B6 + 1*B7 + 1*B8
  B_F3 =~ -1*B1 + -1*B2 + -1*B3 + -1*B4 + 1*B5 + 1*B6 + 1*B7 + 1*B8
  # Rating observed intercepts fixed to zero
  R1 ~ 0*1
  R2 ~ 0*1
  R3 ~ 0*1
  R4 ~ 0*1
  R5 ~ 0*1
  R6 ~ 0*1
  R7 ~ 0*1
  R8 ~ 0*1
  # Binary thresholds fixed at zero (location normalization)
  B1 | 0*t1
  B2 | 0*t1
  B3 | 0*t1
  B4 | 0*t1
  B5 | 0*t1
  B6 | 0*t1
  B7 | 0*t1
  B8 | 0*t1
  # Heterogeneous item residuals
  R1 ~~ R1
  R2 ~~ R2
  R3 ~~ R3
  R4 ~~ R4
  R5 ~~ R5
  R6 ~~ R6
  R7 ~~ R7
  R8 ~~ R8
  B1 ~~ B1
  B2 ~~ B2
  B3 ~~ B3
  B4 ~~ B4
  B5 ~~ B5
  B6 ~~ B6
  B7 ~~ B7
  B8 ~~ B8
  # Format-specific latent means
  R_F0 ~ mR0*1
  B_F0 ~ mB0*1
  R_F1 ~ mR1*1
  B_F1 ~ mB1*1
  R_F2 ~ mR2*1
  B_F2 ~ mB2*1
  R_F3 ~ mR3*1
  B_F3 ~ mB3*1
  # Factor variances/covariances; F1 variance anchored across formats
  R_F0 ~~ vR00*R_F0
  B_F0 ~~ vB00*B_F0
  R_F0 ~~ vR01*R_F1
  B_F0 ~~ vB01*B_F1
  R_F0 ~~ vR02*R_F2
  B_F0 ~~ vB02*B_F2
  R_F0 ~~ vR03*R_F3
  B_F0 ~~ vB03*B_F3
  R_F1 ~~ v11*R_F1
  B_F1 ~~ v11*B_F1
  R_F1 ~~ vR12*R_F2
  B_F1 ~~ vB12*B_F2
  R_F1 ~~ vR13*R_F3
  B_F1 ~~ vB13*B_F3
  R_F2 ~~ vR22*R_F2
  B_F2 ~~ vB22*B_F2
  R_F2 ~~ vR23*R_F3
  B_F2 ~~ vB23*B_F3
  R_F3 ~~ vR33*R_F3
  B_F3 ~~ vB33*B_F3

  # All cross-format factor covariances free
  R_F0 ~~ c00*B_F0
  R_F0 ~~ c01*B_F1
  R_F0 ~~ c02*B_F2
  R_F0 ~~ c03*B_F3
  R_F1 ~~ c10*B_F0
  R_F1 ~~ c11*B_F1
  R_F1 ~~ c12*B_F2
  R_F1 ~~ c13*B_F3
  R_F2 ~~ c20*B_F0
  R_F2 ~~ c21*B_F1
  R_F2 ~~ c22*B_F2
  R_F2 ~~ c23*B_F3
  R_F3 ~~ c30*B_F0
  R_F3 ~~ c31*B_F1
  R_F3 ~~ c32*B_F2
  R_F3 ~~ c33*B_F3
'

model_M1 <- '
  # M1: explicitly specified fixed conjoint measurement model
  R_F0 =~ 1*R1 + 1*R2 + 1*R3 + 1*R4 + 1*R5 + 1*R6 + 1*R7 + 1*R8
  R_F1 =~ -1*R1 + 1*R2 + -1*R3 + 1*R4 + -1*R5 + 1*R6 + -1*R7 + 1*R8
  R_F2 =~ -1*R1 + -1*R2 + 1*R3 + 1*R4 + -1*R5 + -1*R6 + 1*R7 + 1*R8
  R_F3 =~ -1*R1 + -1*R2 + -1*R3 + -1*R4 + 1*R5 + 1*R6 + 1*R7 + 1*R8
  B_F0 =~ 1*B1 + 1*B2 + 1*B3 + 1*B4 + 1*B5 + 1*B6 + 1*B7 + 1*B8
  B_F1 =~ -1*B1 + 1*B2 + -1*B3 + 1*B4 + -1*B5 + 1*B6 + -1*B7 + 1*B8
  B_F2 =~ -1*B1 + -1*B2 + 1*B3 + 1*B4 + -1*B5 + -1*B6 + 1*B7 + 1*B8
  B_F3 =~ -1*B1 + -1*B2 + -1*B3 + -1*B4 + 1*B5 + 1*B6 + 1*B7 + 1*B8
  # Rating observed intercepts fixed to zero
  R1 ~ 0*1
  R2 ~ 0*1
  R3 ~ 0*1
  R4 ~ 0*1
  R5 ~ 0*1
  R6 ~ 0*1
  R7 ~ 0*1
  R8 ~ 0*1
  # Binary thresholds fixed at zero (location normalization)
  B1 | 0*t1
  B2 | 0*t1
  B3 | 0*t1
  B4 | 0*t1
  B5 | 0*t1
  B6 | 0*t1
  B7 | 0*t1
  B8 | 0*t1
  # Heterogeneous item residuals
  R1 ~~ R1
  R2 ~~ R2
  R3 ~~ R3
  R4 ~~ R4
  R5 ~~ R5
  R6 ~~ R6
  R7 ~~ R7
  R8 ~~ R8
  B1 ~~ B1
  B2 ~~ B2
  B3 ~~ B3
  B4 ~~ B4
  B5 ~~ B5
  B6 ~~ B6
  B7 ~~ B7
  B8 ~~ B8
  # Format-specific latent means
  R_F0 ~ mR0*1
  B_F0 ~ mB0*1
  R_F1 ~ m1*1
  B_F1 ~ m1*1
  R_F2 ~ m2*1
  B_F2 ~ m2*1
  R_F3 ~ m3*1
  B_F3 ~ m3*1
  # Factor variances/covariances; F1 variance anchored across formats
  R_F0 ~~ vR00*R_F0
  B_F0 ~~ vB00*B_F0
  R_F0 ~~ vR01*R_F1
  B_F0 ~~ vB01*B_F1
  R_F0 ~~ vR02*R_F2
  B_F0 ~~ vB02*B_F2
  R_F0 ~~ vR03*R_F3
  B_F0 ~~ vB03*B_F3
  R_F1 ~~ v11*R_F1
  B_F1 ~~ v11*B_F1
  R_F1 ~~ vR12*R_F2
  B_F1 ~~ vB12*B_F2
  R_F1 ~~ vR13*R_F3
  B_F1 ~~ vB13*B_F3
  R_F2 ~~ vR22*R_F2
  B_F2 ~~ vB22*B_F2
  R_F2 ~~ vR23*R_F3
  B_F2 ~~ vB23*B_F3
  R_F3 ~~ vR33*R_F3
  B_F3 ~~ vB33*B_F3
  # All cross-format factor covariances free
  R_F0 ~~ c00*B_F0
  R_F0 ~~ c01*B_F1
  R_F0 ~~ c02*B_F2
  R_F0 ~~ c03*B_F3
  R_F1 ~~ c10*B_F0
  R_F1 ~~ c11*B_F1
  R_F1 ~~ c12*B_F2
  R_F1 ~~ c13*B_F3
  R_F2 ~~ c20*B_F0
  R_F2 ~~ c21*B_F1
  R_F2 ~~ c22*B_F2
  R_F2 ~~ c23*B_F3
  R_F3 ~~ c30*B_F0
  R_F3 ~~ c31*B_F1
  R_F3 ~~ c32*B_F2
  R_F3 ~~ c33*B_F3
'

model_M2 <- '
  # M2: explicitly specified fixed conjoint measurement model
  R_F0 =~ 1*R1 + 1*R2 + 1*R3 + 1*R4 + 1*R5 + 1*R6 + 1*R7 + 1*R8
  R_F1 =~ -1*R1 + 1*R2 + -1*R3 + 1*R4 + -1*R5 + 1*R6 + -1*R7 + 1*R8
  R_F2 =~ -1*R1 + -1*R2 + 1*R3 + 1*R4 + -1*R5 + -1*R6 + 1*R7 + 1*R8
  R_F3 =~ -1*R1 + -1*R2 + -1*R3 + -1*R4 + 1*R5 + 1*R6 + 1*R7 + 1*R8
  B_F0 =~ 1*B1 + 1*B2 + 1*B3 + 1*B4 + 1*B5 + 1*B6 + 1*B7 + 1*B8
  B_F1 =~ -1*B1 + 1*B2 + -1*B3 + 1*B4 + -1*B5 + 1*B6 + -1*B7 + 1*B8
  B_F2 =~ -1*B1 + -1*B2 + 1*B3 + 1*B4 + -1*B5 + -1*B6 + 1*B7 + 1*B8
  B_F3 =~ -1*B1 + -1*B2 + -1*B3 + -1*B4 + 1*B5 + 1*B6 + 1*B7 + 1*B8
  # Rating observed intercepts fixed to zero
  R1 ~ 0*1
  R2 ~ 0*1
  R3 ~ 0*1
  R4 ~ 0*1
  R5 ~ 0*1
  R6 ~ 0*1
  R7 ~ 0*1
  R8 ~ 0*1
  # Binary thresholds fixed at zero (location normalization)
  B1 | 0*t1
  B2 | 0*t1
  B3 | 0*t1
  B4 | 0*t1
  B5 | 0*t1
  B6 | 0*t1
  B7 | 0*t1
  B8 | 0*t1
  # Heterogeneous item residuals
  R1 ~~ R1
  R2 ~~ R2
  R3 ~~ R3
  R4 ~~ R4
  R5 ~~ R5
  R6 ~~ R6
  R7 ~~ R7
  R8 ~~ R8
  B1 ~~ B1
  B2 ~~ B2
  B3 ~~ B3
  B4 ~~ B4
  B5 ~~ B5
  B6 ~~ B6
  B7 ~~ B7
  B8 ~~ B8
  # Format-specific latent means
  R_F0 ~ mR0*1
  B_F0 ~ mB0*1
  R_F1 ~ m1*1
  B_F1 ~ m1*1
  R_F2 ~ m2*1
  B_F2 ~ m2*1
  R_F3 ~ m3*1
  B_F3 ~ m3*1
  # Factor variances/covariances; F1 variance anchored across formats
  R_F0 ~~ v00*R_F0
  B_F0 ~~ v00*B_F0
  R_F0 ~~ v01*R_F1
  B_F0 ~~ v01*B_F1
  R_F0 ~~ v02*R_F2
  B_F0 ~~ v02*B_F2
  R_F0 ~~ v03*R_F3
  B_F0 ~~ v03*B_F3
  R_F1 ~~ v11*R_F1
  B_F1 ~~ v11*B_F1
  R_F1 ~~ v12*R_F2
  B_F1 ~~ v12*B_F2
  R_F1 ~~ v13*R_F3
  B_F1 ~~ v13*B_F3
  R_F2 ~~ v22*R_F2
  B_F2 ~~ v22*B_F2
  R_F2 ~~ v23*R_F3
  B_F2 ~~ v23*B_F3
  R_F3 ~~ v33*R_F3
  B_F3 ~~ v33*B_F3
  # All cross-format factor covariances free
  R_F0 ~~ c00*B_F0
  R_F0 ~~ c01*B_F1
  R_F0 ~~ c02*B_F2
  R_F0 ~~ c03*B_F3
  R_F1 ~~ c10*B_F0
  R_F1 ~~ c11*B_F1
  R_F1 ~~ c12*B_F2
  R_F1 ~~ c13*B_F3
  R_F2 ~~ c20*B_F0
  R_F2 ~~ c21*B_F1
  R_F2 ~~ c22*B_F2
  R_F2 ~~ c23*B_F3
  R_F3 ~~ c30*B_F0
  R_F3 ~~ c31*B_F1
  R_F3 ~~ c32*B_F2
  R_F3 ~~ c33*B_F3
'

model_M3 <- '
  # M3: explicitly specified fixed conjoint measurement model
  R_F0 =~ 1*R1 + 1*R2 + 1*R3 + 1*R4 + 1*R5 + 1*R6 + 1*R7 + 1*R8
  B_F0 =~ 1*B1 + 1*B2 + 1*B3 + 1*B4 + 1*B5 + 1*B6 + 1*B7 + 1*B8
  F1 =~ -1*R1 + 1*R2 + -1*R3 + 1*R4 + -1*R5 + 1*R6 + -1*R7 + 1*R8 + -1*B1 + 1*B2 + -1*B3 + 1*B4 + -1*B5 + 1*B6 + -1*B7 + 1*B8
  F2 =~ -1*R1 + -1*R2 + 1*R3 + 1*R4 + -1*R5 + -1*R6 + 1*R7 + 1*R8 + -1*B1 + -1*B2 + 1*B3 + 1*B4 + -1*B5 + -1*B6 + 1*B7 + 1*B8
  F3 =~ -1*R1 + -1*R2 + -1*R3 + -1*R4 + 1*R5 + 1*R6 + 1*R7 + 1*R8 + -1*B1 + -1*B2 + -1*B3 + -1*B4 + 1*B5 + 1*B6 + 1*B7 + 1*B8
  # Rating observed intercepts fixed to zero
  R1 ~ 0*1
  R2 ~ 0*1
  R3 ~ 0*1
  R4 ~ 0*1
  R5 ~ 0*1
  R6 ~ 0*1
  R7 ~ 0*1
  R8 ~ 0*1
  # Binary thresholds fixed at zero (location normalization)
  B1 | 0*t1
  B2 | 0*t1
  B3 | 0*t1
  B4 | 0*t1
  B5 | 0*t1
  B6 | 0*t1
  B7 | 0*t1
  B8 | 0*t1
  # Heterogeneous item residuals
  R1 ~~ R1
  R2 ~~ R2
  R3 ~~ R3
  R4 ~~ R4
  R5 ~~ R5
  R6 ~~ R6
  R7 ~~ R7
  R8 ~~ R8
  B1 ~~ B1
  B2 ~~ B2
  B3 ~~ B3
  B4 ~~ B4
  B5 ~~ B5
  B6 ~~ B6
  B7 ~~ B7
  B8 ~~ B8
  # Shared slope factors F1-F3; separate F0 location factors
  R_F0 ~ 1
  B_F0 ~ 1
  F1 ~ 1
  F2 ~ 1
  F3 ~ 1
  R_F0 ~~ R_F0
  R_F0 ~~ B_F0
  R_F0 ~~ F1
  R_F0 ~~ F2
  R_F0 ~~ F3
  B_F0 ~~ B_F0
  B_F0 ~~ F1
  B_F0 ~~ F2
  B_F0 ~~ F3
  F1 ~~ F1
  F1 ~~ F2
  F1 ~~ F3
  F2 ~~ F2
  F2 ~~ F3
  F3 ~~ F3
'

# Fit each model directly. lavaan warnings remain visible.

# Model M0
fit_M0 <- lavaan::sem(
  model = model_M0, data = dat, ordered = ordered_items,
  estimator = "WLSMV", parameterization = "theta",
  meanstructure = TRUE, std.lv = FALSE,
  auto.fix.first = FALSE, auto.var = FALSE, auto.th = FALSE,
  auto.cov.lv.x = FALSE, auto.cov.y = FALSE
)
summary(fit_M0, fit.measures = TRUE)

# Model M1
fit_M1 <- lavaan::sem(
  model = model_M1, data = dat, ordered = ordered_items,
  estimator = "WLSMV", parameterization = "theta",
  meanstructure = TRUE, std.lv = FALSE,
  auto.fix.first = FALSE, auto.var = FALSE, auto.th = FALSE,
  auto.cov.lv.x = FALSE, auto.cov.y = FALSE
)
summary(fit_M1, fit.measures = TRUE)

# Model M2
fit_M2 <- lavaan::sem(
  model = model_M2, data = dat, ordered = ordered_items,
  estimator = "WLSMV", parameterization = "theta",
  meanstructure = TRUE, std.lv = FALSE,
  auto.fix.first = FALSE, auto.var = FALSE, auto.th = FALSE,
  auto.cov.lv.x = FALSE, auto.cov.y = FALSE
)
summary(fit_M2, fit.measures = TRUE)

# Model M3
fit_M3 <- lavaan::sem(
  model = model_M3, data = dat, ordered = ordered_items,
  estimator = "WLSMV", parameterization = "theta",
  meanstructure = TRUE, std.lv = FALSE,
  auto.fix.first = FALSE, auto.var = FALSE, auto.th = FALSE,
  auto.cov.lv.x = FALSE, auto.cov.y = FALSE
)
summary(fit_M3, fit.measures = TRUE)

models <- list(M0 = model_M0, M1 = model_M1, M2 = model_M2, M3 = model_M3)
fits <- list(M0 = fit_M0, M1 = fit_M1, M2 = fit_M2, M3 = fit_M3)

# Constrained information rank. lavaan may retain both raw free parameters
# for an equality label in lavInspect(fit, "information"). Therefore raw
# rank / raw npar is NOT a valid identification check in such models.
#
# Here all equality constraints in M0-M2 are simple shared parameter labels.
# The transformation T maps one parameter per label group to the original
# free-parameter coordinates. I_reduced = t(T) %*% I %*% T.
# This function intentionally stops on unsupported explicit/nonlinear
# constraints rather than silently reporting a misleading rank.
constrained_information_rank <- function(fit, tol = 1e-8) {
  pt <- lavaan::parTable(fit)
  info <- as.matrix(lavaan::lavInspect(fit, "information"))
  p <- ncol(info)
  stopifnot(nrow(info) == p)

  if (any(pt$op %in% c("<", ">", ":="))) {
    stop("Unsupported nonlinear/inequality constraints in rank diagnostic")
  }
  # lavaan may expand repeated labels into explicit == rows. These are
  # supported only if they correspond to the same-label groups below.
  equality_rows <- pt[pt$op == "==", , drop = FALSE]

  free_rows <- pt[pt$free > 0L, , drop = FALSE]
  ids <- sort(unique(free_rows$free))
  if (!identical(ids, seq_len(p))) {
    stop("Free parameter IDs do not match information matrix dimensions")
  }

  # A blank label is unique to its free parameter; repeated nonblank
  # labels define equality groups.
  lab <- as.character(free_rows$label)
  lab[is.na(lab)] <- ""
  key <- ifelse(nzchar(lab), paste0("label:", lab),
                paste0("free:", free_rows$free))
  # Collapse rows referring to the same raw free ID, if any.
  key_by_id <- vapply(ids, function(i) {
    k <- unique(key[free_rows$free == i])
    if (length(k) != 1L) stop("Inconsistent labels for free parameter ", i)
    k
  }, character(1))
  groups <- match(key_by_id, unique(key_by_id))
  q <- max(groups)
  if (nrow(equality_rows) > 0L && nrow(equality_rows) != p - q)
    stop("Explicit equalities not fully explained by repeated labels")
  T <- matrix(0, nrow = p, ncol = q)
  T[cbind(seq_len(p), groups)] <- 1
  reduced <- crossprod(T, info %*% T)
  reduced <- (reduced + t(reduced)) / 2
  eigenvalues <- eigen(reduced, symmetric = TRUE, only.values = TRUE)$values
  relative <- eigenvalues / max(abs(eigenvalues))
  rank <- sum(relative > tol)
  list(raw_npar = p, raw_rank = qr(info)$rank,
       n_equalities = p - q, effective_npar = q,
       rank = rank, full_rank = rank == q,
       min_relative_eigenvalue = min(relative))
}

rank_results <- list()
for (nm in names(fits)) {
  f <- fits[[nm]]
  cat("\n-----", nm, "-----\n")
  if (is.null(f)) next
  converged <- lavaan::lavInspect(f, "converged")
  cat("Converged:", converged, "\n")
  pt <- lavaan::parTable(f)
  rr <- subset(pt, op == "~~" & lhs == rhs & lhs %in% ordered_items)
  n_binary_resid <- length(unique(rr$free[rr$free > 0L]))
  cat("Binary residual variances with distinct free IDs:",
      n_binary_resid, "/ 8\n")
  if (n_binary_resid != 8L)
    warning(nm, ": Binary residuals are not all independently free")

  rank_results[[nm]] <- tryCatch(
    constrained_information_rank(f),
    error = function(e) { warning(nm, ": ", conditionMessage(e)); NULL }
  )
  d <- rank_results[[nm]]
  if (!is.null(d)) {
    cat("Raw information rank:", d$raw_rank, "/", d$raw_npar,
        "(ignores equality constraints)\n")
    cat("Equality constraints:", d$n_equalities, "\n")
    cat("Constrained information rank:", d$rank, "/",
        d$effective_npar, "\n")
    cat("Minimum relative eigenvalue:",
        format(d$min_relative_eigenvalue, digits = 5), "\n")
    if (!d$full_rank) warning(nm, ": constrained information is rank deficient")
  }
  print(tryCatch(lavaan::fitMeasures(f, c("chisq.scaled", "df.scaled",
    "cfi.scaled", "rmsea.scaled", "srmr")),
    error = function(e) conditionMessage(e)))
}

# Confirm nested restrictions and their intended effective parameter counts.
# M1 adds 3 slope-mean equalities; M2 adds 9 further covariance
# equalities (including F0 terms), beyond the F1 variance scale anchor.
if (all(vapply(rank_results[c("M0", "M1", "M2")],
               function(x) !is.null(x), logical(1)))) {
  neff <- vapply(rank_results[c("M0", "M1", "M2")],
                 function(x) x$effective_npar, integer(1))
  cat("\nEffective parameter counts (M0/M1/M2):",
      paste(neff, collapse = " / "), "\n")
  if (!identical(as.integer(neff[1] - neff[2]), 3L) ||
      !identical(as.integer(neff[2] - neff[3]), 9L))
    warning("Unexpected M0-M1-M2 equality-constraint counts")
}

for (pair in list(c("M0","M1"),c("M1","M2"))) {
  a <- fits[[pair[1]]]; b <- fits[[pair[2]]]
  if (is.null(a)||is.null(b)) next
  if (!all(vapply(list(a,b),lavaan::lavInspect,logical(1),"converged"))) next
  cat("\nWLSMV comparison:",paste(pair,collapse=" vs "),"\n")
  print(tryCatch(lavaan::lavTestLRT(a, b, method = "satorra.2000"),
                 error = function(e) conditionMessage(e)))
}
cat("\nM2 vs M3: no ordinary chi-square difference test (boundary).\n")
saveRDS(list(sim=sim, fits=fits, models=models, rank_diagnostics=rank_results,
  note="F0 excluded from cross-format mean equivalence; M3 shares slopes only"),
  "Demo03_results.rds")
