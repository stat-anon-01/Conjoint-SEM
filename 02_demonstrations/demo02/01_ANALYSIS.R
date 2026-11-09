# DEMONSTRATION 2: Residual dependence in paired comparisons
# 01_ANALYSIS.R — Fit the models used in the article
#
# Run in order from this folder:
#   1. source("01_ANALYSIS.R")   # generates data and fits models
#   2. source("02_TABLE.R")      # exports tables without refitting
#
# Dependencies: MASS (data generation), lavaan (estimation).
# Models and their fixed loadings/restrictions are written explicitly so
# readers can inspect the assumptions; do not substitute generated syntax.
# A: profile-level; B: pair-specific; C: both.
# Running from a fresh R session is recommended.
# ------------------------------------------------------------

if (!requireNamespace('lavaan',quietly=TRUE)) stop('Install lavaan')
source("00_DGP.R")

# B: original S01 SEM. The factor loadings and all thresholds are explicit.
# Factor F1 variance=1 fixes latent scale (NOT the DGP's .25).


# A: product errors only, pair-specific residual variances constrained to zero.
model_A <- '
  F1 =~ -2*P1_2 + -2*P1_4 + -2*P1_6 + -2*P1_8 + 2*P2_3 + 2*P2_5 + 2*P2_7 + -2*P3_4 + -2*P3_6 + -2*P3_8 + 2*P4_5 + 2*P4_7 + -2*P5_6 + -2*P5_8 + 2*P6_7 + -2*P7_8
  F2 =~ -2*P1_3 + -2*P1_4 + -2*P1_7 + -2*P1_8 + -2*P2_3 + -2*P2_4 + -2*P2_7 + -2*P2_8 + 2*P3_5 + 2*P3_6 + 2*P4_5 + 2*P4_6 + -2*P5_7 + -2*P5_8 + -2*P6_7 + -2*P6_8
  F3 =~ -2*P1_5 + -2*P1_6 + -2*P1_7 + -2*P1_8 + -2*P2_5 + -2*P2_6 + -2*P2_7 + -2*P2_8 + -2*P3_5 + -2*P3_6 + -2*P3_7 + -2*P3_8 + -2*P4_5 + -2*P4_6 + -2*P4_7 + -2*P4_8

  P1_2 | 0*t1
  P1_3 | 0*t1
  P1_4 | 0*t1
  P1_5 | 0*t1
  P1_6 | 0*t1
  P1_7 | 0*t1
  P1_8 | 0*t1
  P2_3 | 0*t1
  P2_4 | 0*t1
  P2_5 | 0*t1
  P2_6 | 0*t1
  P2_7 | 0*t1
  P2_8 | 0*t1
  P3_4 | 0*t1
  P3_5 | 0*t1
  P3_6 | 0*t1
  P3_7 | 0*t1
  P3_8 | 0*t1
  P4_5 | 0*t1
  P4_6 | 0*t1
  P4_7 | 0*t1
  P4_8 | 0*t1
  P5_6 | 0*t1
  P5_7 | 0*t1
  P5_8 | 0*t1
  P6_7 | 0*t1
  P6_8 | 0*t1
  P7_8 | 0*t1

  F1 ~ 1
  F2 ~ 1
  F3 ~ 1

  F1 ~~ 1*F1 + F2 + F3
  F2 ~~ F2 + F3
  F3 ~~ F3

  P1_2 ~~ 0*P1_2
  P1_3 ~~ 0*P1_3
  P1_4 ~~ 0*P1_4
  P1_5 ~~ 0*P1_5
  P1_6 ~~ 0*P1_6
  P1_7 ~~ 0*P1_7
  P1_8 ~~ 0*P1_8
  P2_3 ~~ 0*P2_3
  P2_4 ~~ 0*P2_4
  P2_5 ~~ 0*P2_5
  P2_6 ~~ 0*P2_6
  P2_7 ~~ 0*P2_7
  P2_8 ~~ 0*P2_8
  P3_4 ~~ 0*P3_4
  P3_5 ~~ 0*P3_5
  P3_6 ~~ 0*P3_6
  P3_7 ~~ 0*P3_7
  P3_8 ~~ 0*P3_8
  P4_5 ~~ 0*P4_5
  P4_6 ~~ 0*P4_6
  P4_7 ~~ 0*P4_7
  P4_8 ~~ 0*P4_8
  P5_6 ~~ 0*P5_6
  P5_7 ~~ 0*P5_7
  P5_8 ~~ 0*P5_8
  P6_7 ~~ 0*P6_7
  P6_8 ~~ 0*P6_8
  P7_8 ~~ 0*P7_8

  U1 =~ 1*P1_2 + 1*P1_3 + 1*P1_4 + 1*P1_5 + 1*P1_6 + 1*P1_7 + 1*P1_8
  U2 =~ -1*P1_2 + 1*P2_3 + 1*P2_4 + 1*P2_5 + 1*P2_6 + 1*P2_7 + 1*P2_8
  U3 =~ -1*P1_3 + -1*P2_3 + 1*P3_4 + 1*P3_5 + 1*P3_6 + 1*P3_7 + 1*P3_8
  U4 =~ -1*P1_4 + -1*P2_4 + -1*P3_4 + 1*P4_5 + 1*P4_6 + 1*P4_7 + 1*P4_8
  U5 =~ -1*P1_5 + -1*P2_5 + -1*P3_5 + -1*P4_5 + 1*P5_6 + 1*P5_7 + 1*P5_8
  U6 =~ -1*P1_6 + -1*P2_6 + -1*P3_6 + -1*P4_6 + -1*P5_6 + 1*P6_7 + 1*P6_8
  U7 =~ -1*P1_7 + -1*P2_7 + -1*P3_7 + -1*P4_7 + -1*P5_7 + -1*P6_7 + 1*P7_8
  U8 =~ -1*P1_8 + -1*P2_8 + -1*P3_8 + -1*P4_8 + -1*P5_8 + -1*P6_8 + -1*P7_8

  U1 ~~ U1
  U2 ~~ U2
  U3 ~~ U3
  U4 ~~ U4
  U5 ~~ U5
  U6 ~~ U6
  U7 ~~ U7
  U8 ~~ U8
  U1 ~ 0*1
  U2 ~ 0*1
  U3 ~ 0*1
  U4 ~ 0*1
  U5 ~ 0*1
  U6 ~ 0*1
  U7 ~ 0*1
  U8 ~ 0*1
  U1 ~~ 0*F1
  U1 ~~ 0*F2
  U1 ~~ 0*F3
  U2 ~~ 0*F1
  U2 ~~ 0*F2
  U2 ~~ 0*F3
  U3 ~~ 0*F1
  U3 ~~ 0*F2
  U3 ~~ 0*F3
  U4 ~~ 0*F1
  U4 ~~ 0*F2
  U4 ~~ 0*F3
  U5 ~~ 0*F1
  U5 ~~ 0*F2
  U5 ~~ 0*F3
  U6 ~~ 0*F1
  U6 ~~ 0*F2
  U6 ~~ 0*F3
  U7 ~~ 0*F1
  U7 ~~ 0*F2
  U7 ~~ 0*F3
  U8 ~~ 0*F1
  U8 ~~ 0*F2
  U8 ~~ 0*F3
  U1 ~~ 0*U2
  U1 ~~ 0*U3
  U1 ~~ 0*U4
  U1 ~~ 0*U5
  U1 ~~ 0*U6
  U1 ~~ 0*U7
  U1 ~~ 0*U8
  U2 ~~ 0*U3
  U2 ~~ 0*U4
  U2 ~~ 0*U5
  U2 ~~ 0*U6
  U2 ~~ 0*U7
  U2 ~~ 0*U8
  U3 ~~ 0*U4
  U3 ~~ 0*U5
  U3 ~~ 0*U6
  U3 ~~ 0*U7
  U3 ~~ 0*U8
  U4 ~~ 0*U5
  U4 ~~ 0*U6
  U4 ~~ 0*U7
  U4 ~~ 0*U8
  U5 ~~ 0*U6
  U5 ~~ 0*U7
  U5 ~~ 0*U8
  U6 ~~ 0*U7
  U6 ~~ 0*U8
  U7 ~~ 0*U8
'



model_B <- '
  F1 =~ -2*P1_2 + -2*P1_4 + -2*P1_6 + -2*P1_8 + 2*P2_3 + 2*P2_5 + 2*P2_7 + -2*P3_4 + -2*P3_6 + -2*P3_8 + 2*P4_5 + 2*P4_7 + -2*P5_6 + -2*P5_8 + 2*P6_7 + -2*P7_8
  F2 =~ -2*P1_3 + -2*P1_4 + -2*P1_7 + -2*P1_8 + -2*P2_3 + -2*P2_4 + -2*P2_7 + -2*P2_8 + 2*P3_5 + 2*P3_6 + 2*P4_5 + 2*P4_6 + -2*P5_7 + -2*P5_8 + -2*P6_7 + -2*P6_8
  F3 =~ -2*P1_5 + -2*P1_6 + -2*P1_7 + -2*P1_8 + -2*P2_5 + -2*P2_6 + -2*P2_7 + -2*P2_8 + -2*P3_5 + -2*P3_6 + -2*P3_7 + -2*P3_8 + -2*P4_5 + -2*P4_6 + -2*P4_7 + -2*P4_8

  P1_2 | 0*t1
  P1_3 | 0*t1
  P1_4 | 0*t1
  P1_5 | 0*t1
  P1_6 | 0*t1
  P1_7 | 0*t1
  P1_8 | 0*t1
  P2_3 | 0*t1
  P2_4 | 0*t1
  P2_5 | 0*t1
  P2_6 | 0*t1
  P2_7 | 0*t1
  P2_8 | 0*t1
  P3_4 | 0*t1
  P3_5 | 0*t1
  P3_6 | 0*t1
  P3_7 | 0*t1
  P3_8 | 0*t1
  P4_5 | 0*t1
  P4_6 | 0*t1
  P4_7 | 0*t1
  P4_8 | 0*t1
  P5_6 | 0*t1
  P5_7 | 0*t1
  P5_8 | 0*t1
  P6_7 | 0*t1
  P6_8 | 0*t1
  P7_8 | 0*t1

  F1 ~ 1
  F2 ~ 1
  F3 ~ 1

  F1 ~~ 1*F1 + F2 + F3
  F2 ~~ F2 + F3
  F3 ~~ F3

  P1_2 ~~ P1_2
  P1_3 ~~ P1_3
  P1_4 ~~ P1_4
  P1_5 ~~ P1_5
  P1_6 ~~ P1_6
  P1_7 ~~ P1_7
  P1_8 ~~ P1_8
  P2_3 ~~ P2_3
  P2_4 ~~ P2_4
  P2_5 ~~ P2_5
  P2_6 ~~ P2_6
  P2_7 ~~ P2_7
  P2_8 ~~ P2_8
  P3_4 ~~ P3_4
  P3_5 ~~ P3_5
  P3_6 ~~ P3_6
  P3_7 ~~ P3_7
  P3_8 ~~ P3_8
  P4_5 ~~ P4_5
  P4_6 ~~ P4_6
  P4_7 ~~ P4_7
  P4_8 ~~ P4_8
  P5_6 ~~ P5_6
  P5_7 ~~ P5_7
  P5_8 ~~ P5_8
  P6_7 ~~ P6_7
  P6_8 ~~ P6_8
  P7_8 ~~ P7_8
'

# C: add eight independent product-level random errors, with signed loadings.
model_C <- '
  F1 =~ -2*P1_2 + -2*P1_4 + -2*P1_6 + -2*P1_8 + 2*P2_3 + 2*P2_5 + 2*P2_7 + -2*P3_4 + -2*P3_6 + -2*P3_8 + 2*P4_5 + 2*P4_7 + -2*P5_6 + -2*P5_8 + 2*P6_7 + -2*P7_8
  F2 =~ -2*P1_3 + -2*P1_4 + -2*P1_7 + -2*P1_8 + -2*P2_3 + -2*P2_4 + -2*P2_7 + -2*P2_8 + 2*P3_5 + 2*P3_6 + 2*P4_5 + 2*P4_6 + -2*P5_7 + -2*P5_8 + -2*P6_7 + -2*P6_8
  F3 =~ -2*P1_5 + -2*P1_6 + -2*P1_7 + -2*P1_8 + -2*P2_5 + -2*P2_6 + -2*P2_7 + -2*P2_8 + -2*P3_5 + -2*P3_6 + -2*P3_7 + -2*P3_8 + -2*P4_5 + -2*P4_6 + -2*P4_7 + -2*P4_8

  P1_2 | 0*t1
  P1_3 | 0*t1
  P1_4 | 0*t1
  P1_5 | 0*t1
  P1_6 | 0*t1
  P1_7 | 0*t1
  P1_8 | 0*t1
  P2_3 | 0*t1
  P2_4 | 0*t1
  P2_5 | 0*t1
  P2_6 | 0*t1
  P2_7 | 0*t1
  P2_8 | 0*t1
  P3_4 | 0*t1
  P3_5 | 0*t1
  P3_6 | 0*t1
  P3_7 | 0*t1
  P3_8 | 0*t1
  P4_5 | 0*t1
  P4_6 | 0*t1
  P4_7 | 0*t1
  P4_8 | 0*t1
  P5_6 | 0*t1
  P5_7 | 0*t1
  P5_8 | 0*t1
  P6_7 | 0*t1
  P6_8 | 0*t1
  P7_8 | 0*t1

  F1 ~ 1
  F2 ~ 1
  F3 ~ 1

  F1 ~~ 1*F1 + F2 + F3
  F2 ~~ F2 + F3
  F3 ~~ F3

  P1_2 ~~ P1_2
  P1_3 ~~ P1_3
  P1_4 ~~ P1_4
  P1_5 ~~ P1_5
  P1_6 ~~ P1_6
  P1_7 ~~ P1_7
  P1_8 ~~ P1_8
  P2_3 ~~ P2_3
  P2_4 ~~ P2_4
  P2_5 ~~ P2_5
  P2_6 ~~ P2_6
  P2_7 ~~ P2_7
  P2_8 ~~ P2_8
  P3_4 ~~ P3_4
  P3_5 ~~ P3_5
  P3_6 ~~ P3_6
  P3_7 ~~ P3_7
  P3_8 ~~ P3_8
  P4_5 ~~ P4_5
  P4_6 ~~ P4_6
  P4_7 ~~ P4_7
  P4_8 ~~ P4_8
  P5_6 ~~ P5_6
  P5_7 ~~ P5_7
  P5_8 ~~ P5_8
  P6_7 ~~ P6_7
  P6_8 ~~ P6_8
  P7_8 ~~ P7_8

  U1 =~ 1*P1_2 + 1*P1_3 + 1*P1_4 + 1*P1_5 + 1*P1_6 + 1*P1_7 + 1*P1_8
  U2 =~ -1*P1_2 + 1*P2_3 + 1*P2_4 + 1*P2_5 + 1*P2_6 + 1*P2_7 + 1*P2_8
  U3 =~ -1*P1_3 + -1*P2_3 + 1*P3_4 + 1*P3_5 + 1*P3_6 + 1*P3_7 + 1*P3_8
  U4 =~ -1*P1_4 + -1*P2_4 + -1*P3_4 + 1*P4_5 + 1*P4_6 + 1*P4_7 + 1*P4_8
  U5 =~ -1*P1_5 + -1*P2_5 + -1*P3_5 + -1*P4_5 + 1*P5_6 + 1*P5_7 + 1*P5_8
  U6 =~ -1*P1_6 + -1*P2_6 + -1*P3_6 + -1*P4_6 + -1*P5_6 + 1*P6_7 + 1*P6_8
  U7 =~ -1*P1_7 + -1*P2_7 + -1*P3_7 + -1*P4_7 + -1*P5_7 + -1*P6_7 + 1*P7_8
  U8 =~ -1*P1_8 + -1*P2_8 + -1*P3_8 + -1*P4_8 + -1*P5_8 + -1*P6_8 + -1*P7_8

  U1 ~~ U1
  U2 ~~ U2
  U3 ~~ U3
  U4 ~~ U4
  U5 ~~ U5
  U6 ~~ U6
  U7 ~~ U7
  U8 ~~ U8
  U1 ~ 0*1
  U2 ~ 0*1
  U3 ~ 0*1
  U4 ~ 0*1
  U5 ~ 0*1
  U6 ~ 0*1
  U7 ~ 0*1
  U8 ~ 0*1
  U1 ~~ 0*F1
  U1 ~~ 0*F2
  U1 ~~ 0*F3
  U2 ~~ 0*F1
  U2 ~~ 0*F2
  U2 ~~ 0*F3
  U3 ~~ 0*F1
  U3 ~~ 0*F2
  U3 ~~ 0*F3
  U4 ~~ 0*F1
  U4 ~~ 0*F2
  U4 ~~ 0*F3
  U5 ~~ 0*F1
  U5 ~~ 0*F2
  U5 ~~ 0*F3
  U6 ~~ 0*F1
  U6 ~~ 0*F2
  U6 ~~ 0*F3
  U7 ~~ 0*F1
  U7 ~~ 0*F2
  U7 ~~ 0*F3
  U8 ~~ 0*F1
  U8 ~~ 0*F2
  U8 ~~ 0*F3
  U1 ~~ 0*U2
  U1 ~~ 0*U3
  U1 ~~ 0*U4
  U1 ~~ 0*U5
  U1 ~~ 0*U6
  U1 ~~ 0*U7
  U1 ~~ 0*U8
  U2 ~~ 0*U3
  U2 ~~ 0*U4
  U2 ~~ 0*U5
  U2 ~~ 0*U6
  U2 ~~ 0*U7
  U2 ~~ 0*U8
  U3 ~~ 0*U4
  U3 ~~ 0*U5
  U3 ~~ 0*U6
  U3 ~~ 0*U7
  U3 ~~ 0*U8
  U4 ~~ 0*U5
  U4 ~~ 0*U6
  U4 ~~ 0*U7
  U4 ~~ 0*U8
  U5 ~~ 0*U6
  U5 ~~ 0*U7
  U5 ~~ 0*U8
  U6 ~~ 0*U7
  U6 ~~ 0*U8
  U7 ~~ 0*U8
'


# Fit each model directly. lavaan warnings remain visible.

# Model A
fit_A <- lavaan::sem(
  model = model_A, data = paired_demo, ordered = pn,
  estimator = "WLSMV", parameterization = "theta",
  meanstructure = TRUE, std.lv = FALSE,
  auto.fix.first = FALSE, auto.var = FALSE, auto.th = FALSE,
  auto.cov.lv.x = FALSE, auto.cov.y = FALSE
)
summary(fit_A, fit.measures = TRUE)

# Model B
fit_B <- lavaan::sem(
  model = model_B, data = paired_demo, ordered = pn,
  estimator = "WLSMV", parameterization = "theta",
  meanstructure = TRUE, std.lv = FALSE,
  auto.fix.first = FALSE, auto.var = FALSE, auto.th = FALSE,
  auto.cov.lv.x = FALSE, auto.cov.y = FALSE
)
summary(fit_B, fit.measures = TRUE)

# Model C
fit_C <- lavaan::sem(
  model = model_C, data = paired_demo, ordered = pn,
  estimator = "WLSMV", parameterization = "theta",
  meanstructure = TRUE, std.lv = FALSE,
  auto.fix.first = FALSE, auto.var = FALSE, auto.th = FALSE,
  auto.cov.lv.x = FALSE, auto.cov.y = FALSE
)
summary(fit_C, fit.measures = TRUE)

# Compare the three fits.
fit_indices <- c("chisq.scaled", "df.scaled", "pvalue.scaled",
                 "cfi.scaled", "rmsea.scaled", "srmr")
fit_table <- rbind(
  A = lavaan::fitMeasures(fit_A, fit_indices),
  B = lavaan::fitMeasures(fit_B, fit_indices),
  C = lavaan::fitMeasures(fit_C, fit_indices)
)
print(round(fit_table, 4))

# Keep fitted objects accessible to 02_TABLE.R.
results <- list(
  A = list(fit = fit_A, converged = lavaan::lavInspect(fit_A, "converged")),
  B = list(fit = fit_B, converged = lavaan::lavInspect(fit_B, "converged")),
  C = list(fit = fit_C, converged = lavaan::lavInspect(fit_C, "converged"))
)

scale_factor <- 1 / as.numeric(sim$truth$Phi["F1","F1"])
truth_normalized <- list(
  product_error_variances=setNames(sd_U^2*scale_factor,paste0("U",1:8)),
  pair_error_variances=setNames(sd_P^2*scale_factor,pn)
)
# Model C is retained without boundary adjustments. A negative variance is
# reported as a diagnostic and stops the analysis rather than silently changing
# the model. This also prevents stale C_fixed results from entering 02_TABLE.R.
if (!isTRUE(results$C$converged)) stop("Model C did not converge.")
pe_C <- lavaan::parameterEstimates(results$C$fit)
var_C <- subset(pe_C, op == "~~" & lhs == rhs)
if (any(!is.finite(var_C$est))) stop("Model C has non-finite variance estimates.")
neg_C <- subset(var_C, est < -1e-8)
cat("\n===== MODEL C VARIANCE DIAGNOSTICS =====\n")
cat("Negative variances:", nrow(neg_C), "\n")
cat("Minimum variance:", min(var_C$est), "\n")
if (nrow(neg_C) > 0L) {
  print(neg_C[, c("lhs", "est", "se")], row.names = FALSE)
  stop("Model C has negative variance estimates; inspect the unconstrained fit.")
}

# Preserve the structure consumed by 02_TABLE.R, except obsolete C_fixed fields.
# This is an original A/B/C analysis, not a boundary-corrected model.
results <- results[c("A", "B", "C")]
saveRDS(list(
  data = paired_demo,
  results = results,
  sim_truth = sim$truth,
  truth_normalized = truth_normalized,
  scale_factor = scale_factor,
  pn = pn,
  model_A = model_A,
  model_B = model_B,
  model_C = model_C
), "Demo02_results.rds")
cat("Saved Demo02_results.rds (A/B/C only)\n")
