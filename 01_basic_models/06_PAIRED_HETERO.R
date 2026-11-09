# ================================================================
# Supplementary Example S3: Conventional paired-comparison conjoint SEM
# ================================================================
# Requirements: lavaan and MASS.
# Place this script beside 00_DGP.R.
# Run: source("06_PAIRED_HETERO.R")
#
# Eight profiles combine three two-level attributes (effect coded
# -1/+1). Each respondent makes all 28 pairwise choices. Pj_k = 1
# indicates that profile j is preferred over profile k (j < k).
#
# Reference assumptions:
#   (a) additive part-worth utilities without attribute interactions;
#   (b) zero pairwise thresholds (no systematic position bias);
#   (c) independent, equal-variance comparison-specific errors.
#
# Latent preference difference:
#   P*_(i,jk) = (x_j - x_k)' (F1_i, F2_i, F3_i) + e_(i,jk)
#   P_(i,jk)  = 1(P*_(i,jk) > 0)
#
# The baseline F0 cancels from each difference and is not estimated.
# Each comparison error is independent with variance fixed at 1,
# defining the probit response scale. No rating-model anchor is used.
# Factor means, variances, and covariances are estimated freely
# in the conventional unit-residual specification.
# The 28 loading vectors are the differences between profile rows
# in the shared effect-coded design matrix.
# ================================================================

if (!requireNamespace("lavaan", quietly = TRUE)) {
  stop('Install lavaan: install.packages("lavaan")')
}

library(lavaan)

source("00_DGP.R")

paired_data <- sim$paired


# The factor loadings below are deliberately written explicitly so
# readers can see the full lavaan syntax and adapt it to new designs.
# Zero loadings are omitted; the original design matrix is sim$D.
# To extend this model to heterogeneous residual variances,
# replace the explicit fixed residual lines below with free variances,
# and set Var(F1) = 1 as an alternative scale anchor.
# The alternative is shown as commented lavaan syntax.

# ============================================================
# ADDITIONAL EXAMPLE: HETEROSCEDASTIC MODEL
# ============================================================
# This is a standalone HETERO model; no HOMO fit is required.
# This example fits the SAME data, generated with HOMO errors.
# To simulate HETERO errors, set sd_R, sd_B, or sd_P in DGP.
# Categorical latent-response scale: fix Var(F1) = 1
# while freeing individual residual variances (theta parameterization).
# This changes the metric relative to the HOMO fit. Do not compare
# unstandardized parameter estimates across the two scales directly.
# The fixed factor variance is a scale anchor, not a substantive test.

model_paired_hetero <- '

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

  # Heterogeneous independent comparison errors
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

fit_paired_hetero <- lavaan::sem(
  model = model_paired_hetero,
  data = paired_data,
  ordered = names(paired_data),
  estimator = "WLSMV",
  parameterization = "theta",
  meanstructure = TRUE,
  std.lv = FALSE,
  auto.fix.first = FALSE,
  auto.var = FALSE,
  auto.th = FALSE,
  auto.cov.lv.x = FALSE,
  auto.cov.y = FALSE
)

summary(fit_paired_hetero, fit.measures = TRUE)

# Scale interpretation:
# The HOMO model fixes each residual variance at 1; this HETERO
# model instead fixes Var(F1) = 1. In the shared DGP,
# Var(F1) = 0.25 and residual variances are 1. Under the
# HETERO scale normalization, the corresponding residual variance
# is therefore 1 / 0.25 = 4, not 1. The estimated residual
# variances should be interpreted relative to this scale.

