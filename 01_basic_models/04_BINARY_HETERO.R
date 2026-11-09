# ================================================================
# Supplementary Example S2: Conventional binary conjoint SEM
# ================================================================
# Requirements: lavaan and MASS.
# Place this script beside 00_DGP.R.
# Run: source("04_BINARY_HETERO.R")
#
# Eight profiles combine three two-level attributes (effect coded
# -1/+1). Each respondent accepts (1) or rejects (0) each profile.
# F0 is an individual-specific baseline; F1-F3 are individual-specific
# part-worths. All loadings are fixed by the experimental design.
#
# Model: independent heterogeneous residual variances.
# Common data and design are shared with the HOMO example.

if (!requireNamespace("lavaan", quietly = TRUE)) {
  stop('Install lavaan: install.packages("lavaan")')
}
library(lavaan)

source("00_DGP.R")

binary_data <- sim$binary


# Binary latent response:
#   B*_ij = F0_i + x_j1 F1_i + x_j2 F2_i + x_j3 F3_i + e_ij
#   B_ij  = 1(B*_ij > 0)
# Here e_ij has unit variance. Fixing thresholds at zero absorbs
# the common DGP threshold into the mean of F0.
#
# Factor variances and covariances are estimated freely. Unlike the
# original heterogeneous model, this specification fixes all eight
# latent-response residual variances at 1. It therefore does not
# require an additional F1 variance anchor.
# To extend this model to heterogeneous residual variances,
# replace the explicit fixed residual lines below with free variances,
# and set Var(F0) = 1 as an alternative scale anchor.
# The alternative is shown as commented lavaan syntax.

# ============================================================
# ADDITIONAL EXAMPLE: HETEROSCEDASTIC MODEL
# ============================================================
# This is a standalone HETERO model; no HOMO fit is required.
# This example fits the SAME data, generated with HOMO errors.
# To simulate HETERO errors, set sd_R, sd_B, or sd_P in DGP.
# Categorical latent-response scale: fix Var(F0) = 1
# while freeing individual residual variances (theta parameterization).
# This changes the metric relative to the HOMO fit. Do not compare
# unstandardized parameter estimates across the two scales directly.
# The fixed factor variance is a scale anchor, not a substantive test.

model_binary_hetero <- '
  # Experimental-design loadings
  F0 =~     1*B1 +     1*B2 +     1*B3 +     1*B4 +     1*B5 +     1*B6 +     1*B7 +     1*B8
  F1 =~ (-1)*B1 +     1*B2 + (-1)*B3 +     1*B4 + (-1)*B5 +     1*B6 + (-1)*B7 +     1*B8
  F2 =~ (-1)*B1 + (-1)*B2 +     1*B3 +     1*B4 + (-1)*B5 + (-1)*B6 +     1*B7 +     1*B8
  F3 =~ (-1)*B1 + (-1)*B2 + (-1)*B3 + (-1)*B4 +     1*B5 +     1*B6 +     1*B7 +     1*B8

  # Common threshold, represented by zero after location normalization
  B1 | 0*t1
  B2 | 0*t1
  B3 | 0*t1
  B4 | 0*t1
  B5 | 0*t1
  B6 | 0*t1
  B7 | 0*t1
  B8 | 0*t1

  # Baseline and part-worth means
  F0 ~ 1
  F1 ~ 1
  F2 ~ 1
  F3 ~ 1

  # Unrestricted individual-difference covariance matrix
  F0 ~~ 1*F0 + F1 + F2 + F3
  F1 ~~ F1 + F2 + F3
  F2 ~~ F2 + F3
  F3 ~~ F3

  # Heterogeneous: independent, freely estimated residual variances
  B1 ~~ B1
  B2 ~~ B2
  B3 ~~ B3
  B4 ~~ B4
  B5 ~~ B5
  B6 ~~ B6
  B7 ~~ B7
  B8 ~~ B8
'

fit_binary_hetero <- lavaan::sem(
  model = model_binary_hetero,
  data = binary_data,
  ordered = names(binary_data),
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

summary(fit_binary_hetero, fit.measures = TRUE)

# Scale interpretation:
# The HOMO model fixes each residual variance at 1; this HETERO
# model instead fixes Var(F0) = 1. In the shared DGP,
# Var(F0) = 0.25 and residual variances are 1. Under the
# HETERO scale normalization, the corresponding residual variance
# is therefore 1 / 0.25 = 4, not 1. The estimated residual
# variances should be interpreted relative to this scale.



