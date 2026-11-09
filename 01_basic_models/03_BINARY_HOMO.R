# ================================================================
# Supplementary Example S2: Conventional binary conjoint SEM
# ================================================================
# Requirements: lavaan and MASS.
# Place this script beside 00_DGP.R.
# Run: source("03_BINARY_HOMO.R")
#
# Eight profiles combine three two-level attributes (effect coded
# -1/+1). Each respondent accepts (1) or rejects (0) each profile.
# F0 is an individual-specific baseline; F1-F3 are individual-specific
# part-worths. All loadings are fixed by the experimental design.
#
# Reference assumptions:
#   (a) additive preference function (no interactions);
#   (b) a common response threshold across profiles;
#   (c) independent, equal-scale latent-response errors.
#
# Binary outcomes are generated with a common DGP threshold of 0.5
# and independent residuals of SD 1. In estimation, thresholds are
# fixed to zero as a LOCATION normalization. Consequently, the
# estimated F0 mean corresponds to the DGP baseline minus 0.5.
# The binary latent-response SCALE is normalized by fixing the
# residual variance to 1 (probit scale). This is not an empirical
# claim that the binary and rating scales are identical.
# ================================================================

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
model_binary_conventional <- '
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
  F0 ~~ F0 + F1 + F2 + F3
  F1 ~~ F1 + F2 + F3
  F2 ~~ F2 + F3
  F3 ~~ F3

  # Conventional: independent unit-variance residuals
  B1 ~~ 1*B1
  B2 ~~ 1*B2
  B3 ~~ 1*B3
  B4 ~~ 1*B4
  B5 ~~ 1*B5
  B6 ~~ 1*B6
  B7 ~~ 1*B7
  B8 ~~ 1*B8

  # Heterogeneous alternative: replace the eight lines above with:
  # B1 ~~ B1
  # B2 ~~ B2
  # B3 ~~ B3
  # B4 ~~ B4
  # B5 ~~ B5
  # B6 ~~ B6
  # B7 ~~ B7
  # B8 ~~ B8
  # Also replace F0 ~~ F0 + F1 + F2 + F3 with
  # F0 ~~ 1*F0 + F1 + F2 + F3
'


fit_binary_conventional <- lavaan::sem(
  model = model_binary_conventional,
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


summary(fit_binary_conventional, fit.measures = TRUE)


