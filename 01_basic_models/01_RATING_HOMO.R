# ================================================================
# Supplementary Example S1: Conventional rating-based conjoint SEM
# ================================================================
# Requirements: R packages lavaan and MASS.
# Place this script in the same folder as 00_DGP.R.
# Run: source("01_RATING_HOMO.R")
#
# Eight product profiles are formed by three two-level attributes.
# Each participant rates all eight profiles. Effect coding (-1, +1)
# defines three individual-specific part-worth factors (F1-F3).
# F0 is the participant-specific overall evaluation (baseline).
#
# Conventional reference assumptions illustrated here:
#   (a) additive preference function (no attribute interactions);
#   (b) no additional profile-specific intercepts;
#   (c) independent errors with a common residual variance.
#
# The supplied 00_DGP.R generates independent rating
# residuals with a common variance of 1, matching this model.
# Binary and paired-comparison residuals are generated with unit SD;
# they are not used in this rating example.
# ================================================================

if (!requireNamespace("lavaan", quietly = TRUE)) {
  stop('Install lavaan with install.packages("lavaan")')
}
library(lavaan)

source("00_DGP.R")

rating_data <- sim$rating

# The design matrix is:
#   Profile  X1  X2  X3
#      1    -1  -1  -1
#      2    +1  -1  -1
#      3    -1  +1  -1
#      4    +1  +1  -1
#      5    -1  -1  +1
#      6    +1  -1  +1
#      7    -1  +1  +1
#      8    +1  +1  +1
print(sim$X)

# All factor loadings are fixed by the experimental design.
# Rating intercepts are fixed to zero; factor means are estimated.
# This is a location specification, not an assertion that mean ratings
# are zero. F0-F3 variances/covariances are freely estimated.
#
# The repeated label 'theta' constrains all eight residual variances
# to the same positive value. No residual covariance is specified;
# auto.cov.y = FALSE prevents lavaan from adding any automatically.
model_rating_conventional <- '
  # Part-worth representation: U_ij = F0_i + x_j1 F1_i + x_j2 F2_i + x_j3 F3_i
  F0 =~     1*R1 +     1*R2 +     1*R3 +     1*R4 +     1*R5 +     1*R6 +     1*R7 +     1*R8
  F1 =~ (-1)*R1 +     1*R2 + (-1)*R3 +     1*R4 + (-1)*R5 +     1*R6 + (-1)*R7 +     1*R8
  F2 =~ (-1)*R1 + (-1)*R2 +     1*R3 +     1*R4 + (-1)*R5 + (-1)*R6 +     1*R7 +     1*R8
  F3 =~ (-1)*R1 + (-1)*R2 + (-1)*R3 + (-1)*R4 +     1*R5 +     1*R6 +     1*R7 +     1*R8

  # No additional profile-specific intercepts
  R1 ~ 0*1
  R2 ~ 0*1
  R3 ~ 0*1
  R4 ~ 0*1
  R5 ~ 0*1
  R6 ~ 0*1
  R7 ~ 0*1
  R8 ~ 0*1

  # Individual differences in baseline and attribute part-worth means
  F0 ~ 1
  F1 ~ 1
  F2 ~ 1
  F3 ~ 1

  # Unrestricted latent preference variance-covariance matrix
  F0 ~~ F0 + F1 + F2 + F3
  F1 ~~ F1 + F2 + F3
  F2 ~~ F2 + F3
  F3 ~~ F3

  # Conventional assumption: equal, independent residual errors
  R1 ~~ theta*R1
  R2 ~~ theta*R2
  R3 ~~ theta*R3
  R4 ~~ theta*R4
  R5 ~~ theta*R5
  R6 ~~ theta*R6
  R7 ~~ theta*R7
  R8 ~~ theta*R8
'

fit_rating_conventional <- lavaan::sem(
  model = model_rating_conventional,
  data = rating_data,
  estimator = "ML",
  meanstructure = TRUE,
  std.lv = FALSE,
  auto.fix.first = FALSE,
  auto.var = FALSE,
  auto.cov.lv.x = FALSE,
  auto.cov.y = FALSE
)

summary(fit_rating_conventional, fit.measures = TRUE)



