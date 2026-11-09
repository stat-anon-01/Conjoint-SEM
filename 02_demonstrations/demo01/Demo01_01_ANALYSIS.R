# DEMONSTRATION 1: Binary residual variance and threshold restrictions
# Demo01_01_ANALYSIS.R — Fit the four models used in the article.
#
# From this folder, run:
#   source("Demo01_00_DGP.R")       # generate the simulated binary data
#   source("Demo01_01_ANALYSIS.R")  # fit M00, M10, M01, M11
#   source("Demo01_02_TABLE.R")     # summarize the fitted models
#
# Models (first digit: residual variance; second digit: B8 threshold):
#   M00: HOMO, equal thresholds; M10: HETERO, equal thresholds
#   M01: HOMO, B8 threshold free; M11: HETERO, B8 threshold free
# All use lavaan WLSMV with theta parameterization.
# This script assumes the DGP script has already created B in this session.
# ------------------------------------------------------------

library(lavaan)


# M00: HOMO; equal thresholds
model_M00 <- '
  F0 =~  1*B1 +  1*B2 +  1*B3 +  1*B4 +  1*B5 +  1*B6 +  1*B7 +  1*B8
  F1 =~ -1*B1 +  1*B2 + -1*B3 +  1*B4 + -1*B5 +  1*B6 + -1*B7 +  1*B8
  F2 =~ -1*B1 + -1*B2 +  1*B3 +  1*B4 + -1*B5 + -1*B6 +  1*B7 +  1*B8
  F3 =~ -1*B1 + -1*B2 + -1*B3 + -1*B4 +  1*B5 +  1*B6 +  1*B7 +  1*B8

  B1 | 0*t1
  B2 | 0*t1
  B3 | 0*t1
  B4 | 0*t1
  B5 | 0*t1
  B6 | 0*t1
  B7 | 0*t1
  B8 | 0*t1

  F0 ~ 1
  F1 ~ 1
  F2 ~ 1
  F3 ~ 1

  F0 ~~ F0 + F1 + F2 + F3
  F1 ~~ 1*F1 + F2 + F3
  F2 ~~ F2 + F3
  F3 ~~ F3

  B1 ~~ v*B1
  B2 ~~ v*B2
  B3 ~~ v*B3
  B4 ~~ v*B4
  B5 ~~ v*B5
  B6 ~~ v*B6
  B7 ~~ v*B7
  B8 ~~ v*B8
'

# M10: HETERO; equal thresholds
model_M10 <- '
  F0 =~  1*B1 +  1*B2 +  1*B3 +  1*B4 +  1*B5 +  1*B6 +  1*B7 +  1*B8
  F1 =~ -1*B1 +  1*B2 + -1*B3 +  1*B4 + -1*B5 +  1*B6 + -1*B7 +  1*B8
  F2 =~ -1*B1 + -1*B2 +  1*B3 +  1*B4 + -1*B5 + -1*B6 +  1*B7 +  1*B8
  F3 =~ -1*B1 + -1*B2 + -1*B3 + -1*B4 +  1*B5 +  1*B6 +  1*B7 +  1*B8

  B1 | 0*t1
  B2 | 0*t1
  B3 | 0*t1
  B4 | 0*t1
  B5 | 0*t1
  B6 | 0*t1
  B7 | 0*t1
  B8 | 0*t1

  F0 ~ 1
  F1 ~ 1
  F2 ~ 1
  F3 ~ 1

  F0 ~~ F0 + F1 + F2 + F3
  F1 ~~ 1*F1 + F2 + F3
  F2 ~~ F2 + F3
  F3 ~~ F3

  B1 ~~ B1
  B2 ~~ B2
  B3 ~~ B3
  B4 ~~ B4
  B5 ~~ B5
  B6 ~~ B6
  B7 ~~ B7
  B8 ~~ B8
'

# M01: HOMO; B8 threshold free
model_M01 <- '
  F0 =~  1*B1 +  1*B2 +  1*B3 +  1*B4 +  1*B5 +  1*B6 +  1*B7 +  1*B8
  F1 =~ -1*B1 +  1*B2 + -1*B3 +  1*B4 + -1*B5 +  1*B6 + -1*B7 +  1*B8
  F2 =~ -1*B1 + -1*B2 +  1*B3 +  1*B4 + -1*B5 + -1*B6 +  1*B7 +  1*B8
  F3 =~ -1*B1 + -1*B2 + -1*B3 + -1*B4 +  1*B5 +  1*B6 +  1*B7 +  1*B8

  B1 | 0*t1
  B2 | 0*t1
  B3 | 0*t1
  B4 | 0*t1
  B5 | 0*t1
  B6 | 0*t1
  B7 | 0*t1
  B8 | tau8*t1

  F0 ~ 1
  F1 ~ 1
  F2 ~ 1
  F3 ~ 1

  F0 ~~ F0 + F1 + F2 + F3
  F1 ~~ 1*F1 + F2 + F3
  F2 ~~ F2 + F3
  F3 ~~ F3

  B1 ~~ v*B1
  B2 ~~ v*B2
  B3 ~~ v*B3
  B4 ~~ v*B4
  B5 ~~ v*B5
  B6 ~~ v*B6
  B7 ~~ v*B7
  B8 ~~ v*B8
'

# M11: HETERO; B8 threshold free
model_M11 <- '
  F0 =~  1*B1 +  1*B2 +  1*B3 +  1*B4 +  1*B5 +  1*B6 +  1*B7 +  1*B8
  F1 =~ -1*B1 +  1*B2 + -1*B3 +  1*B4 + -1*B5 +  1*B6 + -1*B7 +  1*B8
  F2 =~ -1*B1 + -1*B2 +  1*B3 +  1*B4 + -1*B5 + -1*B6 +  1*B7 +  1*B8
  F3 =~ -1*B1 + -1*B2 + -1*B3 + -1*B4 +  1*B5 +  1*B6 +  1*B7 +  1*B8

  B1 | 0*t1
  B2 | 0*t1
  B3 | 0*t1
  B4 | 0*t1
  B5 | 0*t1
  B6 | 0*t1
  B7 | 0*t1
  B8 | tau8*t1

  F0 ~ 1
  F1 ~ 1
  F2 ~ 1
  F3 ~ 1

  F0 ~~ F0 + F1 + F2 + F3
  F1 ~~ 1*F1 + F2 + F3
  F2 ~~ F2 + F3
  F3 ~~ F3

  B1 ~~ B1
  B2 ~~ B2
  B3 ~~ B3
  B4 ~~ B4
  B5 ~~ B5
  B6 ~~ B6
  B7 ~~ B7
  B8 ~~ B8
'

fit_M00 <- sem(model_M00, data = B, ordered = names(B),
  estimator = 'WLSMV', parameterization = 'theta',
  meanstructure = TRUE, std.lv = FALSE,
  auto.fix.first = FALSE, auto.var = FALSE, auto.th = FALSE,
  auto.cov.lv.x = FALSE, auto.cov.y = FALSE)

fit_M10 <- sem(model_M10, data = B, ordered = names(B),
  estimator = 'WLSMV', parameterization = 'theta',
  meanstructure = TRUE, std.lv = FALSE,
  auto.fix.first = FALSE, auto.var = FALSE, auto.th = FALSE,
  auto.cov.lv.x = FALSE, auto.cov.y = FALSE)

fit_M01 <- sem(model_M01, data = B, ordered = names(B),
  estimator = 'WLSMV', parameterization = 'theta',
  meanstructure = TRUE, std.lv = FALSE,
  auto.fix.first = FALSE, auto.var = FALSE, auto.th = FALSE,
  auto.cov.lv.x = FALSE, auto.cov.y = FALSE)

fit_M11 <- sem(model_M11, data = B, ordered = names(B),
  estimator = 'WLSMV', parameterization = 'theta',
  meanstructure = TRUE, std.lv = FALSE,
  auto.fix.first = FALSE, auto.var = FALSE, auto.th = FALSE,
  auto.cov.lv.x = FALSE, auto.cov.y = FALSE)

# Model fit
fitMeasures(fit_M00, c('chisq.scaled','df.scaled','cfi.scaled','rmsea.scaled','srmr'))
fitMeasures(fit_M10, c('chisq.scaled','df.scaled','cfi.scaled','rmsea.scaled','srmr'))
fitMeasures(fit_M01, c('chisq.scaled','df.scaled','cfi.scaled','rmsea.scaled','srmr'))
fitMeasures(fit_M11, c('chisq.scaled','df.scaled','cfi.scaled','rmsea.scaled','srmr'))

# Confirmatory nested comparisons (WLSMV)
lavTestLRT(fit_M00, fit_M10)
lavTestLRT(fit_M01, fit_M11)
lavTestLRT(fit_M00, fit_M01)
lavTestLRT(fit_M10, fit_M11)

# Estimates
parameterEstimates(fit_M11)
