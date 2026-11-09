# DEMO 1 / 02_TABLE.R
# Run after Demo01_00_DGP.R and Demo01_01_ANALYSIS.R.

library(lavaan)

fits <- list(M00 = fit_M00, M10 = fit_M10,
             M01 = fit_M01, M11 = fit_M11)

# Table 1: Model fit
indices <- c('chisq.scaled', 'df.scaled', 'cfi.scaled',
             'rmsea.scaled', 'srmr')
Table1 <- data.frame(
  Model = c('M00: HOMO / Equal', 'M10: HETERO / Equal',
            'M01: HOMO / B8 free', 'M11: HETERO / B8 free'),
  do.call(rbind, lapply(fits, fitMeasures, fit.measures = indices)),
  row.names = NULL, check.names = FALSE
)
print(Table1, row.names = FALSE)

# Table 1A: WLSMV scaled chi-square difference tests
lrt_00_10 <- lavTestLRT(fit_M00, fit_M10)
lrt_01_11 <- lavTestLRT(fit_M01, fit_M11)
lrt_00_01 <- lavTestLRT(fit_M00, fit_M01)
lrt_10_11 <- lavTestLRT(fit_M10, fit_M11)

Table1A <- data.frame(
  Comparison = c('M00 vs M10', 'M01 vs M11',
                 'M00 vs M01', 'M10 vs M11'),
  Test = c('Residual equality, thresholds equal',
           'Residual equality, B8 threshold free',
           'B8 threshold, homogeneous residuals',
           'B8 threshold, heterogeneous residuals'),
  Delta_Chisq = c(lrt_00_10[2, 'Chisq diff'],
                  lrt_01_11[2, 'Chisq diff'],
                  lrt_00_01[2, 'Chisq diff'],
                  lrt_10_11[2, 'Chisq diff']),
  Delta_df = c(lrt_00_10[2, 'Df diff'],
               lrt_01_11[2, 'Df diff'],
               lrt_00_01[2, 'Df diff'],
               lrt_10_11[2, 'Df diff']),
  p = c(lrt_00_10[2, 'Pr(>Chisq)'],
        lrt_01_11[2, 'Pr(>Chisq)'],
        lrt_00_01[2, 'Pr(>Chisq)'],
        lrt_10_11[2, 'Pr(>Chisq)'])
)
print(Table1A, row.names = FALSE)

# Table 2: Latent factor means with model-based standard errors
# DGP: Var(F1) = .25; fitted Var(F1) = 1, hence scale factor = 2.
# The common DGP threshold (.50) is absorbed into the F0 mean.
# DGP targets: F0 = 2*(.20-.50) = -.60; F1-F3 = 2*.20 = .40.
# Estimates and SEs below are on the fitted (identified) scale.

factors <- paste0('F', 0:3)
DGP_target <- c(-.60, .40, .40, .40)

# Extract estimates and SEs by parameter name, rather than row order.
get_factor_means <- function(fit) {
  pe <- parameterEstimates(fit)
  rows <- pe[pe$op == '~1' & pe$lhs %in% factors,
             c('lhs', 'est', 'se')]
  idx <- match(factors, rows$lhs)
  if (anyNA(idx)) stop('Missing factor means in parameterEstimates().')
  rows[idx, c('est', 'se')]
}

means <- lapply(fits, get_factor_means)

# Numeric version is retained for further analysis or export.
Table2_numeric <- data.frame(Factor = factors, DGP_target = DGP_target)
for (m in names(means)) {
  Table2_numeric[[paste0(m, '_est')]] <- means[[m]]$est
  Table2_numeric[[paste0(m, '_SE')]] <- means[[m]]$se
}

format_est_se <- function(est, se, digits = 3) {
  ifelse(is.na(se),
         sprintf(paste0('%.', digits, 'f (NA)'), est),
         sprintf(paste0('%.', digits, 'f (%.', digits, 'f)'), est, se))
}

Table2 <- data.frame(
  Factor = factors,
  DGP_target = sprintf('%.3f', DGP_target),
  M00 = format_est_se(means$M00$est, means$M00$se),
  M10 = format_est_se(means$M10$est, means$M10$se),
  M01 = format_est_se(means$M01$est, means$M01$se),
  M11 = format_est_se(means$M11$est, means$M11$se),
  check.names = FALSE
)
cat('\nTable 2. Factor means: Estimate (SE)\n')
print(Table2, row.names = FALSE, right = FALSE)
cat('Note. Parentheses contain lavaan model-based standard errors.\n')

# Supplementary inspection: latent variances and covariances.
# Align by both lhs and rhs, avoiding accidental row-order mismatches.
get_latent_cov <- function(fit) {
  pe <- parameterEstimates(fit)
  pe <- pe[pe$op == '~~' & pe$lhs %in% factors & pe$rhs %in% factors,
           c('lhs', 'rhs', 'est', 'se')]
  pe$Parameter <- paste(pe$lhs, pe$rhs, sep = ' ~~ ')
  pe
}
latent <- lapply(fits, get_latent_cov)
parameters <- latent$M11$Parameter
TableS1 <- data.frame(Parameter = parameters, check.names = FALSE)
for (m in names(latent)) {
  idx <- match(parameters, latent[[m]]$Parameter)
  TableS1[[m]] <- latent[[m]]$est[idx]
}
cat('\nSupplementary Table S1. Latent variances and covariances\n')
print(TableS1, row.names = FALSE)
