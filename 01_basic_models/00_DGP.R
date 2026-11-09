# ============================================================
# Supplementary S1 | 00_DGP_CONVENTIONAL.R | Rating, binary, and paired-comparison reference DGP
# Rating (8), Binary (8), Paired comparison (28)
# Fixed conjoint design loadings; independent residuals by response.
# Rating residual variances are equal, matching 01_RATING_HOMO.R.
# Binary and paired errors default to independent unit variances;
# positive SD vectors allow heterogeneous simulation scenarios.
# ============================================================

if (!requireNamespace('MASS', quietly = TRUE)) {
  stop('Install MASS: install.packages("MASS")')
}

make_data <- function(N = 1000L, seed = 12345L,
                             sd_R = rep(1, 8),
                             sd_B = rep(1, 8),
                             sd_P = rep(1, 28)) {
  stopifnot(length(N) == 1L, N > 4L)
  set.seed(seed)

  design <- as.matrix(expand.grid(
    X1 = c(-1, 1), X2 = c(-1, 1), X3 = c(-1, 1)
  ))
  rownames(design) <- paste0('Profile', 1:8)
  X <- cbind(F0 = 1, design)

  mu <- c(F0 = .20, F1 = .20, F2 = .20, F3 = .20)
  Phi <- matrix(c(
    1, .30, .20, .10,
    .30, 1, .25, .20,
    .20, .25, 1.50, .30,
    .10, .20, .30, 2
  ), nrow = 4, byrow = TRUE) * .25
  dimnames(Phi) <- list(names(mu), names(mu))
  stopifnot(min(eigen(Phi, symmetric = TRUE)$values) > 0)

  F <- MASS::mvrnorm(n = N, mu = mu, Sigma = Phi)
  colnames(F) <- names(mu)
  utility <- F %*% t(X)

  # RATING: homoscedastic, independent residuals (Var = 1).
  # BINARY: independent unit-variance errors. PAIRED: independent unit-variance errors.
  # The same latent F generates R, B, and P; no format-specific factors.
  # For heteroscedastic extensions, supply positive SD vectors of
  # lengths 8, 8, and 28. Defaults reproduce the conventional DGP.
  stopifnot(length(sd_R) == 8L, length(sd_B) == 8L,
            length(sd_P) == 28L, all(is.finite(c(sd_R, sd_B, sd_P))),
            all(c(sd_R, sd_B, sd_P) > 0))
  tau_B <- rep(.50, 8)  # DGP threshold; estimation may set tau=0 as location normalization

  pairs <- t(utils::combn(1:8, 2))  # 28 x 2; first profile preferred if 1
  pair_names <- paste0('P', pairs[, 1], '_', pairs[, 2])
  # PAIRED: 28 independent comparison-specific errors with common SD 1.
  # This matches the conventional paired-comparison reference model.
  names(sd_R) <- paste0('R', 1:8)
  names(sd_B) <- paste0('B', 1:8)
  names(sd_P) <- pair_names
  names(tau_B) <- names(sd_B)

  e_R <- sweep(matrix(rnorm(N * 8), N, 8), 2, sd_R, '*')
  e_B <- sweep(matrix(rnorm(N * 8), N, 8), 2, sd_B, '*')
  e_P <- sweep(matrix(rnorm(N * 28), N, 28), 2, sd_P, '*')

  rating <- as.data.frame(utility + e_R)
  names(rating) <- names(sd_R)
  binary <- as.data.frame(1L * (sweep(utility + e_B, 2, tau_B, '-') > 0))
  names(binary) <- names(sd_B)

  # IMPORTANT: each pair has an independent comparison-level error.
  # Not constructed by subtracting profile-level errors.
  D <- X[pairs[, 1], , drop = FALSE] - X[pairs[, 2], , drop = FALSE]
  rownames(D) <- pair_names
  # D[, 'F0'] == 0: paired responses cannot identify F0.
  pair_utility <- F %*% t(D)
  paired <- as.data.frame(1L * (pair_utility + e_P > 0))
  names(paired) <- pair_names

  truth <- list(
    mu = mu, Phi = Phi, cor = cov2cor(Phi),
    residual_var_R = sd_R^2,
    residual_var_B = sd_B^2,
    residual_var_P = sd_P^2,
    sd_R = sd_R, sd_B = sd_B, sd_P = sd_P,
    threshold_B = tau_B,
    # When fitting B with thresholds fixed to zero, F0 absorbs -tau.
    mu_B_tau0 = setNames(c(mu["F0"] - tau_B[1], unname(mu[-1])), names(mu)),
    threshold_P = setNames(rep(0, 28), pair_names)
  )

  list(N = N, seed = seed, X = X, design = design, pairs = pairs,
       D = D, F = F, utility = utility, rating = rating,
       binary = binary, paired = paired, truth = truth,
       data = data.frame(rating, binary, paired, check.names = FALSE))
}

# Run once. Defaults: equal residual SDs of 1.
# To simulate heteroscedastic errors, replace sd_B / sd_P below
# with positive vectors of lengths 8 / 28, respectively.
# Keep model specification and DGP assumptions conceptually distinct.
sim <- make_data(
  N = 1000L, seed = 12345L,
  sd_R = rep(1, 8),
  sd_B = rep(1, 8),
  sd_P = rep(1, 28)
)

N <- sim$N
X <- sim$X
D <- sim$D
pairs <- sim$pairs
rating <- sim$rating
binary <- sim$binary
paired <- sim$paired
truth <- sim$truth
mu_pop <- truth$mu
Phi_pop <- truth$Phi


