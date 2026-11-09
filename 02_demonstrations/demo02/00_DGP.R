# DEMONSTRATION 2: Residual dependence in paired comparisons
# 00_DGP.R — Data generation only
#
# Generates 1,000 respondents, eight effect-coded profiles and 28
# binary paired comparisons. Model C is the generating mechanism:
# profile-level disturbances plus independent comparison-specific errors.
# The intermediate rating/binary data are retained to preserve the original
# random-number sequence; only paired_demo is fitted in 01_ANALYSIS.R.
# Do not change the seeds if reproducing the published tables.
#
# Run 01_ANALYSIS.R from this directory; it sources this file.
# ------------------------------------------------------------

if (!requireNamespace('MASS', quietly = TRUE)) {
  stop('Install MASS: install.packages("MASS")')
}

make_study1_data <- function(N = 1000L, seed = 12345L) {
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

  # HETERO: each observed response has its own residual SD.
  # The same latent F generates R, B, and P; no format-specific factors.
  sd_R <- c(.80, .80, .60, .60, 1, 1, 1.20, 1.20)
  sd_B <- c(.80, .80, .60, .60, 1, 1, 1.20, 1.20)
  tau_B <- rep(.50, 8)  # DGP threshold; estimation may set tau=0 as location normalization

  pairs <- t(utils::combn(1:8, 2))  # 28 x 2; first profile preferred if 1
  pair_names <- paste0('P', pairs[, 1], '_', pairs[, 2])
  sd_P <- c(
    .80, 1.00, .70, 1.20, .90, 1.10, .60,
    1.30, .75, 1.15, .85, 1.25, .95,
    .65, 1.05, .90, 1.20, .80,
    1.10, .70, 1.30, .85,
    1.00, .60, 1.15,
    .75, 1.25,
    .95
  )
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

# Run once for Study 1. Change N/seed here only.
sim <- make_study1_data(N = 1000L, seed = 12345L)
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

stopifnot(ncol(rating) == 8L, ncol(binary) == 8L,
          ncol(paired) == 28L, all(D[, 'F0'] == 0),
          all(unlist(binary) %in% 0:1),
          all(unlist(paired) %in% 0:1))

cat('Intermediate DGP generated (Binary tau =', truth$threshold_B[1], '):', N,
    'respondents; 8R + 8B + 28P.\n')

# Demo02-specific paired comparison DGP (Model C).
set.seed(20261103)
N <- sim$N
pairs <- sim$pairs
pn <- colnames(sim$paired)
sd_P <- unname(sim$truth$sd_P)
sd_U <- c(.80,.90,.70,1.00,.75,.95,.85,1.10)

G <- matrix(0,28,8,dimnames=list(pn,paste0("U",1:8)))
for (r in seq_len(28)) {
  G[r,pairs[r,1]] <- 1
  G[r,pairs[r,2]] <- -1
}
stopifnot(qr(G)$rank == 7L, identical(rownames(G),pn))

U <- sweep(matrix(rnorm(N*8),N,8),2,sd_U,"*")
E <- sweep(matrix(rnorm(N*28),N,28),2,sd_P,"*")
base_util <- sim$F %*% t(sim$D)
paired_demo <- as.data.frame(1L * (base_util + U %*% t(G) + E > 0))
names(paired_demo) <- pn
stopifnot(ncol(paired_demo)==28L, all(unlist(paired_demo) %in% 0:1))

# Population covariance decomposition, not estimated quantities.
Theta_A <- G %*% diag(sd_U^2,8) %*% t(G)
Theta_B <- diag(sd_P^2,28)
Theta_C <- Theta_A + Theta_B
cat("Population residual covariance ranks: ",
    paste(c(A=qr(Theta_A)$rank,B=qr(Theta_B)$rank,
            C=qr(Theta_C)$rank),collapse=" / "),"\n")


