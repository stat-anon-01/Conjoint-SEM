# Demo01 / 00_DGP.R
# Recovered verbatim DGP logic from Demo01_Binary_Threshold.R.
# All design, mu, Phi, residual SD, seed, N, delta and RNG order preserved.
# Sourcing this file creates B, tau_B, delta, X, mu, Phi, F.

if (!requireNamespace('MASS', quietly=TRUE)) stop('Install MASS first.')
set.seed(20261017)
N <- 1000L
delta <- 0.80  # B8-specific threshold deviation

# ------------------------------------------------------------
# 1. DGP: same 2^3 design, latent moments, and binary SDs as S01
# ------------------------------------------------------------
design <- as.matrix(expand.grid(
  X1 = c(-1, 1), X2 = c(-1, 1), X3 = c(-1, 1)
))
X <- cbind(F0 = 1, design)
rownames(X) <- paste0('B', 1:8)

mu <- c(F0 = .20, F1 = .20, F2 = .20, F3 = .20)
Phi <- matrix(c(
  1, .30, .20, .10,
  .30, 1, .25, .20,
  .20, .25, 1.50, .30,
  .10, .20, .30, 2
), 4, 4, byrow = TRUE) * .25

F <- MASS::mvrnorm(N, mu = mu, Sigma = Phi)
colnames(F) <- names(mu)
sd_B <- c(.80, .80, .80, .80, 1.60, 1.60, 1.60, 1.60)
tau_B <- rep(.50, 8)
tau_B[8] <- tau_B[8]+ delta
names(tau_B) <- paste0('B', 1:8)

utility <- F %*% t(X)
e_B <- sweep(matrix(rnorm(N * 8), N, 8), 2, sd_B, '*')
B <- as.data.frame(1L * (sweep(utility + e_B, 2, tau_B, '-') > 0))
names(B) <- paste0('B', 1:8)



