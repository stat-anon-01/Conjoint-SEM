# DEMONSTRATION 3: Preference comparability across response formats
# 00_DGP.R — Data-generation function only
#
# make_study2_data() generates rating and binary responses for the same
# respondents. Population part-worth means/covariances agree across formats,
# but respondent-level preferences correlate imperfectly (rho = 0.70).
# 01_ANALYSIS.R calls this function with N = 1000 and seed = 12345.
# The optional paired-comparison format is not used in this demonstration.
# ------------------------------------------------------------

make_study2_data <- function(
    N = 1000L,
    formats = c("R", "B"),
    rho = 0.70,
    seed = NULL) {

  if (!requireNamespace("MASS", quietly = TRUE))
    stop("Install package MASS.")
  if (!is.null(seed)) set.seed(seed)
  formats <- match.arg(formats, c("R", "B", "P"), several.ok = TRUE)
  if (!identical(formats, c("R", "B")) &&
      !identical(formats, c("R", "B", "P")))
    stop("Use formats=c('R','B') or c('R','B','P').")
  if (length(N) != 1L || N < 2 || N != as.integer(N)) stop("Invalid N")
  if (length(rho) != 1L || !is.finite(rho) || rho < 0 || rho >= 1)
    stop("For non-common scenarios rho must be in [0,1).")

  fn <- paste0("F", 0:3)
  mu <- setNames(rep(0.20, 4), fn)
  Phi <- matrix(c(1,.30,.20,.10,
                  .30,1,.25,.20,
                  .20,.25,1.50,.30,
                  .10,.20,.30,2), 4, 4, byrow=TRUE) * .25
  dimnames(Phi) <- list(fn, fn)
  stopifnot(min(eigen(Phi, symmetric=TRUE)$values) > 0)

  design <- as.matrix(expand.grid(X1=c(-1,1), X2=c(-1,1), X3=c(-1,1)))
  rownames(design) <- paste0("Profile", 1:8)
  X <- cbind(F0=1, design)
  pairs <- utils::combn(1:8, 2)
  pair_names <- apply(pairs, 2, function(z) paste0("P", z[1], "_", z[2]))
  pair_design <- X[pairs[1,],,drop=FALSE] - X[pairs[2,],,drop=FALSE]
  rownames(pair_design) <- pair_names

  # HETERO measurement-error settings; these are DGP values, not identification constraints.
  sigma_R <- c(.80,.80,.60,.60,1,1,1.20,1.20)
  sigma_B <- sigma_R
  pair_sd <- c(.80,1,.70,1.20,.90,1.10,.60,
               1.30,.75,1.15,.85,1.25,.95,
               .65,1.05,.90,1.20,.80,
               1.10,.70,1.30,.85,
               1,.60,1.15,.75,1.25,.95)
  names(pair_sd) <- pair_names
  tau <- .50  # Binary DGP threshold; estimation may fix thresholds to zero.

  draw_F <- function() {
    z <- MASS::mvrnorm(N, mu=rep(0,4), Sigma=Phi)
    colnames(z) <- fn
    z
  }
  # One DGP only: identical marginal latent distributions,
  # but format-specific individual preferences (rho < 1).
  # Cov(F_format_a, F_format_b) = rho * Phi for a != b.
  common <- draw_F()
  F <- lapply(formats, function(fmt) {
    z <- sqrt(rho) * common + sqrt(1 - rho) * draw_F()
    f <- sweep(z, 2, mu, `+`)
    colnames(f) <- fn
    f
  })
  names(F) <- formats

  eta <- lapply(F, function(f) f %*% t(X))
  observed <- list()
  if ("R" %in% formats) {
    y <- eta$R + sweep(matrix(stats::rnorm(N*8), N, 8), 2, sigma_R, `*`)
    observed$R <- as.data.frame(y)
    names(observed$R) <- paste0("R", 1:8)
  }
  if ("B" %in% formats) {
    y <- eta$B + sweep(matrix(stats::rnorm(N*8), N, 8), 2, sigma_B, `*`)
    observed$B <- as.data.frame(1L * (y > tau))
    names(observed$B) <- paste0("B", 1:8)
  }
  if ("P" %in% formats) {
    delta <- eta$P[,pairs[1,],drop=FALSE] - eta$P[,pairs[2,],drop=FALSE]
    y <- delta + sweep(matrix(stats::rnorm(N*28), N, 28), 2, pair_sd, `*`)
    observed$P <- as.data.frame(1L * (y > 0))
    names(observed$P) <- pair_names
  }

  # Population truth under the single Study 2 DGP.
  # If Binary thresholds are fixed to zero in the fitted model,
  # its location parameter is F0_B - tau (before scale normalization).
  mu_B_tau0 <- mu
  mu_B_tau0["F0"] <- mu["F0"] - tau
  true_means <- setNames(replicate(length(formats), mu, simplify=FALSE), formats)
  true_cov <- kronecker(
    matrix(rho, length(formats), length(formats)) +
      diag(1-rho, length(formats)), Phi)
  dn <- unlist(lapply(formats, function(fmt) paste0(fmt,"_",fn)))
  dimnames(true_cov) <- list(dn, dn)

  dat <- do.call(cbind, observed)
  # cbind() prefixes columns with list names (e.g., R.R1, B.B1).
  # Restore the variable names used by the explicit lavaan models.
  names(dat) <- unlist(lapply(observed, names), use.names = FALSE)
  rownames(dat) <- NULL
  stopifnot(!anyDuplicated(names(dat)))
  stopifnot(ncol(dat) == if (length(formats)==2) 16L else 44L,
            all(as.matrix(dat[,grep("^(B|P)",names(dat)),drop=FALSE]) %in% 0:1))
  list(data=dat, observed=observed, F=F,
       truth=list(means=true_means, covariance=true_cov,
                  Phi_base=Phi, mu_base=mu, rho=rho,
                  sigma_R=sigma_R, sigma_B=sigma_B,
                  pair_sd=pair_sd, tau=tau,
                  mu_B_tau0=mu_B_tau0),
       design=design, X=X, pairs=pairs,
       pair_names=pair_names, pair_design=pair_design,
       formats=formats, N=N)
}

# No lavaan model is generated here.
# Examples (uncomment to run):
# d_RB  <- make_study2_data(1000, c("R","B"), seed=12345)
# d_RBP <- make_study2_data(1000, c("R","B","P"), seed=123)
# dat <- d_RBP$data
