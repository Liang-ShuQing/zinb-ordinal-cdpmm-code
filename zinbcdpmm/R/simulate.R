# Data generators
# Source: 模拟研究/prog1_cdpmm_joint_vs_separate.R (algorithm bodies preserved).

generate_correlated_covariates <- function(N, p, rho = 0.5) {
  Sigma_X <- outer(1:p, 1:p, function(j, k) rho^abs(j - k))
  mvtnorm::rmvnorm(N, mean = rep(0, p), sigma = Sigma_X)
}

generate_scenario2_covariates <- function(N) {
  # Drop I(U>0) to avoid collinearity with U (improves alpha1 coverage).
  # zero: U, Z, B  -> 3 cov (+ intercept -> p_zero = 4)
  x1_zero <- runif(N, -sqrt(3), sqrt(3))
  x2_zero <- rnorm(N, 0, 1)
  x3_zero <- rbinom(N, 1, 0.5)
  X_zero_cov <- cbind(x1_zero, x2_zero, x3_zero)

  # count: U, Z, B, C  -> 4 cov (+ intercept -> p_count = 5)
  x1_count <- runif(N, -sqrt(3), sqrt(3))
  x2_count <- rnorm(N, 0, 1)
  x3_count <- rbinom(N, 1, 0.5)
  x4_count <- sample(0:2, N, replace = TRUE)
  X_count_cov <- cbind(x1_count, x2_count, x3_count, x4_count)

  # ordinal: U, Z  -> 2 cov (+ intercept -> p_ordinal = 3)
  x1_ord <- runif(N, -sqrt(3), sqrt(3))
  x2_ord <- rnorm(N, 0, 1)
  X_ordinal_cov <- cbind(x1_ord, x2_ord)

  list(
    X_ordinal_cov = X_ordinal_cov,
    X_zero_cov = X_zero_cov,
    X_count_cov = X_count_cov
  )
}

# Fixed-effect truths: Scenario 1 keeps AR(1) dims 5/6/4;
# Scenario 2 (no I(U>0)) uses dims 4/5/3 with the I(U>0) slope removed.
true_fixed_effects <- function(scenario) {
  if (as.integer(scenario)[1] == 1L) {
    list(
      alpha = c(0.3, -0.8, 0.6, -0.4, 0.5),
      beta = c(0.8, 0.6, -0.7, 0.4, -0.5, 0.3),
      gamma = c(-0.2, 0.9, -0.6, 0.4)
    )
  } else {
    list(
      alpha = c(0.3, -0.8, -0.4, 0.5),
      beta = c(0.8, 0.6, 0.4, -0.5, 0.3),
      gamma = c(-0.2, 0.9, 0.4)
    )
  }
}

# Centered finite mixture of normals for random effects (CDPMM-style truth)
# Returns list(b, Sigma, pi, mu, Omega_list)
generate_mixture_re <- function(n, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  q <- 3L
  # Uncentered means and weights (will be recentered)
  # Column order: (zero, count, ordinal) = (b1, b2, b3)
  pi_true <- c(0.50, 0.30, 0.20)
  mu_star <- rbind(
    c( 0.8,  0.6,  1.2),
    c(-1.0, -0.4, -0.9),
    c( 0.3, -0.9, -0.6)
  )
  # Permute existing Omegas from old (ord, zero, count) = (1,2,3) to (zero, count, ord)
  idx <- c(2L, 3L, 1L)
  Omega_list_old <- list(
    matrix(c(0.6, 0.25, 0.20,
             0.25, 0.7, 0.22,
             0.20, 0.22, 0.5), 3, 3, byrow = TRUE),
    matrix(c(0.9, 0.35, 0.30,
             0.35, 0.8, 0.28,
             0.30, 0.28, 0.7), 3, 3, byrow = TRUE),
    matrix(c(0.5, 0.15, 0.18,
             0.15, 0.6, 0.20,
             0.18, 0.20, 0.8), 3, 3, byrow = TRUE)
  )
  Omega_list <- lapply(Omega_list_old, function(O) O[idx, idx, drop = FALSE])
  mu_bar <- as.numeric(crossprod(pi_true, mu_star))
  mu <- sweep(mu_star, 2, mu_bar, "-")  # zero mean
  
  Sigma <- matrix(0, q, q)
  for (g in seq_along(pi_true)) {
    Sigma <- Sigma + pi_true[g] * (Omega_list[[g]] + tcrossprod(mu[g, ]))
  }
  
  L <- sample.int(length(pi_true), n, replace = TRUE, prob = pi_true)
  b <- matrix(NA_real_, n, q)
  for (i in seq_len(n)) {
    g <- L[i]
    b[i, ] <- as.numeric(mvtnorm::rmvnorm(1, mu[g, ], Omega_list[[g]]))
  }
  # numerical recenter (tiny drift from finite n)
  b <- sweep(b, 2, colMeans(b), "-")
  
  list(
    b = b,
    Sigma = Sigma,
    pi = pi_true,
    mu = mu,
    Omega_list = Omega_list,
    L = L
  )
}

# Multivariate-t random effects (mean 0).
# scale_Sigma is the t-distribution scale; population Cov = df/(df-2) * scale_Sigma (df>2).
generate_mvt_re <- function(n, df = 4, rho = 0.5, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (df <= 2) stop("df must be > 2 for finite covariance.")
  R <- matrix(c(
    1.0, rho, rho,
    rho, 1.0, rho,
    rho, rho, 1.0
  ), 3, 3, byrow = TRUE)
  # Choose scale so that Cov diagonals = 1 (same marginal variance as normal case)
  scale_Sigma <- ((df - 2) / df) * R
  Sigma <- (df / (df - 2)) * scale_Sigma  # = R, by construction
  # Scale-mixture construction (avoids rmvt edge cases):
  # b | u ~ N(0, scale_Sigma / u), u ~ Chi^2(df)/df
  Z <- mvtnorm::rmvnorm(n, sigma = scale_Sigma)
  u <- stats::rchisq(n, df = df) / df
  b <- Z / sqrt(u)
  list(
    b = b,
    Sigma = Sigma,
    scale_Sigma = scale_Sigma,
    df = df
  )
}

generate_data <- function(n = 100, nis = 10, seed = 2026,
                          random_nis = FALSE, nis_range = 1:20,
                          scenario = 1, C = 5,
                          re_dist = "normal",
                          mvt_df = 4,
                          rho = NULL,
                          zi_intercept = NULL) {
  set.seed(seed)
  if (random_nis) {
    nis <- sample(nis_range, n, replace = TRUE)
  } else if (length(nis) == 1) {
    nis <- rep(nis, n)
  }
  if (length(nis) != n) stop("nis must be a scalar or a vector of length n.")
  id <- rep(1:n, times = nis)
  N <- length(id)
  # Common RE correlation for normal/mvt truths (default 0.5)
  rho_use <- if (is.null(rho)) 0.5 else as.numeric(rho)[1]
  if (!is.finite(rho_use)) stop("rho must be finite.")
  # Keep Sigma PD when |rho| hits the boundary of the plotting grid
  rho_gen <- max(min(rho_use, 0.999), -0.999)

  if (scenario == 1) {
    X_zero_cov <- generate_correlated_covariates(N, p = 4, rho = 0.5)
    X_count_cov <- generate_correlated_covariates(N, p = 5, rho = 0.5)
    X_ordinal_cov <- generate_correlated_covariates(N, p = 3, rho = 0.5)
  } else if (scenario == 2) {
    covs <- generate_scenario2_covariates(N)
    X_ordinal_cov <- covs$X_ordinal_cov
    X_zero_cov <- covs$X_zero_cov
    X_count_cov <- covs$X_count_cov
  } else {
    stop("scenario must be 1 or 2.")
  }

  X_ordinal <- cbind(1, X_ordinal_cov)
  X_zero <- cbind(1, X_zero_cov)
  X_count <- cbind(1, X_count_cov)

  fe <- true_fixed_effects(scenario)
  alpha_true <- fe$alpha
  beta_true <- fe$beta
  gamma_true <- fe$gamma
  if (ncol(X_zero) != length(alpha_true) ||
      ncol(X_count) != length(beta_true) ||
      ncol(X_ordinal) != length(gamma_true)) {
    stop("Covariate dimension does not match true fixed-effect length.")
  }
  r_true <- 2
  # Optional override of zero-inflation intercept (smaller -> higher zero rate)
  if (!is.null(zi_intercept)) {
    zi_intercept <- as.numeric(zi_intercept)[1]
    if (!is.finite(zi_intercept)) stop("zi_intercept must be finite.")
    alpha_true[1] <- zi_intercept
  }

  if (identical(re_dist, "mixture")) {
    mix <- generate_mixture_re(n)
    b_true <- mix$b
    Sigma_true <- mix$Sigma
    re_info <- list(pi = mix$pi, mu = mix$mu, Omega_list = mix$Omega_list)
  } else if (identical(re_dist, "mvt")) {
    mt <- generate_mvt_re(n, df = mvt_df, rho = rho_gen)
    b_true <- mt$b
    Sigma_true <- mt$Sigma
    re_info <- list(df = mt$df, scale_Sigma = mt$scale_Sigma)
  } else if (identical(re_dist, "normal")) {
    Sigma_true <- matrix(c(
      1.0, rho_gen, rho_gen,
      rho_gen, 1.0, rho_gen,
      rho_gen, rho_gen, 1.0
    ), 3, 3, byrow = TRUE)
    b_true <- mvtnorm::rmvnorm(n, sigma = Sigma_true)
    re_info <- NULL
  } else {
    stop("re_dist must be 'normal', 'mixture', or 'mvt'.")
  }

  b1_true <- b_true[, 1]  # zero
  b2_true <- b_true[, 2]  # count
  b3_true <- b_true[, 3]  # ordinal
  colnames(b_true) <- c("zero", "count", "ord")

  # y1 = ZINB count: alpha + b1 (zero), beta + b2 (count among at-risk)
  eta_zero <- as.numeric(X_zero %*% alpha_true + rep(b1_true, times = nis))
  pi_at_risk <- inv_logit(eta_zero)
  u_true <- rbinom(N, 1, pi_at_risk)

  y1 <- rep(0L, N)
  pos_idx <- which(u_true == 1)
  if (length(pos_idx) > 0) {
    eta_count <- as.numeric(X_count[pos_idx, ] %*% beta_true + rep(b2_true, times = nis)[pos_idx])
    phi <- inv_logit(eta_count)
    mu <- r_true * phi / (1 - phi)
    mu <- pmax(mu, 1e-10)
    y1[pos_idx] <- rnbinom(length(pos_idx), size = r_true, mu = mu)
  }

  # y2 = ordinal: gamma + b3
  eta_ord <- as.numeric(X_ordinal %*% gamma_true + rep(b3_true, times = nis))
  l_true <- eta_ord + rlogis(N, 0, 1)
  delta_true <- seq(0, 3, length.out = C - 1)
  y2 <- as.integer(cut(
    l_true,
    breaks = c(-Inf, delta_true, Inf),
    labels = FALSE,
    include.lowest = TRUE
  ))

  list(
    id = id, N = N, n = n, nis = nis, C = C,
    X_ordinal = X_ordinal, X_zero = X_zero, X_count = X_count,
    y1 = y1, y2 = y2,
    b_true = b_true,
    alpha_true = alpha_true,
    beta_true = beta_true,
    gamma_true = gamma_true,
    delta_true = delta_true,
    r_true = r_true,
    Sigma_true = Sigma_true,
    re_dist = re_dist,
    re_info = re_info
  )
}

