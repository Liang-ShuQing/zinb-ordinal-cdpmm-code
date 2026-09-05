# Program 2 (standalone): CDPMM joint vs Gaussian joint (same data per replicate)
# Does NOT fit separate/univariate models.
# Reports DIC / WAIC / LOOIC win-rate %: P(CDPMM joint < Gaussian joint)
# CDPMM atoms: strict normal--inverse-Wishart (NIW), full Omega_g
############################################################
# 0. Local study settings (edit these for local runs)
#    Results -> sibling folder 模拟研究结果/
#    Optional env overrides (SCENARIO, N, N_SIM, CHAIN, ...) used on HPC.
############################################################
set.seed(2026)
n <- 100
rho1 <- 0.5
C <- 5
scenario <- 2
# Random-effects truth: "normal" | "mixture" | "mvt"
re_dist <- "normal"
mvt_df <- 4
# This program always fits BOTH CDPMM and Gaussian joint models on the same data
# (no separate/univariate models).
# Paper-scale defaults (slow locally — reduce n_sim/chain for smoke tests)
n_sim <- 128
chain <- 5000
burn <- 2000
thin <- 5
nis_fixed <- 10
random_nis <- TRUE
nis_range <- 6:9
delta_min <- -10
delta_max <- 10
# CDPMM truncation / hyperparameters (strict full-covariance NIW)
G_mix <- 8
tau_a1 <- 2
tau_a2 <- 4
zeta0_sd2 <- 10
kappa0 <- 1
iw_nu0 <- 6
iw_S0_scale <- 1

.env_int <- function(key, default) {
  v <- suppressWarnings(as.integer(Sys.getenv(key, unset = "")))
  if (length(v) == 1L && is.finite(v)) v else default
}
.env_chr <- function(key, default) {
  v <- Sys.getenv(key, unset = "")
  if (nzchar(v)) v else default
}
scenario <- .env_int("SCENARIO", scenario)
n <- .env_int("N", n)
n_sim <- .env_int("N_SIM", n_sim)
chain <- .env_int("CHAIN", chain)
burn <- .env_int("BURN", burn)
thin <- .env_int("THIN", thin)
re_dist <- .env_chr("RE_DIST", re_dist)

# Keep caller OUT_DIR if already set (e.g. plot_re_density.R); else CLI / sibling 结果
if (!exists("OUT_DIR") || is.null(OUT_DIR) || !nzchar(as.character(OUT_DIR)[1])) {
  args_cli <- commandArgs(trailingOnly = TRUE)
  OUT_DIR <- if (length(args_cli) >= 1) {
    args_cli[[1]]
  } else {
    env_out <- Sys.getenv("OUT_DIR", unset = "")
    if (nzchar(env_out)) {
      env_out
    } else {
      code_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)
      parent <- dirname(code_dir)
      base <- basename(code_dir)
      if (grepl("[\u4e00-\u9fff]", base)) {
        file.path(parent, paste0(base, "结果"))
      } else {
        file.path(parent, paste0(base, "_results"))
      }
    }
  }
}

library(BayesLogit)
library(mvtnorm)
library(MCMCpack)
library(truncnorm)
library(parallel)
library(doParallel)
library(foreach)
library(loo)
library(coda)

############################################################
# 1. Utility functions
############################################################

dnbinom_zero <- function(r, mu) {
  prob <- (r / (r + mu))^r
  prob[is.na(prob)] <- 1
  prob <- pmin(pmax(prob, 1e-10), 1 - 1e-10)
  prob
}

log_sum_exp <- function(logx) {
  m <- max(logx)
  m + log(sum(exp(logx - m)))
}

inv_logit <- function(x) plogis(x)

safe_cor <- function(x, y) {
  out <- suppressWarnings(cor(x, y))
  if (is.na(out) || !is.finite(out)) out <- NA_real_
  out
}

ordinal_loglik_vec <- function(y, eta, delta) {
  cuts <- c(-Inf, delta, Inf)
  ll <- numeric(length(y))
  for (i in seq_along(y)) {
    k <- y[i]
    prob <- plogis(cuts[k + 1] - eta[i]) - plogis(cuts[k] - eta[i])
    prob <- pmin(pmax(prob, 1e-10), 1 - 1e-10)
    ll[i] <- log(prob)
  }
  ll
}

update_latent_ordinal <- function(y, eta, delta, omega1, C) {
  N <- length(y)
  l <- numeric(N)
  for (i in seq_len(N)) {
    k <- y[i]
    sd_i <- 1 / sqrt(max(omega1[i], 1e-10))
    if (k == 1) {
      l[i] <- rtruncnorm(1, a = -Inf, b = delta[1], mean = eta[i], sd = sd_i)
    } else if (k == C) {
      l[i] <- rtruncnorm(1, a = delta[C - 1], b = Inf, mean = eta[i], sd = sd_i)
    } else {
      l[i] <- rtruncnorm(1, a = delta[k - 1], b = delta[k], mean = eta[i], sd = sd_i)
    }
  }
  l[!is.finite(l)] <- eta[!is.finite(l)]
  l
}

enforce_delta_order <- function(delta, fix_first = TRUE) {
  delta_new <- delta
  if (fix_first) {
    delta_new[1] <- 0
    start_k <- 2L
  } else {
    start_k <- 1L
  }
  for (k in start_k:length(delta_new)) {
    if (k == start_k && !fix_first) {
      next
    }
    min_val <- if (k == 1) -Inf else delta_new[k - 1] + 1e-4
    if (!is.finite(min_val) || delta_new[k] <= min_val) {
      delta_new[k] <- min_val + 0.5
    }
  }
  delta_new
}

update_thresholds <- function(y, l, delta, C, delta_min = -10, delta_max = 10,
                              fix_first = TRUE) {
  delta_new <- delta
  if (fix_first) {
    delta_new[1] <- 0
  }
  start_k <- if (fix_first) 2L else 1L
  for (k in start_k:(C - 1)) {
    idx_c <- which(y == k)
    idx_c1 <- which(y == k + 1)

    lower_candidates <- c(delta_min)
    if (k > 1) lower_candidates <- c(lower_candidates, delta_new[k - 1])
    if (length(idx_c) > 0) lower_candidates <- c(lower_candidates, max(l[idx_c]))
    lower_bound <- max(lower_candidates)

    upper_candidates <- c(delta_max)
    if (k < C - 1) upper_candidates <- c(upper_candidates, delta_new[k + 1])
    if (length(idx_c1) > 0) upper_candidates <- c(upper_candidates, min(l[idx_c1]))
    upper_bound <- min(upper_candidates)

    if (lower_bound < upper_bound) {
      delta_new[k] <- runif(1, lower_bound, upper_bound)
    } else {
      delta_new[k] <- (lower_bound + upper_bound) / 2
    }
  }
  enforce_delta_order(delta_new, fix_first = fix_first)
}

compute_dic_from_loglik <- function(loglik_mat, theta_bar_loglik_fun) {
  dev_s <- -2 * rowSums(loglik_mat)
  D_bar <- mean(dev_s, na.rm = TRUE)
  D_hat <- theta_bar_loglik_fun()
  pD <- D_bar - D_hat
  DIC <- D_hat + 2 * pD
  list(DIC = DIC, D_bar = D_bar, D_hat = D_hat, pD = pD)
}

compute_waic <- function(loglik_mat) {
  waic_res <- loo::waic(loglik_mat)
  waic_val <- waic_res$estimates["waic", "Estimate"]
  list(value = waic_val, obj = waic_res)
}

compute_looic <- function(loglik_mat) {
  loo_res <- loo::loo(loglik_mat)
  looic_val <- loo_res$estimates["looic", "Estimate"]
  list(value = looic_val, obj = loo_res)
}

calculate_stats <- function(estimates, true_value, ci_lower, ci_upper) {
  n_params <- length(true_value)
  bias <- rep(NA_real_, n_params)
  rmse <- rep(NA_real_, n_params)
  ci_lower_mean <- rep(NA_real_, n_params)
  ci_upper_mean <- rep(NA_real_, n_params)
  cp <- rep(NA_real_, n_params)

  for (j in seq_len(n_params)) {
    bias[j] <- mean(estimates[, j] - true_value[j], na.rm = TRUE)
    rmse[j] <- sqrt(mean((estimates[, j] - true_value[j])^2, na.rm = TRUE))
    ci_lower_mean[j] <- mean(ci_lower[, j], na.rm = TRUE)
    ci_upper_mean[j] <- mean(ci_upper[, j], na.rm = TRUE)
    cp[j] <- mean(ci_lower[, j] <= true_value[j] & true_value[j] <= ci_upper[, j], na.rm = TRUE)
  }

  list(
    bias = bias,
    rmse = rmse,
    ci_lower_mean = ci_lower_mean,
    ci_upper_mean = ci_upper_mean,
    cp = cp
  )
}

calculate_stats_varying_true <- function(estimates, true_mat, ci_lower, ci_upper) {
  n_params <- ncol(estimates)
  bias <- rep(NA_real_, n_params)
  rmse <- rep(NA_real_, n_params)
  ci_lower_mean <- rep(NA_real_, n_params)
  ci_upper_mean <- rep(NA_real_, n_params)
  cp <- rep(NA_real_, n_params)

  for (j in seq_len(n_params)) {
    bias[j] <- mean(estimates[, j] - true_mat[, j], na.rm = TRUE)
    rmse[j] <- sqrt(mean((estimates[, j] - true_mat[, j])^2, na.rm = TRUE))
    ci_lower_mean[j] <- mean(ci_lower[, j], na.rm = TRUE)
    ci_upper_mean[j] <- mean(ci_upper[, j], na.rm = TRUE)
    cp[j] <- mean(ci_lower[, j] <= true_mat[, j] & true_mat[, j] <= ci_upper[, j], na.rm = TRUE)
  }

  list(
    bias = bias,
    rmse = rmse,
    ci_lower_mean = ci_lower_mean,
    ci_upper_mean = ci_upper_mean,
    cp = cp
  )
}

safe_col_var <- function(x) {
  apply(x, 2, function(z) stats::var(z, na.rm = TRUE))
}

# Stick-breaking weights from nu in (0,1)^{G-1}, with nu_G forced to 1
stick_break_weights <- function(nu) {
  G <- length(nu)
  pi <- numeric(G)
  rem <- 1
  for (g in seq_len(G - 1L)) {
    pi[g] <- nu[g] * rem
    rem <- rem * (1 - nu[g])
  }
  pi[G] <- rem
  s <- sum(pi)
  if (!is.finite(s) || s <= 0) rep(1 / G, G) else pi / s
}

# Implied zero-mean mixture covariance with full component covariances
# Omega_list: length-G list of q x q matrices; mu: G x q (already centered)
cdpmm_implied_Sigma <- function(pi, mu, Omega_list) {
  q <- ncol(mu)
  Sigma <- matrix(0, q, q)
  for (g in seq_along(pi)) {
    if (pi[g] < 1e-16) next
    Sigma <- Sigma + pi[g] * (Omega_list[[g]] + tcrossprod(mu[g, ]))
  }
  eig <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
  if (any(!is.finite(eig)) || min(eig) < 1e-8) {
    Sigma <- Sigma + diag(1e-6, q)
  }
  Sigma
}

safe_solve <- function(A, eps = 1e-8) {
  A <- (A + t(A)) / 2
  eig <- eigen(A, symmetric = TRUE, only.values = TRUE)$values
  if (any(!is.finite(eig)) || min(eig) < eps) {
    A <- A + diag(eps - min(c(eig[is.finite(eig)], 0), na.rm = TRUE) + eps, nrow(A))
  }
  solve(A)
}

safe_riwish <- function(nu, S) {
  S <- (S + t(S)) / 2
  eig <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  if (any(!is.finite(eig)) || min(eig) <= 1e-10) {
    S <- S + diag(1e-4, nrow(S))
  }
  MCMCpack::riwish(max(nu, nrow(S) + 1), S)
}

# ---- Univariate CDPMM = 1D NIW (NIG) for separate models ----
# Prior df matched to joint 3D IW diagonal marginal:
#   Omega ~ IW_3(nu0, s0 I)  =>  Omega_ii ~ IW_1(nu0 - 3 + 1, s0)
# so 1D uses nu0_1d = iw_nu0 - 2 (with iw_nu0=6 => 4; E[omega]=s0/(4-2)=0.5
# same as E[Omega_ii]=s0/(6-3-1)=0.5). Do NOT reuse raw iw_nu0 in 1D.
cdpmm_implied_var_1d <- function(pi, mu, omega) {
  v <- sum(pi * (omega + mu^2))
  max(v, 1e-8)
}

init_cdpmm_1d <- function(n, G = G_mix) {
  G <- max(2L, as.integer(G))
  q_joint <- 3L
  nu0 <- max(iw_nu0 - q_joint + 1L, 3L)  # edge-match to joint IW
  S0 <- iw_S0_scale
  tau <- 1
  nu <- c(rbeta(G - 1L, 1, tau), 1)
  pi_w <- stick_break_weights(nu)
  zeta <- rnorm(1, 0, sqrt(zeta0_sd2))
  omega <- vapply(seq_len(G), function(g) {
    max(as.numeric(safe_riwish(nu0, matrix(S0, 1, 1))), 1e-8)
  }, numeric(1))
  # mu*|omega ~ N(zeta, omega/kappa0)
  mu_star <- rnorm(G, zeta, sqrt(omega / kappa0))
  mu_bar <- sum(pi_w * mu_star)
  mu <- mu_star - mu_bar
  L <- sample.int(G, n, replace = TRUE, prob = pi_w)
  list(
    G = G, n = n, nu0 = nu0, S0 = S0,
    tau = tau, nu = nu, pi_w = pi_w,
    zeta = zeta,
    mu_star = mu_star, mu = mu, omega = omega, L = L
  )
}

# Strict 1D NIW updates (Tang stick-breaking/centering retained)
update_cdpmm_1d <- function(st, b) {
  G <- st$G
  n <- st$n
  log_pi <- log(pmax(st$pi_w, 1e-300))
  log_dens <- matrix(log_pi, n, G, byrow = TRUE)
  for (g in seq_len(G)) {
    og <- max(st$omega[g], 1e-8)
    resid <- b - st$mu[g]
    log_dens[, g] <- log_dens[, g] -
      0.5 * log(2 * base::pi * og) - 0.5 * (resid^2) / og
  }
  log_dens <- log_dens - apply(log_dens, 1, max)
  dens <- exp(log_dens)
  dens <- dens / pmax(rowSums(dens), 1e-300)
  st$L <- vapply(seq_len(n), function(j) sample.int(G, 1, prob = dens[j, ]), integer(1))

  n_g <- tabulate(st$L, nbins = G)
  for (g in seq_len(G - 1L)) {
    a_g <- 1 + n_g[g]
    b_g <- st$tau + sum(n_g[(g + 1):G])
    st$nu[g] <- rbeta(1, a_g, max(b_g, 1e-8))
  }
  st$nu[G] <- 1
  st$pi_w <- stick_break_weights(st$nu)

  sum_log <- sum(log(pmax(1 - st$nu[seq_len(G - 1L)], 1e-12)))
  st$tau <- rgamma(1, shape = tau_a1 + (G - 1),
                   rate = max(tau_a2 - sum_log, 1e-8))

  # atoms: NIW then zeta | {mu*, omega}
  for (g in seq_len(G)) {
    idx_g <- which(st$L == g)
    ng <- length(idx_g)
    if (ng == 0) {
      st$omega[g] <- max(as.numeric(safe_riwish(st$nu0, matrix(st$S0, 1, 1))), 1e-8)
      st$mu_star[g] <- rnorm(1, st$zeta, sqrt(st$omega[g] / kappa0))
    } else {
      b_g <- b[idx_g]
      bbar <- mean(b_g)
      scat <- sum((b_g - bbar)^2)
      kn <- kappa0 + ng
      df_post <- st$nu0 + ng
      S_post <- st$S0 + scat + (kappa0 * ng / kn) * (bbar - st$zeta)^2
      st$omega[g] <- max(as.numeric(safe_riwish(df_post, matrix(S_post, 1, 1))), 1e-8)
      m_n <- (kappa0 * st$zeta + ng * bbar) / kn
      st$mu_star[g] <- rnorm(1, m_n, sqrt(st$omega[g] / kn))
    }
  }

  # zeta | mu*, omega  (prior N(0, zeta0_sd2); mu*_g|zeta,omega ~ N(zeta, omega/kappa0))
  prec <- 1 / zeta0_sd2
  num <- 0
  for (g in seq_len(G)) {
    w <- kappa0 / max(st$omega[g], 1e-8)
    prec <- prec + w
    num <- num + w * st$mu_star[g]
  }
  st$zeta <- rnorm(1, num / prec, sqrt(1 / prec))

  mu_bar <- sum(st$pi_w * st$mu_star)
  st$mu <- st$mu_star - mu_bar

  # pure mixture variance (Tang-style reporting)
  st$sigma2 <- cdpmm_implied_var_1d(st$pi_w, st$mu, st$omega)
  st$nclust <- sum(n_g > 0)
  st
}

############################################################
############################################################
# 2. Data generator
############################################################

generate_correlated_covariates <- function(N, p, rho = 0.5) {
  Sigma_X <- outer(1:p, 1:p, function(j, k) rho^abs(j - k))
  rmvnorm(N, mean = rep(0, p), sigma = Sigma_X)
}

generate_scenario2_covariates <- function(N) {
  # Drop I(U>0): zero U,Z,B; count U,Z,B,C; ordinal U,Z
  x1_zero <- runif(N, -sqrt(3), sqrt(3))
  x2_zero <- rnorm(N, 0, 1)
  x3_zero <- rbinom(N, 1, 0.5)
  X_zero_cov <- cbind(x1_zero, x2_zero, x3_zero)

  x1_count <- runif(N, -sqrt(3), sqrt(3))
  x2_count <- rnorm(N, 0, 1)
  x3_count <- rbinom(N, 1, 0.5)
  x4_count <- sample(0:2, N, replace = TRUE)
  X_count_cov <- cbind(x1_count, x2_count, x3_count, x4_count)

  x1_ord <- runif(N, -sqrt(3), sqrt(3))
  x2_ord <- rnorm(N, 0, 1)
  X_ordinal_cov <- cbind(x1_ord, x2_ord)

  list(
    X_ordinal_cov = X_ordinal_cov,
    X_zero_cov = X_zero_cov,
    X_count_cov = X_count_cov
  )
}

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
    b[i, ] <- as.numeric(rmvnorm(1, mu[g, ], Omega_list[[g]]))
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
                          mvt_df = 4) {
  set.seed(seed)
  if (random_nis) {
    nis <- sample(nis_range, n, replace = TRUE)
  } else if (length(nis) == 1) {
    nis <- rep(nis, n)
  }
  if (length(nis) != n) stop("nis must be a scalar or a vector of length n.")
  id <- rep(1:n, times = nis)
  N <- length(id)

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

  if (identical(re_dist, "mixture")) {
    mix <- generate_mixture_re(n)
    b_true <- mix$b
    Sigma_true <- mix$Sigma
    re_info <- list(pi = mix$pi, mu = mix$mu, Omega_list = mix$Omega_list)
  } else if (identical(re_dist, "mvt")) {
    mt <- generate_mvt_re(n, df = mvt_df, rho = rho1)
    b_true <- mt$b
    Sigma_true <- mt$Sigma
    re_info <- list(df = mt$df, scale_Sigma = mt$scale_Sigma)
  } else if (identical(re_dist, "normal")) {
    Sigma_true <- matrix(c(
      1.0, rho1, rho1,
      rho1, 1.0, rho1,
      rho1, rho1, 1.0
    ), 3, 3, byrow = TRUE)
    b_true <- rmvnorm(n, sigma = Sigma_true)
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

############################################################
# 3. Joint model
#    re_prior = "cdpmm": centered DPMM (Tang et al., 2014) + full NIW atoms
#    re_prior = "gaussian": b_i ~ N(0, Sigma), Sigma ~ IW(nu0, S0)
############################################################
fit_joint_model <- function(dat, chain = 5000, burn = 2000, thin = 5,
                            delta_min = -10, delta_max = 10,
                            G = G_mix, init = NULL,
                            re_prior = "cdpmm") {
  re_prior <- tolower(as.character(re_prior)[1])
  if (!re_prior %in% c("cdpmm", "gaussian")) {
    stop("fit_joint_model: re_prior must be 'cdpmm' or 'gaussian'")
  }
  use_cdpmm <- identical(re_prior, "cdpmm")

  id <- dat$id
  N <- dat$N
  n <- dat$n
  nis <- dat$nis
  C <- dat$C

  X_ordinal <- dat$X_ordinal
  X_zero <- dat$X_zero
  X_count <- dat$X_count
  y1 <- dat$y1
  y2 <- dat$y2

  p_ordinal <- ncol(X_ordinal)
  p_zero <- ncol(X_zero)
  p_count <- ncol(X_count)
  q_re <- 3L

  alpha0 <- rep(0, p_zero)
  beta0 <- rep(0, p_count)
  gamma0 <- rep(0, p_ordinal)

  T0a <- diag(0.001, p_zero)
  T0b <- diag(0.001, p_count)
  T0g <- diag(0.001, p_ordinal)

  shape_r <- 0.01
  rate_r <- 0.01

  alpha <- rep(0, p_zero)
  beta <- rep(0, p_count)
  gamma_ord <- rep(0, p_ordinal)
  delta <- seq(0, 3, length.out = C - 1)

  b <- matrix(rnorm(n * q_re), n, q_re)
  r <- 1.0
  l <- rep(0, N)
  omega_ord <- rep(1, N)

  S0_iw <- diag(iw_S0_scale, q_re)
  nu0_iw <- max(iw_nu0, q_re + 2)

  # CDPMM state (only used when use_cdpmm)
  G <- max(2L, as.integer(G))
  zeta0 <- rep(0, q_re)
  Psi0_inv <- diag(1 / zeta0_sd2, q_re)
  tau <- 1
  nu <- c(rbeta(G - 1L, 1, tau), 1)
  pi_w <- stick_break_weights(nu)
  zeta <- as.numeric(rmvnorm(1, zeta0, diag(zeta0_sd2, q_re)))
  mu_star <- matrix(0, G, q_re)
  Omega_list <- vector("list", G)
  Omega_inv_list <- vector("list", G)
  mu <- matrix(0, G, q_re)
  L <- rep(1L, n)
  n_g <- rep(0L, G)

  # Homogeneous Gaussian RE covariance (only used when !use_cdpmm)
  Sigma <- diag(1, q_re)
  Sigma_inv <- diag(1, q_re)

  if (use_cdpmm) {
    for (g in seq_len(G)) {
      Omega_list[[g]] <- safe_riwish(nu0_iw, S0_iw)
      Omega_inv_list[[g]] <- safe_solve(Omega_list[[g]])
      mu_star[g, ] <- as.numeric(rmvnorm(1, zeta, Omega_list[[g]] / kappa0))
    }
    mu_bar <- as.numeric(crossprod(pi_w, mu_star))
    mu <- sweep(mu_star, 2, mu_bar, "-")
    L <- sample.int(G, n, replace = TRUE, prob = pi_w)
  } else {
    Sigma <- safe_riwish(nu0_iw, S0_iw)
    Sigma_inv <- safe_solve(Sigma)
  }

  if (!is.null(init)) {
    if (!is.null(init$alpha)) {
      stopifnot(length(init$alpha) == p_zero)
      alpha <- as.numeric(init$alpha)
    }
    if (!is.null(init$beta)) {
      stopifnot(length(init$beta) == p_count)
      beta <- as.numeric(init$beta)
    }
    if (!is.null(init$gamma)) {
      stopifnot(length(init$gamma) == p_ordinal)
      gamma_ord <- as.numeric(init$gamma)
    }
    if (!is.null(init$delta)) {
      stopifnot(length(init$delta) == C - 1)
      delta <- as.numeric(init$delta)
    }
    if (!is.null(init$r)) {
      r <- max(as.numeric(init$r)[1], 1e-4)
    }
    if (use_cdpmm && !is.null(init$tau)) {
      tau <- max(as.numeric(init$tau)[1], 1e-4)
      nu <- c(rbeta(G - 1L, 1, tau), 1)
      pi_w <- stick_break_weights(nu)
    }

    if (use_cdpmm) {
      atoms_set <- FALSE
      if (!is.null(init$Omega_list)) {
        stopifnot(is.list(init$Omega_list), length(init$Omega_list) == G)
        for (g in seq_len(G)) {
          Og <- as.matrix(init$Omega_list[[g]])
          stopifnot(nrow(Og) == q_re, ncol(Og) == q_re)
          eig0 <- eigen(Og, symmetric = TRUE, only.values = TRUE)$values
          if (any(!is.finite(eig0)) || min(eig0) < 1e-8) Og <- Og + diag(1e-4, q_re)
          Omega_list[[g]] <- Og
          Omega_inv_list[[g]] <- safe_solve(Og)
        }
        atoms_set <- TRUE
      } else {
        Sig0 <- init$Omega
        if (is.null(Sig0)) Sig0 <- init$Sigma
        if (!is.null(Sig0)) {
          Sig0 <- as.matrix(Sig0)
          stopifnot(nrow(Sig0) == q_re, ncol(Sig0) == q_re)
          eig0 <- eigen(Sig0, symmetric = TRUE, only.values = TRUE)$values
          if (any(!is.finite(eig0)) || min(eig0) < 1e-8) {
            Sig0 <- Sig0 + diag(1e-4, q_re)
          }
          for (g in seq_len(G)) {
            Omega_list[[g]] <- Sig0
            Omega_inv_list[[g]] <- safe_solve(Sig0)
          }
          atoms_set <- TRUE
        }
      }

      if (!is.null(init$mu_star)) {
        ms <- init$mu_star
        if (is.vector(ms) && length(ms) == q_re) {
          mu_star <- matrix(ms, G, q_re, byrow = TRUE)
        } else {
          stopifnot(is.matrix(ms), nrow(ms) == G, ncol(ms) == q_re)
          mu_star <- ms
        }
        atoms_set <- TRUE
      } else if (atoms_set && is.null(init$Omega_list) &&
                 (!is.null(init$Omega) || !is.null(init$Sigma))) {
        mu_star <- matrix(0, G, q_re)
      }

      if (atoms_set) {
        zeta <- rep(0, q_re)
        mu_bar <- as.numeric(crossprod(pi_w, mu_star))
        mu <- sweep(mu_star, 2, mu_bar, "-")
      }

      if (!is.null(init$L)) {
        stopifnot(length(init$L) == n)
        L <- as.integer(init$L)
      } else if (atoms_set) {
        L <- sample.int(G, n, replace = TRUE, prob = pi_w)
      }

      if (!is.null(init$b)) {
        if (is.vector(init$b) && length(init$b) == q_re) {
          b <- matrix(init$b, n, q_re, byrow = TRUE)
        } else {
          stopifnot(is.matrix(init$b), nrow(init$b) == n, ncol(init$b) == q_re)
          b <- init$b
        }
      } else if (atoms_set) {
        for (i in seq_len(n)) {
          g <- L[i]
          b[i, ] <- as.numeric(rmvnorm(1, mu[g, ], Omega_list[[g]]))
        }
      }
    } else {
      Sig0 <- init$Sigma
      if (is.null(Sig0)) Sig0 <- init$Omega
      if (!is.null(Sig0)) {
        Sig0 <- as.matrix(Sig0)
        stopifnot(nrow(Sig0) == q_re, ncol(Sig0) == q_re)
        eig0 <- eigen(Sig0, symmetric = TRUE, only.values = TRUE)$values
        if (any(!is.finite(eig0)) || min(eig0) < 1e-8) {
          Sig0 <- Sig0 + diag(1e-4, q_re)
        }
        Sigma <- Sig0
        Sigma_inv <- safe_solve(Sigma)
      }
      if (!is.null(init$b)) {
        if (is.vector(init$b) && length(init$b) == q_re) {
          b <- matrix(init$b, n, q_re, byrow = TRUE)
        } else {
          stopifnot(is.matrix(init$b), nrow(init$b) == n, ncol(init$b) == q_re)
          b <- init$b
        }
      }
    }
  }

  b1 <- b[, 1]
  b2 <- b[, 2]
  b3 <- b[, 3]
  u_est <- as.numeric(y1 > 0)

  id_index <- split(seq_len(N), id)

  save_every <- floor((chain - burn) / thin)
  Alpha_store <- matrix(NA, save_every, p_zero)
  Beta_store <- matrix(NA, save_every, p_count)
  Gamma_store <- matrix(NA, save_every, p_ordinal)
  Delta_store <- matrix(NA, save_every, C - 1)
  Sigma_store <- matrix(NA, save_every, 9)
  R_store <- rep(NA, save_every)
  Rho_store <- matrix(NA, save_every, 3)
  Tau_store <- rep(NA, save_every)
  Nclust_store <- rep(NA, save_every)

  loglik_y1 <- matrix(NA, save_every, N)
  loglik_y2 <- matrix(NA, save_every, N)

  zero_rep <- rep(NA, save_every)
  cor_rep <- rep(NA, save_every)
  B_sum <- matrix(0, n, q_re)
  n_b_saved <- 0L

  for (iter in seq_len(chain)) {
    # ordinal (y2): gamma + b3
    eta_ord <- as.numeric(X_ordinal %*% gamma_ord + rep(b3, times = nis))
    l <- update_latent_ordinal(y2, eta_ord, delta, omega_ord, C)
    psi_pg <- pmin(pmax(l - eta_ord, -50), 50)
    omega_ord <- rpg(N, 2, psi_pg)

    # zero-inflation (y1): alpha + b1
    eta_zero <- as.numeric(X_zero %*% alpha + rep(b1, times = nis))
    eta_zero <- pmin(pmax(eta_zero, -10), 10)
    pi_at_risk <- inv_logit(eta_zero)
    pi_at_risk <- pmin(pmax(pi_at_risk, 1e-10), 1 - 1e-10)

    # count (y1 among at-risk): beta + b2
    eta_count_all <- as.numeric(X_count %*% beta + rep(b2, times = nis))
    eta_count_all <- pmin(pmax(eta_count_all, -10), 10)
    phi <- inv_logit(eta_count_all)
    phi <- pmin(pmax(phi, 1e-10), 1 - 1e-10)

    mu_nb <- r * phi / (1 - phi)
    mu_nb <- pmax(mu_nb, 1e-10)
    q_nb0 <- dnbinom_zero(r, mu_nb)

    for (i in seq_len(N)) {
      if (y1[i] == 0) {
        log_p1 <- log(pi_at_risk[i]) + log(q_nb0[i])
        log_p0 <- log(1 - pi_at_risk[i])
        denom <- log_sum_exp(c(log_p1, log_p0))
        theta <- exp(log_p1 - denom)
        theta <- pmin(pmax(theta, 1e-10), 1 - 1e-10)
        u_est[i] <- rbinom(1, 1, theta)
      } else {
        u_est[i] <- 1
      }
    }

    pos_idx <- which(u_est == 1)
    n_pos <- length(pos_idx)

    omega_zero <- rpg(N, 1, eta_zero)
    z_alpha <- (u_est - 0.5) / pmax(omega_zero, 1e-10)

    if (n_pos > 0) {
      y1_pos <- y1[pos_idx]
      X_pos <- X_count[pos_idx, , drop = FALSE]
      b2_pos <- rep(b2, times = nis)[pos_idx]

      eta_count_pos <- as.numeric(X_pos %*% beta + b2_pos)
      eta_count_pos <- pmin(pmax(eta_count_pos, -10), 10)
      w_count <- rpg(n_pos, y1_pos + r, eta_count_pos)
      z_beta <- (y1_pos - r) / (2 * pmax(w_count, 1e-10))
    }

    V_alpha <- solve(crossprod(sqrt(omega_zero) * X_zero) + T0a)
    m_alpha <- V_alpha %*% (T0a %*% alpha0 +
                              crossprod(X_zero, omega_zero * (z_alpha - rep(b1, times = nis))))
    alpha <- as.numeric(rmvnorm(1, m_alpha, V_alpha))

    if (n_pos > 0) {
      V_beta <- solve(crossprod(sqrt(w_count) * X_pos) + T0b)
      m_beta <- V_beta %*% (T0b %*% beta0 +
                              crossprod(X_pos, w_count * (z_beta - b2_pos)))
      beta <- as.numeric(rmvnorm(1, m_beta, V_beta))
    }

    V_gamma <- solve(crossprod(sqrt(omega_ord) * X_ordinal) + T0g)
    m_gamma <- V_gamma %*% (T0g %*% gamma0 +
                              crossprod(X_ordinal, omega_ord * (l - rep(b3, times = nis))))
    gamma_ord <- as.numeric(rmvnorm(1, m_gamma, V_gamma))

    delta <- update_thresholds(y2, l, delta, C, delta_min, delta_max, fix_first = TRUE)

    if (n_pos > 0) {
      k_crp <- rep(0, n_pos)
      for (j in seq_len(n_pos)) {
        if (y1_pos[j] > 0) {
          probs <- r / (r + 0:(y1_pos[j] - 1))
          probs <- pmin(pmax(probs, 1e-10), 1 - 1e-10)
          k_crp[j] <- sum(rbinom(y1_pos[j], 1, probs))
        }
      }

      eta_count_current <- as.numeric(X_pos %*% beta + b2_pos)
      psi_current <- inv_logit(eta_count_current)
      psi_current <- pmin(pmax(psi_current, 1e-10), 1 - 1e-10)

      log_term <- sum(log(1 - psi_current))
      if (!is.finite(log_term)) log_term <- 0

      r_shape_post <- shape_r + sum(k_crp)
      r_rate_post <- rate_r - log_term
      r_rate_post <- max(r_rate_post, 1e-10)
      r <- rgamma(1, shape = r_shape_post, rate = r_rate_post)
    }

    # ---- random effects b = (b1 zero, b2 count, b3 ordinal) ----
    for (j in seq_len(n)) {
      idx <- id_index[[as.character(j)]]
      if (is.null(idx)) idx <- id_index[[j]]
      if (use_cdpmm) {
        g <- L[j]
        prior_prec <- Omega_inv_list[[g]]
        prior_mean <- mu[g, ]
      } else {
        prior_prec <- Sigma_inv
        prior_mean <- rep(0, q_re)
      }

      Xo_j <- X_ordinal[idx, , drop = FALSE]
      Xz_j <- X_zero[idx, , drop = FALSE]

      o_ord <- omega_ord[idx]
      d_ord <- sum(o_ord * (l[idx] - as.numeric(Xo_j %*% gamma_ord)))
      w_ord_sum <- sum(o_ord)

      o_zero <- omega_zero[idx]
      d_zero <- sum(o_zero * (z_alpha[idx] - as.numeric(Xz_j %*% alpha)))
      w_zero_sum <- sum(o_zero)

      w_count_sum <- 0
      d_count <- 0

      idx_pos <- idx[idx %in% pos_idx]
      if (length(idx_pos) > 0 && n_pos > 0) {
        Xc_j <- X_count[idx_pos, , drop = FALSE]
        pos_in_pos <- match(idx_pos, pos_idx)
        wg_j <- w_count[pos_in_pos]
        zg_j <- z_beta[pos_in_pos]
        w_count_sum <- sum(wg_j)
        d_count <- sum(wg_j * (zg_j - as.numeric(Xc_j %*% beta)))
      }

      W_j <- diag(c(w_zero_sum, w_count_sum, w_ord_sum), q_re)
      d_j <- c(d_zero, d_count, d_ord)

      post_prec <- prior_prec + W_j
      post_var <- safe_solve(post_prec)
      post_mean <- as.numeric(post_var %*% (prior_prec %*% prior_mean + d_j))
      b[j, ] <- as.numeric(rmvnorm(1, post_mean, post_var))
    }

    b1 <- b[, 1]
    b2 <- b[, 2]
    b3 <- b[, 3]

    if (use_cdpmm) {
      # ---- CDPMM blocked Gibbs updates (full-covariance NIW) ----
      log_pi <- log(pmax(pi_w, 1e-300))
      log_dens <- matrix(log_pi, n, G, byrow = TRUE)
      for (g in seq_len(G)) {
        R <- tryCatch(chol(Omega_list[[g]]), error = function(e) NULL)
        if (is.null(R)) {
          Omega_list[[g]] <- Omega_list[[g]] + diag(1e-4, q_re)
          R <- chol(Omega_list[[g]])
          Omega_inv_list[[g]] <- safe_solve(Omega_list[[g]])
        }
        resid <- sweep(b, 2, mu[g, ], "-")
        z_std <- t(backsolve(R, t(resid), transpose = TRUE))
        log_dens[, g] <- log_dens[, g] -
          sum(log(diag(R))) - 0.5 * q_re * log(2 * base::pi) -
          0.5 * rowSums(z_std^2)
      }
      log_dens <- log_dens - apply(log_dens, 1, max)
      dens <- exp(log_dens)
      dens <- dens / pmax(rowSums(dens), 1e-300)
      L <- vapply(seq_len(n), function(j) sample.int(G, 1, prob = dens[j, ]), integer(1))

      n_g <- tabulate(L, nbins = G)
      for (g in seq_len(G - 1L)) {
        a_g <- 1 + n_g[g]
        b_g <- tau + sum(n_g[(g + 1):G])
        nu[g] <- rbeta(1, a_g, max(b_g, 1e-8))
      }
      nu[G] <- 1
      pi_w <- stick_break_weights(nu)

      sum_log <- sum(log(pmax(1 - nu[seq_len(G - 1L)], 1e-12)))
      tau_rate <- tau_a2 - sum_log
      tau <- rgamma(1, shape = tau_a1 + (G - 1), rate = max(tau_rate, 1e-8))

      for (g in seq_len(G)) {
        idx_g <- which(L == g)
        ng <- length(idx_g)
        if (ng == 0) {
          Omega_list[[g]] <- safe_riwish(nu0_iw, S0_iw)
          Omega_inv_list[[g]] <- safe_solve(Omega_list[[g]])
          mu_star[g, ] <- as.numeric(rmvnorm(1, zeta, Omega_list[[g]] / kappa0))
        } else {
          b_g <- b[idx_g, , drop = FALSE]
          bbar <- colMeans(b_g)
          b_c <- sweep(b_g, 2, bbar, "-")
          Scat <- crossprod(b_c)
          kn <- kappa0 + ng
          df_post <- nu0_iw + ng
          S_post <- S0_iw + Scat + (kappa0 * ng / kn) * tcrossprod(bbar - zeta)
          Omega_list[[g]] <- safe_riwish(df_post, S_post)
          Omega_inv_list[[g]] <- safe_solve(Omega_list[[g]])
          m_n <- (kappa0 * zeta + ng * bbar) / kn
          mu_star[g, ] <- as.numeric(rmvnorm(1, m_n, Omega_list[[g]] / kn))
        }
      }

      Prec_zeta <- Psi0_inv
      num_zeta <- Psi0_inv %*% zeta0
      for (g in seq_len(G)) {
        Wg <- kappa0 * Omega_inv_list[[g]]
        Prec_zeta <- Prec_zeta + Wg
        num_zeta <- num_zeta + Wg %*% mu_star[g, ]
      }
      V_zeta <- safe_solve(Prec_zeta)
      zeta <- as.numeric(rmvnorm(1, V_zeta %*% num_zeta, V_zeta))

      mu_bar <- as.numeric(crossprod(pi_w, mu_star))
      mu <- sweep(mu_star, 2, mu_bar, "-")

      Sigma <- cdpmm_implied_Sigma(pi_w, mu, Omega_list)
      eig <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
      if (any(!is.finite(eig)) || min(eig) < 1e-8) {
        Sigma <- Sigma + diag(1e-4, q_re)
      }
    } else {
      # ---- Homogeneous Gaussian: Sigma | b ~ IW(nu0 + n, S0 + b'b) ----
      Scat <- crossprod(b)
      Sigma <- safe_riwish(nu0_iw + n, S0_iw + Scat)
      eig <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
      if (any(!is.finite(eig)) || min(eig) < 1e-8) {
        Sigma <- Sigma + diag(1e-4, q_re)
      }
      Sigma_inv <- safe_solve(Sigma)
    }

    if (iter > burn && ((iter - burn) %% thin == 0)) {
      s <- (iter - burn) / thin
      Alpha_store[s, ] <- alpha
      Beta_store[s, ] <- beta
      Gamma_store[s, ] <- gamma_ord
      Delta_store[s, ] <- delta
      Sigma_store[s, ] <- c(Sigma)
      R_store[s] <- r
      # Rho order: zero-count, zero-ordinal, count-ordinal (Sigma[1,2], [1,3], [2,3])
      Rho_store[s, ] <- c(
        Sigma[1, 2] / sqrt(Sigma[1, 1] * Sigma[2, 2]),
        Sigma[1, 3] / sqrt(Sigma[1, 1] * Sigma[3, 3]),
        Sigma[2, 3] / sqrt(Sigma[2, 2] * Sigma[3, 3])
      )
      Tau_store[s] <- if (use_cdpmm) tau else NA_real_
      Nclust_store[s] <- if (use_cdpmm) sum(n_g > 0) else 1L
      B_sum <- B_sum + b
      n_b_saved <- n_b_saved + 1L

      eta_ord_s <- as.numeric(X_ordinal %*% gamma_ord + rep(b3, times = nis))
      loglik_y2[s, ] <- ordinal_loglik_vec(y2, eta_ord_s, delta)

      eta_zero_s <- as.numeric(X_zero %*% alpha + rep(b1, times = nis))
      pi_s <- inv_logit(eta_zero_s)

      eta_count_s <- as.numeric(X_count %*% beta + rep(b2, times = nis))
      phi_s <- inv_logit(eta_count_s)
      mu_s <- r * phi_s / (1 - phi_s)
      mu_s <- pmax(mu_s, 1e-10)

      loglik_y1[s, ] <- ifelse(
        y1 == 0,
        log((1 - pi_s) + pi_s * dnbinom(y1, size = r, mu = mu_s)),
        log(pi_s) + dnbinom(y1, size = r, mu = mu_s, log = TRUE)
      )

      zero_rep[s] <- mean(y1 == 0)
      cor_rep[s] <- safe_cor(y2, as.numeric(y1 > 0))
    }
  }

  b_est <- if (n_b_saved > 0L) B_sum / n_b_saved else b

  alpha_mean <- colMeans(Alpha_store, na.rm = TRUE)
  beta_mean <- colMeans(Beta_store, na.rm = TRUE)
  gamma_mean <- colMeans(Gamma_store, na.rm = TRUE)
  delta_mean <- colMeans(Delta_store, na.rm = TRUE)
  r_mean <- mean(R_store, na.rm = TRUE)

  theta_bar_loglik_fun <- function() {
    eta_ord <- as.numeric(X_ordinal %*% gamma_mean)
    ll_ord <- ordinal_loglik_vec(y2, eta_ord, delta_mean)

    eta_zero <- as.numeric(X_zero %*% alpha_mean)
    pi_s <- inv_logit(eta_zero)

    eta_count <- as.numeric(X_count %*% beta_mean)
    phi_s <- inv_logit(eta_count)
    mu_s <- r_mean * phi_s / (1 - phi_s)
    mu_s <- pmax(mu_s, 1e-10)

    ll_zinb <- ifelse(
      y1 == 0,
      log((1 - pi_s) + pi_s * dnbinom(y1, size = r_mean, mu = mu_s)),
      log(pi_s) + dnbinom(y1, size = r_mean, mu = mu_s, log = TRUE)
    )

    -2 * sum(ll_ord + ll_zinb)
  }

  loglik_total <- loglik_y1 + loglik_y2
  dic_res <- compute_dic_from_loglik(loglik_total, theta_bar_loglik_fun)
  waic_res <- compute_waic(loglik_total)
  loo_res <- compute_looic(loglik_total)

  list(
    alpha_samples = Alpha_store,
    beta_samples = Beta_store,
    gamma_samples = Gamma_store,
    delta_samples = Delta_store,
    r_samples = R_store,
    Sigma_samples = Sigma_store,
    tau_samples = Tau_store,
    nclust_samples = Nclust_store,
    b_est = b_est,
    alpha_est = alpha_mean,
    beta_est = beta_mean,
    gamma_est = gamma_mean,
    delta_est = delta_mean,
    r_est = r_mean,
    Sigma_est = matrix(colMeans(Sigma_store, na.rm = TRUE), 3, 3),
    Rho_est = colMeans(Rho_store, na.rm = TRUE),
    tau_est = mean(Tau_store, na.rm = TRUE),
    nclust_est = mean(Nclust_store, na.rm = TRUE),
    alpha_ci_lower = apply(Alpha_store, 2, function(x) quantile(x, 0.025, na.rm = TRUE)),
    alpha_ci_upper = apply(Alpha_store, 2, function(x) quantile(x, 0.975, na.rm = TRUE)),
    beta_ci_lower = apply(Beta_store, 2, function(x) quantile(x, 0.025, na.rm = TRUE)),
    beta_ci_upper = apply(Beta_store, 2, function(x) quantile(x, 0.975, na.rm = TRUE)),
    gamma_ci_lower = apply(Gamma_store, 2, function(x) quantile(x, 0.025, na.rm = TRUE)),
    gamma_ci_upper = apply(Gamma_store, 2, function(x) quantile(x, 0.975, na.rm = TRUE)),
    delta_ci_lower = apply(Delta_store, 2, function(x) quantile(x, 0.025, na.rm = TRUE)),
    delta_ci_upper = apply(Delta_store, 2, function(x) quantile(x, 0.975, na.rm = TRUE)),
    r_ci_lower = quantile(R_store, 0.025, na.rm = TRUE),
    r_ci_upper = quantile(R_store, 0.975, na.rm = TRUE),
    Sigma_ci_lower = matrix(apply(Sigma_store, 2, function(x) quantile(x, 0.025, na.rm = TRUE)), 3, 3),
    Sigma_ci_upper = matrix(apply(Sigma_store, 2, function(x) quantile(x, 0.975, na.rm = TRUE)), 3, 3),
    Rho_ci_lower = apply(Rho_store, 2, function(x) quantile(x, 0.025, na.rm = TRUE)),
    Rho_ci_upper = apply(Rho_store, 2, function(x) quantile(x, 0.975, na.rm = TRUE)),
    dic = dic_res$DIC,
    waic = waic_res$value,
    looic = loo_res$value,
    loglik = loglik_total,
    ppc = list(
      obs_zero_prop = mean(y1 == 0),
      rep_zero_prop_mean = mean(zero_rep, na.rm = TRUE),
      rep_zero_prop_ci = quantile(zero_rep, c(0.025, 0.975), na.rm = TRUE),
      obs_cor_y1_y2pos = safe_cor(y2, as.numeric(y1 > 0)),
      obs_cor_y1_y2 = safe_cor(y2, y1),
      rep_cor_y1_y2pos_mean = mean(cor_rep, na.rm = TRUE),
      rep_cor_y1_y2pos_ci = quantile(cor_rep, c(0.025, 0.975), na.rm = TRUE)
    )
  )
}

############################################################
# 6. Repeated simulation settings
############################################################
# Set SKIP_MAIN_SIM <- TRUE before source() to load helpers only
# (e.g. plot_re_density.R).
if (!exists("SKIP_MAIN_SIM") || !isTRUE(SKIP_MAIN_SIM)) {

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

n_cores <- suppressWarnings(as.integer(Sys.getenv(
  c("N_CORES", "SLURM_NTASKS_PER_NODE", "SLURM_CPUS_PER_TASK")
)))
n_cores <- n_cores[is.finite(n_cores) & n_cores >= 1L]
n_cores <- if (length(n_cores)) n_cores[[1]] else max(1L, parallel::detectCores())
n_cores <- min(n_cores, n_sim)
registerDoParallel(cores = n_cores)

cat(paste0(
  "\n=== Prog2: CDPMM joint vs Gaussian joint | 重复模拟 ", n_sim, " 次",
  "，scenario=", scenario,
  "，re_dist=", re_dist,
  "，chain=", chain, "/", burn,
  " thin=", thin, "，cores=", n_cores,
  "\nOUT_DIR=", OUT_DIR, " ===\n"
))

############################################################
# 7. Repeated simulation
############################################################
sim_results <- foreach(i = 1:n_sim,
                       .packages = c("BayesLogit", "mvtnorm", "MCMCpack", "truncnorm", "loo", "coda")) %dopar% {
                         tryCatch({
                           dat <- generate_data(
                             n = n, nis = nis_fixed, seed = 2026 + i,
                             random_nis = random_nis, nis_range = nis_range,
                             scenario = scenario, C = C,
                             re_dist = re_dist, mvt_df = mvt_df
                           )
                           joint_cdpmm <- fit_joint_model(
                             dat, chain = chain, burn = burn, thin = thin,
                             delta_min = delta_min, delta_max = delta_max,
                             re_prior = "cdpmm"
                           )
                           joint_gauss <- fit_joint_model(
                             dat, chain = chain, burn = burn, thin = thin,
                             delta_min = delta_min, delta_max = delta_max,
                             re_prior = "gaussian"
                           )
                           list(dat = dat, joint_cdpmm = joint_cdpmm,
                                joint_gauss = joint_gauss)
                         }, error = function(e) {
                           cat(paste0("第 ", i, " 次模拟失败: ", e$message, "\n"))
                           NULL
                         })
                       }

stopImplicitCluster()

sim_results <- sim_results[!sapply(sim_results, is.null)]
n_success <- length(sim_results)

cat(paste0("\n=== 成功完成 ", n_success, "/", n_sim, " 次模拟 ===\n"))
if (n_success == 0) {
  stop("All simulation replications failed. Please check the model fitting functions or reduce the number of cores.")
}

############################################################
# 8. True values
############################################################
fe_true <- true_fixed_effects(scenario)
alpha_true <- fe_true$alpha
beta_true <- fe_true$beta
gamma_true <- fe_true$gamma
delta_true <- seq(0, 3, length.out = C - 1)
r_true <- 2

# Default (normal / mvt): unit diagonals and common correlation rho1
Sigma_true <- matrix(c(
  1.0, rho1, rho1,
  rho1, 1.0, rho1,
  rho1, rho1, 1.0
), 3, 3, byrow = TRUE)

if (isTRUE(re_dist == "mixture")) {
  Sigma_true <- generate_mixture_re(n = 2)$Sigma
  cat("\n=== Mixture-truth implied Sigma ===\n")
  print(round(Sigma_true, 4))
} else if (isTRUE(re_dist == "mvt")) {
  cat("\n=== Multivariate-t truth (df=", mvt_df,
      ") population Sigma (unit vars, corr=", rho1, ") ===\n", sep = "")
  print(round(Sigma_true, 4))
} else {
  cat("\n=== Normal-truth Sigma ===\n")
  print(round(Sigma_true, 4))
}

rho_true <- c(
  Sigma_true[1, 2] / sqrt(Sigma_true[1, 1] * Sigma_true[2, 2]),
  Sigma_true[1, 3] / sqrt(Sigma_true[1, 1] * Sigma_true[3, 3]),
  Sigma_true[2, 3] / sqrt(Sigma_true[2, 2] * Sigma_true[3, 3])
)


make_joint_param_summary <- function(res_list, key) {
  n_ok <- length(res_list)
  p_a <- length(alpha_true); p_b <- length(beta_true); p_g <- length(gamma_true)
  alpha_est <- matrix(NA, n_ok, p_a)
  alpha_lo <- matrix(NA, n_ok, p_a); alpha_hi <- matrix(NA, n_ok, p_a)
  beta_est <- matrix(NA, n_ok, p_b)
  beta_lo <- matrix(NA, n_ok, p_b); beta_hi <- matrix(NA, n_ok, p_b)
  gamma_est <- matrix(NA, n_ok, p_g)
  gamma_lo <- matrix(NA, n_ok, p_g); gamma_hi <- matrix(NA, n_ok, p_g)
  r_est <- rep(NA, n_ok); r_lo <- rep(NA, n_ok); r_hi <- rep(NA, n_ok)
  Sigma_est <- array(NA, dim = c(n_ok, 3, 3))
  Sigma_lo <- array(NA, dim = c(n_ok, 3, 3)); Sigma_hi <- array(NA, dim = c(n_ok, 3, 3))
  rho_est <- matrix(NA, n_ok, 3)
  rho_lo <- matrix(NA, n_ok, 3); rho_hi <- matrix(NA, n_ok, 3)
  crit <- matrix(NA, n_ok, 3)
  for (i in seq_len(n_ok)) {
    fit <- res_list[[i]][[key]]
    alpha_est[i, ] <- fit$alpha_est; alpha_lo[i, ] <- fit$alpha_ci_lower; alpha_hi[i, ] <- fit$alpha_ci_upper
    beta_est[i, ] <- fit$beta_est; beta_lo[i, ] <- fit$beta_ci_lower; beta_hi[i, ] <- fit$beta_ci_upper
    gamma_est[i, ] <- fit$gamma_est; gamma_lo[i, ] <- fit$gamma_ci_lower; gamma_hi[i, ] <- fit$gamma_ci_upper
    r_est[i] <- fit$r_est; r_lo[i] <- fit$r_ci_lower; r_hi[i] <- fit$r_ci_upper
    Sigma_est[i, , ] <- fit$Sigma_est
    Sigma_lo[i, , ] <- fit$Sigma_ci_lower; Sigma_hi[i, , ] <- fit$Sigma_ci_upper
    rho_est[i, ] <- fit$Rho_est; rho_lo[i, ] <- fit$Rho_ci_lower; rho_hi[i, ] <- fit$Rho_ci_upper
    crit[i, ] <- c(fit$dic, fit$waic, fit$looic)
  }
  a_s <- calculate_stats(alpha_est, alpha_true, alpha_lo, alpha_hi)
  b_s <- calculate_stats(beta_est, beta_true, beta_lo, beta_hi)
  g_s <- calculate_stats(gamma_est, gamma_true, gamma_lo, gamma_hi)
  r_s <- calculate_stats(matrix(r_est, ncol = 1), matrix(r_true, ncol = 1),
                         matrix(r_lo, ncol = 1), matrix(r_hi, ncol = 1))
  Sigma_bias <- apply(Sigma_est, c(2, 3), mean, na.rm = TRUE) - Sigma_true
  Sigma_rmse <- sqrt(apply((Sigma_est - array(rep(Sigma_true, each = n_ok),
                                              dim = c(n_ok, 3, 3)))^2, c(2, 3), mean, na.rm = TRUE))
  Sigma_cp <- matrix(NA, 3, 3)
  for (a in 1:3) for (b in 1:3) {
    Sigma_cp[a, b] <- mean(Sigma_lo[, a, b] <= Sigma_true[a, b] &
                             Sigma_true[a, b] <= Sigma_hi[, a, b], na.rm = TRUE)
  }
  rho_s <- calculate_stats(rho_est, rho_true, rho_lo, rho_hi)
  summary_df <- rbind(
    data.frame(Parameter = paste0("alpha", seq_len(p_a)), True_Value = alpha_true,
               Bias = round(a_s$bias, 4), RMSE = round(a_s$rmse, 4), CP = round(a_s$cp, 3)),
    data.frame(Parameter = paste0("beta", seq_len(p_b)), True_Value = beta_true,
               Bias = round(b_s$bias, 4), RMSE = round(b_s$rmse, 4), CP = round(b_s$cp, 3)),
    data.frame(Parameter = paste0("gamma", seq_len(p_g)), True_Value = gamma_true,
               Bias = round(g_s$bias, 4), RMSE = round(g_s$rmse, 4), CP = round(g_s$cp, 3)),
    data.frame(Parameter = "r", True_Value = r_true,
               Bias = round(r_s$bias, 4), RMSE = round(r_s$rmse, 4), CP = round(r_s$cp, 3)),
    data.frame(Parameter = c("Sigma11", "Sigma22", "Sigma33"), True_Value = diag(Sigma_true),
               Bias = round(diag(Sigma_bias), 4), RMSE = round(diag(Sigma_rmse), 4),
               CP = round(diag(Sigma_cp), 3)),
    data.frame(Parameter = paste0("rho", 1:3), True_Value = rho_true,
               Bias = round(rho_s$bias, 4), RMSE = round(rho_s$rmse, 4), CP = round(rho_s$cp, 3))
  )
  list(summary = summary_df, crit = crit,
       avg_cp = mean(c(a_s$cp, b_s$cp, g_s$cp, r_s$cp, diag(Sigma_cp), rho_s$cp), na.rm = TRUE),
       avg_rmse = mean(c(a_s$rmse, b_s$rmse, g_s$rmse, r_s$rmse, diag(Sigma_rmse), rho_s$rmse), na.rm = TRUE))
}

sum_c <- make_joint_param_summary(sim_results, "joint_cdpmm")
sum_g <- make_joint_param_summary(sim_results, "joint_gauss")

cat("\n================ Prog2: CDPMM joint vs Gaussian joint ================\n")
cat("CDPMM avg CP=", round(sum_c$avg_cp, 3), " RMSE=", round(sum_c$avg_rmse, 3), "\n")
cat("Gauss avg CP=", round(sum_g$avg_cp, 3), " RMSE=", round(sum_g$avg_rmse, 3), "\n")
cat("\n--- CDPMM joint ---\n"); print(sum_c$summary, row.names = FALSE)
cat("\n--- Gaussian joint ---\n"); print(sum_g$summary, row.names = FALSE)

crit_c <- sum_c$crit
crit_g <- sum_g$crit
win_dic <- mean(crit_c[, 1] < crit_g[, 1], na.rm = TRUE) * 100
win_waic <- mean(crit_c[, 2] < crit_g[, 2], na.rm = TRUE) * 100
win_looic <- mean(crit_c[, 3] < crit_g[, 3], na.rm = TRUE) * 100
criteria_winrate <- data.frame(
  Comparison = "CDPMM_joint_better_than_Gaussian_joint",
  DIC_WinPct = round(win_dic, 1),
  WAIC_WinPct = round(win_waic, 1),
  LOOIC_WinPct = round(win_looic, 1),
  DIC_Mean_CDPMM = mean(crit_c[, 1], na.rm = TRUE),
  DIC_Mean_Gauss = mean(crit_g[, 1], na.rm = TRUE),
  WAIC_Mean_CDPMM = mean(crit_c[, 2], na.rm = TRUE),
  WAIC_Mean_Gauss = mean(crit_g[, 2], na.rm = TRUE),
  LOOIC_Mean_CDPMM = mean(crit_c[, 3], na.rm = TRUE),
  LOOIC_Mean_Gauss = mean(crit_g[, 3], na.rm = TRUE)
)
cat("\n================ 准则胜率 (CDPMM joint 优于 Gaussian joint 的百分比) ================\n")
print(criteria_winrate, row.names = FALSE)

# Posterior variance ratio Gauss / CDPMM (>1 => CDPMM more precise)
safe_col_var <- function(M) {
  if (is.null(dim(M))) return(stats::var(M, na.rm = TRUE))
  apply(M, 2, stats::var, na.rm = TRUE)
}
vr_list <- lapply(seq_len(n_success), function(i) {
  jc <- sim_results[[i]]$joint_cdpmm
  jg <- sim_results[[i]]$joint_gauss
  c(
    safe_col_var(jg$alpha_samples) / safe_col_var(jc$alpha_samples),
    safe_col_var(jg$beta_samples) / safe_col_var(jc$beta_samples),
    safe_col_var(jg$gamma_samples) / safe_col_var(jc$gamma_samples),
    stats::var(jg$r_samples, na.rm = TRUE) / stats::var(jc$r_samples, na.rm = TRUE),
    stats::var(jg$Sigma_samples[, 1], na.rm = TRUE) / stats::var(jc$Sigma_samples[, 1], na.rm = TRUE),
    stats::var(jg$Sigma_samples[, 5], na.rm = TRUE) / stats::var(jc$Sigma_samples[, 5], na.rm = TRUE),
    stats::var(jg$Sigma_samples[, 9], na.rm = TRUE) / stats::var(jc$Sigma_samples[, 9], na.rm = TRUE)
  )
})
vr_mat <- do.call(rbind, vr_list)
vr_names <- c(paste0("alpha", seq_along(alpha_true)),
              paste0("beta", seq_along(beta_true)),
              paste0("gamma", seq_along(gamma_true)),
              "r", "sigma2_b1", "sigma2_b2", "sigma2_b3")
variance_ratio_gauss_over_cdpmm <- data.frame(
  Parameter = vr_names,
  Variance_Ratio_Mean = colMeans(vr_mat, na.rm = TRUE),
  Variance_Ratio_TrimMean = apply(vr_mat, 2, mean, trim = 0.1, na.rm = TRUE)
)
cat("\n================ 后验方差比 (Gaussian / CDPMM; >1 表示 CDPMM 更精) ================\n")
print(variance_ratio_gauss_over_cdpmm, row.names = FALSE)

tag <- paste0("scen", scenario, "_", re_dist, "_cdpmm_vs_gauss_n", n, "_chain", chain)
write.csv(sum_c$summary, file.path(OUT_DIR, paste0("joint_param_summary_cdpmm_", tag, ".csv")), row.names = FALSE)
write.csv(sum_g$summary, file.path(OUT_DIR, paste0("joint_param_summary_gauss_", tag, ".csv")), row.names = FALSE)
write.csv(criteria_winrate, file.path(OUT_DIR, paste0("criteria_winrate_", tag, ".csv")), row.names = FALSE)
write.csv(variance_ratio_gauss_over_cdpmm,
          file.path(OUT_DIR, paste0("variance_ratio_gauss_over_cdpmm_", tag, ".csv")), row.names = FALSE)
saveRDS(list(
  settings = list(scenario = scenario, re_dist = re_dist,
                  n = n, n_sim = n_sim, chain = chain, burn = burn, thin = thin,
                  n_success = n_success),
  joint_cdpmm = sum_c$summary, joint_gauss = sum_g$summary,
  criteria_winrate = criteria_winrate,
  variance_ratio_gauss_over_cdpmm = variance_ratio_gauss_over_cdpmm
), file.path(OUT_DIR, paste0("sim_summary_", tag, ".rds")))
cat("Wrote Prog2 summaries to: ", OUT_DIR, "\n", sep = "")


} # end SKIP_MAIN_SIM guard
