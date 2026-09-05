# Program 1 (standalone): CDPMM joint model vs separate models
# Ordinal + ZINB; joint RE prior = multivariate CDPMM; separate = univariate CDPMM
# Reports DIC / WAIC / LOOIC win-rate %: P(joint < sum of separate)
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
# Tag only (this program is CDPMM-only; no Gaussian joint prior)
re_prior <- "cdpmm"
# Paper-scale defaults (slow locally — reduce n_sim/chain for smoke tests)
n_sim <- 200
chain <- 10000
burn <- 5000
thin <- 5
nis_fixed <- 10
random_nis <- TRUE
# Advisor: avoid Uniform{1..20} (too unbalanced); use 6--9 with min 6
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
G_mix <- .env_int("G_MIX", G_mix)
if (!is.finite(G_mix) || G_mix < 2L) G_mix <- 8L

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
  # Common RE correlation for normal/mvt truths (defaults to global rho1)
  rho_use <- if (is.null(rho)) rho1 else as.numeric(rho)[1]
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
# 3. Joint model (CDPMM random effects only)
#    centered DPMM (Tang et al., 2014) + full NIW atoms
############################################################
fit_joint_model <- function(dat, chain = 5000, burn = 2000, thin = 5,
                            delta_min = -10, delta_max = 10,
                            G = G_mix, init = NULL,
                            alpha0 = NULL, beta0 = NULL, gamma0 = NULL,
                            prior_var = 1000,
                            r_update = c("crt", "mh"),
                            mh_r_sd = 0.20) {
  r_update <- match.arg(r_update)
  mh_r_sd <- as.numeric(mh_r_sd)[1]
  if (!is.finite(mh_r_sd) || mh_r_sd <= 0) stop("mh_r_sd must be positive")
  mh_accept <- 0L
  mh_prop <- 0L
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

  # FE prior: N(mean, prior_var * I)  <=>  precision T0 = (1/prior_var) I
  # Default prior_var=1000 matches paper T0=0.001 I and zero means.
  if (is.null(alpha0)) alpha0 <- rep(0, p_zero)
  if (is.null(beta0)) beta0 <- rep(0, p_count)
  if (is.null(gamma0)) gamma0 <- rep(0, p_ordinal)
  alpha0 <- as.numeric(alpha0)
  beta0 <- as.numeric(beta0)
  gamma0 <- as.numeric(gamma0)
  if (length(alpha0) != p_zero) stop("alpha0 length mismatch")
  if (length(beta0) != p_count) stop("beta0 length mismatch")
  if (length(gamma0) != p_ordinal) stop("gamma0 length mismatch")
  prior_var <- as.numeric(prior_var)[1]
  if (!is.finite(prior_var) || prior_var <= 0) stop("prior_var must be positive")
  prec0 <- 1 / prior_var
  T0a <- diag(prec0, p_zero)
  T0b <- diag(prec0, p_count)
  T0g <- diag(prec0, p_ordinal)

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

  # CDPMM state
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
  Sigma <- diag(1, q_re)

  for (g in seq_len(G)) {
    Omega_list[[g]] <- safe_riwish(nu0_iw, S0_iw)
    Omega_inv_list[[g]] <- safe_solve(Omega_list[[g]])
    mu_star[g, ] <- as.numeric(rmvnorm(1, zeta, Omega_list[[g]] / kappa0))
  }
  mu_bar <- as.numeric(crossprod(pi_w, mu_star))
  mu <- sweep(mu_star, 2, mu_bar, "-")
  L <- sample.int(G, n, replace = TRUE, prob = pi_w)

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
    if (!is.null(init$tau)) {
      tau <- max(as.numeric(init$tau)[1], 1e-4)
      nu <- c(rbeta(G - 1L, 1, tau), 1)
      pi_w <- stick_break_weights(nu)
    }

    {
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
      eta_count_current <- as.numeric(X_pos %*% beta + b2_pos)
      psi_current <- inv_logit(eta_count_current)
      psi_current <- pmin(pmax(psi_current, 1e-10), 1 - 1e-10)

      if (identical(r_update, "crt")) {
        k_crp <- rep(0, n_pos)
        for (j in seq_len(n_pos)) {
          if (y1_pos[j] > 0) {
            probs <- r / (r + 0:(y1_pos[j] - 1))
            probs <- pmin(pmax(probs, 1e-10), 1 - 1e-10)
            k_crp[j] <- sum(rbinom(y1_pos[j], 1, probs))
          }
        }

        log_term <- sum(log(1 - psi_current))
        if (!is.finite(log_term)) log_term <- 0

        r_shape_post <- shape_r + sum(k_crp)
        r_rate_post <- rate_r - log_term
        r_rate_post <- max(r_rate_post, 1e-10)
        r <- rgamma(1, shape = r_shape_post, rate = r_rate_post)
      } else {
        # Log-random-walk MH on r using the NB full conditional (no CRT).
        # Target density on r>0; proposal is Gaussian RW on log(r).
        log_post_r <- function(rr) {
          if (!is.finite(rr) || rr <= 0) return(-Inf)
          lp <- (shape_r - 1) * log(rr) - rate_r * rr
          lp <- lp + sum(lgamma(y1_pos + rr) - lgamma(rr))
          lp <- lp + rr * sum(log(1 - psi_current))
          if (!is.finite(lp)) return(-Inf)
          lp
        }
        log_r_cur <- log(max(r, 1e-8))
        log_r_prop <- log_r_cur + mh_r_sd * rnorm(1)
        r_prop <- exp(log_r_prop)
        mh_prop <- mh_prop + 1L
        log_alpha <- log_post_r(r_prop) - log_post_r(r) + log(r_prop) - log(r)
        if (is.finite(log_alpha) && log(runif(1)) < log_alpha) {
          r <- r_prop
          mh_accept <- mh_accept + 1L
        }
      }
    }

    # ---- random effects b = (b1 zero, b2 count, b3 ordinal) ----
    for (j in seq_len(n)) {
      idx <- id_index[[as.character(j)]]
      if (is.null(idx)) idx <- id_index[[j]]
      g <- L[j]
      prior_prec <- Omega_inv_list[[g]]
      prior_mean <- mu[g, ]

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

    {
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
      Tau_store[s] <- tau
      Nclust_store[s] <- sum(n_g > 0)
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
    r_update = r_update,
    mh_accept_rate = if (mh_prop > 0L) mh_accept / mh_prop else NA_real_,
    mh_n_prop = mh_prop,
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
# 4. Ordinal-only model (univariate CDPMM; aligned with joint)
############################################################
fit_ordinal_model <- function(dat, chain = 5000, burn = 2000, thin = 5,
                              delta_min = -10, delta_max = 10,
                              G = G_mix) {
  id <- dat$id
  N <- dat$N
  n <- dat$n
  nis <- dat$nis
  C <- dat$C

  X <- dat$X_ordinal
  y <- dat$y2
  p <- ncol(X)

  gamma0 <- rep(0, p)
  T0 <- diag(0.001, p)

  gamma_ord <- rep(0, p)
  delta <- seq(0, 3, length.out = C - 1)
  b <- rnorm(n)
  cdp <- init_cdpmm_1d(n, G = G)
  id_index <- split(seq_len(N), id)
  l <- rep(0, N)
  omega <- rep(1, N)

  save_every <- floor((chain - burn) / thin)
  Gamma_store <- matrix(NA, save_every, p)
  Delta_store <- matrix(NA, save_every, C - 1)
  Sigma_store <- rep(NA, save_every)
  loglik <- matrix(NA, save_every, N)

  for (iter in seq_len(chain)) {
    eta <- as.numeric(X %*% gamma_ord + rep(b, times = nis))
    l <- update_latent_ordinal(y, eta, delta, omega, C)
    psi_pg <- pmin(pmax(l - eta, -50), 50)
    omega <- rpg(N, 2, psi_pg)

    V_g <- solve(crossprod(sqrt(omega) * X) + T0)
    m_g <- V_g %*% (T0 %*% gamma0 + crossprod(X, omega * (l - rep(b, times = nis))))
    gamma_ord <- as.numeric(rmvnorm(1, m_g, V_g))

    delta <- update_thresholds(y, l, delta, C, delta_min, delta_max, fix_first = TRUE)

    for (j in seq_len(n)) {
      idx <- id_index[[as.character(j)]]
      if (is.null(idx)) idx <- id_index[[j]]
      wj <- omega[idx]
      lj <- l[idx]
      Xj <- X[idx, , drop = FALSE]
      g <- cdp$L[j]
      prior_prec <- 1 / max(cdp$omega[g], 1e-8)
      prior_mean <- cdp$mu[g]
      post_var <- 1 / (prior_prec + sum(wj))
      post_mean <- post_var * (prior_prec * prior_mean +
                                 sum(wj * (lj - as.numeric(Xj %*% gamma_ord))))
      b[j] <- rnorm(1, post_mean, sqrt(post_var))
    }

    cdp <- update_cdpmm_1d(cdp, b)

    if (iter > burn && ((iter - burn) %% thin == 0)) {
      s <- (iter - burn) / thin
      Gamma_store[s, ] <- gamma_ord
      Delta_store[s, ] <- delta
      Sigma_store[s] <- cdp$sigma2

      eta_s <- as.numeric(X %*% gamma_ord + rep(b, times = nis))
      loglik[s, ] <- ordinal_loglik_vec(y, eta_s, delta)
    }
  }

  gamma_mean <- colMeans(Gamma_store, na.rm = TRUE)
  delta_mean <- colMeans(Delta_store, na.rm = TRUE)

  theta_bar_loglik_fun <- function() {
    eta <- as.numeric(X %*% gamma_mean)
    -2 * sum(ordinal_loglik_vec(y, eta, delta_mean))
  }

  dic_res <- compute_dic_from_loglik(loglik, theta_bar_loglik_fun)
  waic_res <- compute_waic(loglik)
  loo_res <- compute_looic(loglik)

  list(
    gamma_samples = Gamma_store,
    delta_samples = Delta_store,
    sigma2_samples = Sigma_store,
    sigma2_est = mean(Sigma_store, na.rm = TRUE),
    sigma2_ci_lower = quantile(Sigma_store, 0.025, na.rm = TRUE),
    sigma2_ci_upper = quantile(Sigma_store, 0.975, na.rm = TRUE),
    gamma_est = gamma_mean,
    delta_est = delta_mean,
    gamma_ci_lower = apply(Gamma_store, 2, function(x) quantile(x, 0.025, na.rm = TRUE)),
    gamma_ci_upper = apply(Gamma_store, 2, function(x) quantile(x, 0.975, na.rm = TRUE)),
    delta_ci_lower = apply(Delta_store, 2, function(x) quantile(x, 0.025, na.rm = TRUE)),
    delta_ci_upper = apply(Delta_store, 2, function(x) quantile(x, 0.975, na.rm = TRUE)),
    dic = dic_res$DIC,
    waic = waic_res$value,
    looic = loo_res$value,
    loglik = loglik
  )
}

# 5. ZINB-only model
# Two independent univariate CDPMMs for (b1 zero, b2 count); no b1-b2 dependence
############################################################
fit_zinb_model <- function(dat, chain = 5000, burn = 2000, thin = 5,
                           G = G_mix) {
  id <- dat$id
  N <- dat$N
  n <- dat$n
  nis <- dat$nis
  
  Xz <- dat$X_zero
  Xc <- dat$X_count
  y1 <- dat$y1
  y2 <- dat$y2
  
  pz <- ncol(Xz)
  pc <- ncol(Xc)
  
  alpha0 <- rep(0, pz)
  beta0 <- rep(0, pc)
  T0a <- diag(0.001, pz)
  T0b <- diag(0.001, pc)
  
  alpha <- rep(0, pz)
  beta <- rep(0, pc)
  r <- 1.0
  
  # independent CDPMM random effects (separate processes for zero and count)
  b1 <- rnorm(n)
  b2 <- rnorm(n)
  cdp1 <- init_cdpmm_1d(n, G = G)
  cdp2 <- init_cdpmm_1d(n, G = G)
  id_index <- split(seq_len(N), id)
  
  save_every <- floor((chain - burn) / thin)
  
  u_est <- as.numeric(y1 > 0)
  Alpha_store <- matrix(NA, save_every, pz)
  Beta_store <- matrix(NA, save_every, pc)
  R_store <- rep(NA, save_every)
  sigma2_b2_store <- rep(NA, save_every)  # zero RE (b1)
  sigma2_b3_store <- rep(NA, save_every)  # count RE (b2)
  loglik <- matrix(NA, save_every, N)
  zero_rep <- rep(NA, save_every)
  cor_rep <- rep(NA, save_every)

  for (iter in 1:chain) {
    eta_zero <- as.numeric(Xz %*% alpha + rep(b1, times = nis))
    eta_zero <- pmin(pmax(eta_zero, -10), 10)
    pi_s <- inv_logit(eta_zero)
    pi_s <- pmin(pmax(pi_s, 1e-10), 1 - 1e-10)
    
    eta_count <- as.numeric(Xc %*% beta + rep(b2, times = nis))
    eta_count <- pmin(pmax(eta_count, -10), 10)
    phi_s <- inv_logit(eta_count)
    phi_s <- pmin(pmax(phi_s, 1e-10), 1 - 1e-10)
    
    mu_s <- r * phi_s / (1 - phi_s)
    mu_s <- pmax(mu_s, 1e-10)
    
    q_nb0 <- dnbinom_zero(r, mu_s)
    
    for (i in 1:N) {
      if (y1[i] == 0) {
        log_p1 <- log(pi_s[i]) + log(q_nb0[i])
        log_p0 <- log(1 - pi_s[i])
        theta <- exp(log_p1 - log_sum_exp(c(log_p1, log_p0)))
        theta <- pmin(pmax(theta, 1e-10), 1 - 1e-10)
        u_est[i] <- rbinom(1, 1, theta)
      } else {
        u_est[i] <- 1
      }
    }
    
    pos_idx <- which(u_est == 1)
    n_pos <- length(pos_idx)
    
    omega_zero <- rpg(N, 1, eta_zero)
    z_a <- (u_est - 0.5) / pmax(omega_zero, 1e-10)
    
    if (n_pos > 0) {
      y1_pos <- y1[pos_idx]
      Xc_pos <- Xc[pos_idx, , drop = FALSE]
      b2_pos <- rep(b2, times = nis)[pos_idx]
      eta_count_pos <- as.numeric(Xc_pos %*% beta + b2_pos)
      eta_count_pos <- pmin(pmax(eta_count_pos, -10), 10)
      w_count <- rpg(n_pos, y1_pos + r, eta_count_pos)
      z_b <- (y1_pos - r) / (2 * pmax(w_count, 1e-10))
    }
    
    V_a <- solve(crossprod(sqrt(omega_zero) * Xz) + T0a)
    m_a <- V_a %*% (T0a %*% alpha0 + crossprod(Xz, omega_zero * (z_a - rep(b1, times = nis))))
    alpha <- as.numeric(rmvnorm(1, m_a, V_a))
    
    if (n_pos > 0) {
      V_b <- solve(crossprod(sqrt(w_count) * Xc_pos) + T0b)
      m_b <- V_b %*% (T0b %*% beta0 + crossprod(Xc_pos, w_count * (z_b - b2_pos)))
      beta <- as.numeric(rmvnorm(1, m_b, V_b))
    }
    
    if (n_pos > 0) {
      k_crp <- rep(0, n_pos)
      for (j in 1:n_pos) {
        if (y1_pos[j] > 0) {
          probs <- r / (r + 0:(y1_pos[j] - 1))
          probs <- pmin(pmax(probs, 1e-10), 1 - 1e-10)
          k_crp[j] <- sum(rbinom(y1_pos[j], 1, probs))
        }
      }
      
      eta_count_current <- as.numeric(Xc_pos %*% beta + b2_pos)
      psi_current <- exp(eta_count_current) / (1 + exp(eta_count_current))
      psi_current <- pmin(pmax(psi_current, 1e-10), 1 - 1e-10)
      
      log_term <- sum(log(1 - psi_current))
      if (!is.finite(log_term)) log_term <- 0
      
      r <- rgamma(1, shape = 0.01 + sum(k_crp), rate = max(0.01 - log_term, 1e-10))
    }
    
    # independent CDPMM updates for b1 (zero) and b2 (count)
    for (j in 1:n) {
      idx <- id_index[[as.character(j)]]
      if (is.null(idx)) idx <- id_index[[j]]
      Xz_j <- Xz[idx, , drop = FALSE]
      oz <- omega_zero[idx]
      za <- z_a[idx]
      g1 <- cdp1$L[j]
      prior_prec1 <- 1 / max(cdp1$omega[g1], 1e-8)
      prior_mean1 <- cdp1$mu[g1]
      post_var_b1 <- 1 / (prior_prec1 + sum(oz))
      post_mean_b1 <- post_var_b1 * (prior_prec1 * prior_mean1 +
                                       sum(oz * (za - as.numeric(Xz_j %*% alpha))))
      b1[j] <- rnorm(1, post_mean_b1, sqrt(post_var_b1))
    }
    
    for (j in 1:n) {
      idx <- id_index[[as.character(j)]]
      if (is.null(idx)) idx <- id_index[[j]]
      g2 <- cdp2$L[j]
      prior_prec2 <- 1 / max(cdp2$omega[g2], 1e-8)
      prior_mean2 <- cdp2$mu[g2]
      idx_pos <- idx[idx %in% pos_idx]
      if (length(idx_pos) > 0 && n_pos > 0) {
        pos_in_pos <- match(idx_pos, pos_idx)
        wg_j <- w_count[pos_in_pos]
        zb_j <- z_b[pos_in_pos]
        Xc_pos_j <- Xc[idx_pos, , drop = FALSE]
        post_var_b2 <- 1 / (prior_prec2 + sum(wg_j))
        post_mean_b2 <- post_var_b2 * (prior_prec2 * prior_mean2 +
                                         sum(wg_j * (zb_j - as.numeric(Xc_pos_j %*% beta))))
        b2[j] <- rnorm(1, post_mean_b2, sqrt(post_var_b2))
      } else {
        b2[j] <- rnorm(1, prior_mean2, sqrt(1 / prior_prec2))
      }
    }
    
    cdp1 <- update_cdpmm_1d(cdp1, b1)
    cdp2 <- update_cdpmm_1d(cdp2, b2)
    
    if (iter > burn && ((iter - burn) %% thin == 0)) {
      s <- (iter - burn) / thin
      Alpha_store[s, ] <- alpha
      Beta_store[s, ] <- beta
      R_store[s] <- r
      sigma2_b2_store[s] <- cdp1$sigma2
      sigma2_b3_store[s] <- cdp2$sigma2
      
      eta_zero_s <- as.numeric(Xz %*% alpha + rep(b1, times = nis))
      pi_s <- inv_logit(eta_zero_s)
      
      eta_count_s <- as.numeric(Xc %*% beta + rep(b2, times = nis))
      phi_s <- inv_logit(eta_count_s)
      
      mu_s <- r * phi_s / (1 - phi_s)
      mu_s <- pmax(mu_s, 1e-10)
      
      loglik[s, ] <- ifelse(
        y1 == 0,
        log((1 - pi_s) + pi_s * dnbinom(y1, size = r, mu = mu_s)),
        log(pi_s) + dnbinom(y1, size = r, mu = mu_s, log = TRUE)
      )
      
      zero_rep[s] <- mean(y1 == 0)
      cor_rep[s] <- safe_cor(y2, as.numeric(y1 > 0))
    }
  }
  
  alpha_mean <- colMeans(Alpha_store, na.rm = TRUE)
  beta_mean <- colMeans(Beta_store, na.rm = TRUE)
  r_mean <- mean(R_store, na.rm = TRUE)
  
  theta_bar_loglik_fun <- function() {
    eta_zero <- as.numeric(Xz %*% alpha_mean)
    pi_s <- inv_logit(eta_zero)
    
    eta_count <- as.numeric(Xc %*% beta_mean)
    phi_s <- inv_logit(eta_count)
    
    mu_s <- r_mean * phi_s / (1 - phi_s)
    mu_s <- pmax(mu_s, 1e-10)
    
    ll <- ifelse(
      y1 == 0,
      log((1 - pi_s) + pi_s * dnbinom(y1, size = r_mean, mu = mu_s)),
      log(pi_s) + dnbinom(y1, size = r_mean, mu = mu_s, log = TRUE)
    )
    -2 * sum(ll)
  }
  
  dic_res <- compute_dic_from_loglik(loglik, theta_bar_loglik_fun)
  waic_res <- compute_waic(loglik)
  loo_res  <- compute_looic(loglik)
  
  list(
    alpha_samples = Alpha_store,
    beta_samples = Beta_store,
    r_samples = R_store,
    sigma2_b2_samples = sigma2_b2_store,
    sigma2_b3_samples = sigma2_b3_store,
    alpha_est = alpha_mean,
    beta_est = beta_mean,
    r_est = r_mean,
    sigma2_b2_est = mean(sigma2_b2_store, na.rm = TRUE),
    sigma2_b3_est = mean(sigma2_b3_store, na.rm = TRUE),
    sigma2_b2_ci_lower = quantile(sigma2_b2_store, 0.025, na.rm = TRUE),
    sigma2_b2_ci_upper = quantile(sigma2_b2_store, 0.975, na.rm = TRUE),
    sigma2_b3_ci_lower = quantile(sigma2_b3_store, 0.025, na.rm = TRUE),
    sigma2_b3_ci_upper = quantile(sigma2_b3_store, 0.975, na.rm = TRUE),
    alpha_ci_lower = apply(Alpha_store, 2, function(x) quantile(x, 0.025, na.rm = TRUE)),
    alpha_ci_upper = apply(Alpha_store, 2, function(x) quantile(x, 0.975, na.rm = TRUE)),
    beta_ci_lower = apply(Beta_store, 2, function(x) quantile(x, 0.025, na.rm = TRUE)),
    beta_ci_upper = apply(Beta_store, 2, function(x) quantile(x, 0.975, na.rm = TRUE)),
    r_ci_lower = quantile(R_store, 0.025, na.rm = TRUE),
    r_ci_upper = quantile(R_store, 0.975, na.rm = TRUE),
    dic = dic_res$DIC,
    waic = waic_res$value,
    looic = loo_res$value,
    loglik = loglik,
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
  "\n=== Prog1: CDPMM joint vs separate | 重复模拟 ", n_sim, " 次",
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
                           joint_fit <- fit_joint_model(
                             dat, chain = chain, burn = burn, thin = thin,
                             delta_min = delta_min, delta_max = delta_max
                           )
                           ordinal_fit <- fit_ordinal_model(
                             dat, chain = chain, burn = burn, thin = thin,
                             delta_min = delta_min, delta_max = delta_max
                           )
                           zinb_fit <- fit_zinb_model(
                             dat, chain = chain, burn = burn, thin = thin
                           )
                           list(dat = dat, joint = joint_fit,
                                ordinal = ordinal_fit, zinb = zinb_fit)
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


cat("True rho (zero-count, zero-ordinal, count-ordinal):",
    paste(round(rho_true, 3), collapse = ", "), "\n")

# Separate-model variance truths = corresponding diagonal of population Sigma
# RE ordering: (zero, count, ordinal) = (b1, b2, b3) = Sigma[1,1], [2,2], [3,3]
sigma2_b1_true <- Sigma_true[1, 1]
sigma2_b2_true <- Sigma_true[2, 2]
sigma2_b3_true <- Sigma_true[3, 3]

############################################################
# 9. Collect estimates across replications
############################################################
alpha_est_joint <- matrix(NA, n_success, length(alpha_true))
alpha_ci_lower_joint <- matrix(NA, n_success, length(alpha_true))
alpha_ci_upper_joint <- matrix(NA, n_success, length(alpha_true))

delta_est_joint <- matrix(NA, n_success, C - 1)
delta_ci_lower_joint <- matrix(NA, n_success, C - 1)
delta_ci_upper_joint <- matrix(NA, n_success, C - 1)
delta_true_mat <- matrix(NA, n_success, C - 1)

beta_est_joint <- matrix(NA, n_success, length(beta_true))
beta_ci_lower_joint <- matrix(NA, n_success, length(beta_true))
beta_ci_upper_joint <- matrix(NA, n_success, length(beta_true))

gamma_est_joint <- matrix(NA, n_success, length(gamma_true))
gamma_ci_lower_joint <- matrix(NA, n_success, length(gamma_true))
gamma_ci_upper_joint <- matrix(NA, n_success, length(gamma_true))

r_est_joint <- rep(NA, n_success)
r_ci_lower_joint <- rep(NA, n_success)
r_ci_upper_joint <- rep(NA, n_success)

Sigma_est_joint <- array(NA, dim = c(n_success, 3, 3))
Sigma_ci_lower_joint <- array(NA, dim = c(n_success, 3, 3))
Sigma_ci_upper_joint <- array(NA, dim = c(n_success, 3, 3))

rho_est_joint <- matrix(NA, n_success, 3)
rho_ci_lower_joint <- matrix(NA, n_success, 3)
rho_ci_upper_joint <- matrix(NA, n_success, 3)

gamma_est_ord <- matrix(NA, n_success, length(gamma_true))
gamma_ci_lower_ord <- matrix(NA, n_success, length(gamma_true))
gamma_ci_upper_ord <- matrix(NA, n_success, length(gamma_true))

delta_est_ord <- matrix(NA, n_success, C - 1)
delta_ci_lower_ord <- matrix(NA, n_success, C - 1)
delta_ci_upper_ord <- matrix(NA, n_success, C - 1)

sigma2_est_ord <- rep(NA, n_success)
sigma2_ci_lower_ord <- rep(NA, n_success)
sigma2_ci_upper_ord <- rep(NA, n_success)

alpha_est_zinb <- matrix(NA, n_success, length(alpha_true))
alpha_ci_lower_zinb <- matrix(NA, n_success, length(alpha_true))
alpha_ci_upper_zinb <- matrix(NA, n_success, length(alpha_true))

beta_est_zinb <- matrix(NA, n_success, length(beta_true))
beta_ci_lower_zinb <- matrix(NA, n_success, length(beta_true))
beta_ci_upper_zinb <- matrix(NA, n_success, length(beta_true))

r_est_zinb <- rep(NA, n_success)
r_ci_lower_zinb <- rep(NA, n_success)
r_ci_upper_zinb <- rep(NA, n_success)

sigma2_b2_est_zinb <- rep(NA, n_success)
sigma2_b2_ci_lower_zinb <- rep(NA, n_success)
sigma2_b2_ci_upper_zinb <- rep(NA, n_success)

sigma2_b3_est_zinb <- rep(NA, n_success)
sigma2_b3_ci_lower_zinb <- rep(NA, n_success)
sigma2_b3_ci_upper_zinb <- rep(NA, n_success)

model_compare_joint <- matrix(NA, n_success, 3)
model_compare_ordinal <- matrix(NA, n_success, 3)
model_compare_zinb <- matrix(NA, n_success, 3)

ppc_joint_zero <- rep(NA, n_success)
ppc_zinb_zero <- rep(NA, n_success)
ppc_joint_cor <- rep(NA, n_success)
ppc_zinb_cor <- rep(NA, n_success)

vr_alpha_list <- vector("list", n_success)
vr_beta_list <- vector("list", n_success)
vr_gamma_list <- vector("list", n_success)
vr_r_list <- rep(NA, n_success)
vr_sigma2_ord_list <- rep(NA, n_success)
vr_sigma2_b2_list <- rep(NA, n_success)
vr_sigma2_b3_list <- rep(NA, n_success)

for (i in seq_len(n_success)) {
  res <- sim_results[[i]]

  alpha_est_joint[i, ] <- res$joint$alpha_est
  alpha_ci_lower_joint[i, ] <- res$joint$alpha_ci_lower
  alpha_ci_upper_joint[i, ] <- res$joint$alpha_ci_upper

  delta_est_joint[i, ] <- res$joint$delta_est
  delta_ci_lower_joint[i, ] <- res$joint$delta_ci_lower
  delta_ci_upper_joint[i, ] <- res$joint$delta_ci_upper
  delta_true_mat[i, ] <- res$dat$delta_true

  beta_est_joint[i, ] <- res$joint$beta_est
  beta_ci_lower_joint[i, ] <- res$joint$beta_ci_lower
  beta_ci_upper_joint[i, ] <- res$joint$beta_ci_upper

  gamma_est_joint[i, ] <- res$joint$gamma_est
  gamma_ci_lower_joint[i, ] <- res$joint$gamma_ci_lower
  gamma_ci_upper_joint[i, ] <- res$joint$gamma_ci_upper

  r_est_joint[i] <- res$joint$r_est
  r_ci_lower_joint[i] <- res$joint$r_ci_lower
  r_ci_upper_joint[i] <- res$joint$r_ci_upper

  Sigma_est_joint[i, , ] <- res$joint$Sigma_est
  Sigma_ci_lower_joint[i, , ] <- res$joint$Sigma_ci_lower
  Sigma_ci_upper_joint[i, , ] <- res$joint$Sigma_ci_upper

  rho_est_joint[i, ] <- res$joint$Rho_est
  rho_ci_lower_joint[i, ] <- res$joint$Rho_ci_lower
  rho_ci_upper_joint[i, ] <- res$joint$Rho_ci_upper

  gamma_est_ord[i, ] <- res$ordinal$gamma_est
  gamma_ci_lower_ord[i, ] <- res$ordinal$gamma_ci_lower
  gamma_ci_upper_ord[i, ] <- res$ordinal$gamma_ci_upper

  delta_est_ord[i, ] <- res$ordinal$delta_est
  delta_ci_lower_ord[i, ] <- res$ordinal$delta_ci_lower
  delta_ci_upper_ord[i, ] <- res$ordinal$delta_ci_upper

  sigma2_est_ord[i] <- res$ordinal$sigma2_est
  sigma2_ci_lower_ord[i] <- res$ordinal$sigma2_ci_lower
  sigma2_ci_upper_ord[i] <- res$ordinal$sigma2_ci_upper

  alpha_est_zinb[i, ] <- res$zinb$alpha_est
  alpha_ci_lower_zinb[i, ] <- res$zinb$alpha_ci_lower
  alpha_ci_upper_zinb[i, ] <- res$zinb$alpha_ci_upper

  beta_est_zinb[i, ] <- res$zinb$beta_est
  beta_ci_lower_zinb[i, ] <- res$zinb$beta_ci_lower
  beta_ci_upper_zinb[i, ] <- res$zinb$beta_ci_upper

  r_est_zinb[i] <- res$zinb$r_est
  r_ci_lower_zinb[i] <- res$zinb$r_ci_lower
  r_ci_upper_zinb[i] <- res$zinb$r_ci_upper

  sigma2_b2_est_zinb[i] <- res$zinb$sigma2_b2_est
  sigma2_b2_ci_lower_zinb[i] <- res$zinb$sigma2_b2_ci_lower
  sigma2_b2_ci_upper_zinb[i] <- res$zinb$sigma2_b2_ci_upper

  sigma2_b3_est_zinb[i] <- res$zinb$sigma2_b3_est
  sigma2_b3_ci_lower_zinb[i] <- res$zinb$sigma2_b3_ci_lower
  sigma2_b3_ci_upper_zinb[i] <- res$zinb$sigma2_b3_ci_upper

  model_compare_joint[i, ] <- c(res$joint$dic, res$joint$waic, res$joint$looic)
  model_compare_ordinal[i, ] <- c(res$ordinal$dic, res$ordinal$waic, res$ordinal$looic)
  model_compare_zinb[i, ] <- c(res$zinb$dic, res$zinb$waic, res$zinb$looic)

  ppc_joint_zero[i] <- res$joint$ppc$rep_zero_prop_mean
  ppc_zinb_zero[i] <- res$zinb$ppc$rep_zero_prop_mean

  ppc_joint_cor[i] <- res$joint$ppc$rep_cor_y1_y2pos_mean
  ppc_zinb_cor[i] <- res$zinb$ppc$rep_cor_y1_y2pos_mean

  vr_alpha_list[[i]] <- safe_col_var(res$zinb$alpha_samples) / safe_col_var(res$joint$alpha_samples)
  vr_beta_list[[i]] <- safe_col_var(res$zinb$beta_samples) / safe_col_var(res$joint$beta_samples)
  vr_gamma_list[[i]] <- safe_col_var(res$ordinal$gamma_samples) / safe_col_var(res$joint$gamma_samples)
  vr_r_list[i] <- stats::var(res$zinb$r_samples, na.rm = TRUE) / stats::var(res$joint$r_samples, na.rm = TRUE)
  # sigma2: zero / count / ordinal = Sigma[1,1], [2,2], [3,3]
  vr_sigma2_ord_list[i] <- stats::var(res$ordinal$sigma2_samples, na.rm = TRUE) / stats::var(res$joint$Sigma_samples[, 9], na.rm = TRUE)
  vr_sigma2_b2_list[i] <- stats::var(res$zinb$sigma2_b2_samples, na.rm = TRUE) / stats::var(res$joint$Sigma_samples[, 1], na.rm = TRUE)
  vr_sigma2_b3_list[i] <- stats::var(res$zinb$sigma2_b3_samples, na.rm = TRUE) / stats::var(res$joint$Sigma_samples[, 5], na.rm = TRUE)
}

############################################################
# 10. Parameter summaries
############################################################
cat("\n================ 参数估计汇总（联合模型） ================\n")

alpha_stats_joint <- calculate_stats(alpha_est_joint, alpha_true, alpha_ci_lower_joint, alpha_ci_upper_joint)
alpha_summary_joint <- data.frame(
  Parameter = paste0("alpha", seq_along(alpha_true)),
  True_Value = alpha_true,
  Bias = round(alpha_stats_joint$bias, 4),
  RMSE = round(alpha_stats_joint$rmse, 4),
  CI_Lower_Mean = round(alpha_stats_joint$ci_lower_mean, 4),
  CI_Upper_Mean = round(alpha_stats_joint$ci_upper_mean, 4),
  CP = round(alpha_stats_joint$cp, 3)
)
print(alpha_summary_joint, row.names = FALSE)

delta_stats_joint <- calculate_stats(delta_est_joint, delta_true, delta_ci_lower_joint, delta_ci_upper_joint)
delta_summary_joint <- data.frame(
  Parameter = paste0("delta", seq_len(C - 1)),
  True_Value = delta_true,
  Bias = round(delta_stats_joint$bias, 4),
  RMSE = round(delta_stats_joint$rmse, 4),
  CI_Lower_Mean = round(delta_stats_joint$ci_lower_mean, 4),
  CI_Upper_Mean = round(delta_stats_joint$ci_upper_mean, 4),
  CP = round(delta_stats_joint$cp, 3)
)
print(delta_summary_joint, row.names = FALSE)

beta_stats_joint <- calculate_stats(beta_est_joint, beta_true, beta_ci_lower_joint, beta_ci_upper_joint)
beta_summary_joint <- data.frame(
  Parameter = paste0("beta", seq_along(beta_true)),
  True_Value = beta_true,
  Bias = round(beta_stats_joint$bias, 4),
  RMSE = round(beta_stats_joint$rmse, 4),
  CI_Lower_Mean = round(beta_stats_joint$ci_lower_mean, 4),
  CI_Upper_Mean = round(beta_stats_joint$ci_upper_mean, 4),
  CP = round(beta_stats_joint$cp, 3)
)
print(beta_summary_joint, row.names = FALSE)

gamma_stats_joint <- calculate_stats(gamma_est_joint, gamma_true, gamma_ci_lower_joint, gamma_ci_upper_joint)
gamma_summary_joint <- data.frame(
  Parameter = paste0("gamma", seq_along(gamma_true)),
  True_Value = gamma_true,
  Bias = round(gamma_stats_joint$bias, 4),
  RMSE = round(gamma_stats_joint$rmse, 4),
  CI_Lower_Mean = round(gamma_stats_joint$ci_lower_mean, 4),
  CI_Upper_Mean = round(gamma_stats_joint$ci_upper_mean, 4),
  CP = round(gamma_stats_joint$cp, 3)
)
print(gamma_summary_joint, row.names = FALSE)

r_stats_joint <- calculate_stats(
  matrix(r_est_joint, ncol = 1),
  matrix(r_true, ncol = 1),
  matrix(r_ci_lower_joint, ncol = 1),
  matrix(r_ci_upper_joint, ncol = 1)
)
r_summary_joint <- data.frame(
  Parameter = "r",
  True_Value = r_true,
  Bias = round(r_stats_joint$bias, 4),
  RMSE = round(r_stats_joint$rmse, 4),
  CI_Lower_Mean = round(r_stats_joint$ci_lower_mean, 4),
  CI_Upper_Mean = round(r_stats_joint$ci_upper_mean, 4),
  CP = round(r_stats_joint$cp, 3)
)
print(r_summary_joint, row.names = FALSE)

cat("\n================ 随机效应方差汇总（联合模型） ================\n")
Sigma_bias <- matrix(NA, 3, 3)
Sigma_rmse <- matrix(NA, 3, 3)
Sigma_ci_lower_mean <- matrix(NA, 3, 3)
Sigma_ci_upper_mean <- matrix(NA, 3, 3)
Sigma_cp <- matrix(NA, 3, 3)
Sigma_summary_df <- data.frame()
for (i in 1:3) {
  for (j in 1:3) {
    Sigma_bias[i, j] <- mean(Sigma_est_joint[, i, j] - Sigma_true[i, j], na.rm = TRUE)
    Sigma_rmse[i, j] <- sqrt(mean((Sigma_est_joint[, i, j] - Sigma_true[i, j])^2, na.rm = TRUE))
    Sigma_ci_lower_mean[i, j] <- mean(Sigma_ci_lower_joint[, i, j], na.rm = TRUE)
    Sigma_ci_upper_mean[i, j] <- mean(Sigma_ci_upper_joint[, i, j], na.rm = TRUE)
    Sigma_cp[i, j] <- mean(Sigma_ci_lower_joint[, i, j] <= Sigma_true[i, j] & Sigma_true[i, j] <= Sigma_ci_upper_joint[, i, j], na.rm = TRUE)
    Sigma_summary_df <- rbind(Sigma_summary_df, data.frame(
      Parameter = paste0("Sigma[", i, ",", j, "]"),
      True_Value = round(Sigma_true[i, j], 3),
      Bias = round(Sigma_bias[i, j], 4),
      RMSE = round(Sigma_rmse[i, j], 4),
      CI_Lower_Mean = round(Sigma_ci_lower_mean[i, j], 4),
      CI_Upper_Mean = round(Sigma_ci_upper_mean[i, j], 4),
      CP = round(Sigma_cp[i, j], 3)
    ))
  }
}
print(Sigma_summary_df, row.names = FALSE)

rho_stats_joint <- calculate_stats(rho_est_joint, rho_true, rho_ci_lower_joint, rho_ci_upper_joint)
rho_summary_joint <- data.frame(
  Parameter = c("rho12", "rho13", "rho23"),
  True_Value = round(rho_true, 3),
  Bias = round(rho_stats_joint$bias, 4),
  RMSE = round(rho_stats_joint$rmse, 4),
  CI_Lower_Mean = round(rho_stats_joint$ci_lower_mean, 4),
  CI_Upper_Mean = round(rho_stats_joint$ci_upper_mean, 4),
  CP = round(rho_stats_joint$cp, 3)
)
print(rho_summary_joint, row.names = FALSE)

cat("\n================ 参数估计汇总（有序单独模型） ================\n")
gamma_stats_ord <- calculate_stats(gamma_est_ord, gamma_true, gamma_ci_lower_ord, gamma_ci_upper_ord)
gamma_summary_ord <- data.frame(
  Parameter = paste0("gamma", seq_along(gamma_true)),
  True_Value = gamma_true,
  Bias = round(gamma_stats_ord$bias, 4),
  RMSE = round(gamma_stats_ord$rmse, 4),
  CI_Lower_Mean = round(gamma_stats_ord$ci_lower_mean, 4),
  CI_Upper_Mean = round(gamma_stats_ord$ci_upper_mean, 4),
  CP = round(gamma_stats_ord$cp, 3)
)
print(gamma_summary_ord, row.names = FALSE)

delta_stats_ord <- calculate_stats(delta_est_ord, delta_true, delta_ci_lower_ord, delta_ci_upper_ord)
delta_summary_ord <- data.frame(
  Parameter = paste0("delta", seq_len(C - 1)),
  True_Value = delta_true,
  Bias = round(delta_stats_ord$bias, 4),
  RMSE = round(delta_stats_ord$rmse, 4),
  CI_Lower_Mean = round(delta_stats_ord$ci_lower_mean, 4),
  CI_Upper_Mean = round(delta_stats_ord$ci_upper_mean, 4),
  CP = round(delta_stats_ord$cp, 3)
)
print(delta_summary_ord, row.names = FALSE)

sigma_stats_ord <- calculate_stats(
  matrix(sigma2_est_ord, ncol = 1),
  matrix(sigma2_b3_true, ncol = 1),
  matrix(sigma2_ci_lower_ord, ncol = 1),
  matrix(sigma2_ci_upper_ord, ncol = 1)
)
sigma_summary_ord <- data.frame(
  Parameter = "sigma2_b3",
  True_Value = sigma2_b3_true,
  Bias = round(sigma_stats_ord$bias, 4),
  RMSE = round(sigma_stats_ord$rmse, 4),
  CI_Lower_Mean = round(sigma_stats_ord$ci_lower_mean, 4),
  CI_Upper_Mean = round(sigma_stats_ord$ci_upper_mean, 4),
  CP = round(sigma_stats_ord$cp, 3)
)
print(sigma_summary_ord, row.names = FALSE)

cat("\n================ 参数估计汇总（ZINB单独模型） ================\n")
alpha_stats_zinb <- calculate_stats(alpha_est_zinb, alpha_true, alpha_ci_lower_zinb, alpha_ci_upper_zinb)
alpha_summary_zinb <- data.frame(
  Parameter = paste0("alpha", seq_along(alpha_true)),
  True_Value = alpha_true,
  Bias = round(alpha_stats_zinb$bias, 4),
  RMSE = round(alpha_stats_zinb$rmse, 4),
  CI_Lower_Mean = round(alpha_stats_zinb$ci_lower_mean, 4),
  CI_Upper_Mean = round(alpha_stats_zinb$ci_upper_mean, 4),
  CP = round(alpha_stats_zinb$cp, 3)
)
print(alpha_summary_zinb, row.names = FALSE)

beta_stats_zinb <- calculate_stats(beta_est_zinb, beta_true, beta_ci_lower_zinb, beta_ci_upper_zinb)
beta_summary_zinb <- data.frame(
  Parameter = paste0("beta", seq_along(beta_true)),
  True_Value = beta_true,
  Bias = round(beta_stats_zinb$bias, 4),
  RMSE = round(beta_stats_zinb$rmse, 4),
  CI_Lower_Mean = round(beta_stats_zinb$ci_lower_mean, 4),
  CI_Upper_Mean = round(beta_stats_zinb$ci_upper_mean, 4),
  CP = round(beta_stats_zinb$cp, 3)
)
print(beta_summary_zinb, row.names = FALSE)

r_stats_zinb <- calculate_stats(
  matrix(r_est_zinb, ncol = 1),
  matrix(r_true, ncol = 1),
  matrix(r_ci_lower_zinb, ncol = 1),
  matrix(r_ci_upper_zinb, ncol = 1)
)
r_summary_zinb <- data.frame(
  Parameter = "r",
  True_Value = r_true,
  Bias = round(r_stats_zinb$bias, 4),
  RMSE = round(r_stats_zinb$rmse, 4),
  CI_Lower_Mean = round(r_stats_zinb$ci_lower_mean, 4),
  CI_Upper_Mean = round(r_stats_zinb$ci_upper_mean, 4),
  CP = round(r_stats_zinb$cp, 3)
)
print(r_summary_zinb, row.names = FALSE)

sigma2_b2_stats_zinb <- calculate_stats(
  matrix(sigma2_b2_est_zinb, ncol = 1),
  matrix(sigma2_b1_true, ncol = 1),
  matrix(sigma2_b2_ci_lower_zinb, ncol = 1),
  matrix(sigma2_b2_ci_upper_zinb, ncol = 1)
)
sigma2_b2_summary_zinb <- data.frame(
  Parameter = "sigma2_b1",
  True_Value = sigma2_b1_true,
  Bias = round(sigma2_b2_stats_zinb$bias, 4),
  RMSE = round(sigma2_b2_stats_zinb$rmse, 4),
  CI_Lower_Mean = round(sigma2_b2_stats_zinb$ci_lower_mean, 4),
  CI_Upper_Mean = round(sigma2_b2_stats_zinb$ci_upper_mean, 4),
  CP = round(sigma2_b2_stats_zinb$cp, 3)
)
print(sigma2_b2_summary_zinb, row.names = FALSE)

sigma2_b3_stats_zinb <- calculate_stats(
  matrix(sigma2_b3_est_zinb, ncol = 1),
  matrix(sigma2_b2_true, ncol = 1),
  matrix(sigma2_b3_ci_lower_zinb, ncol = 1),
  matrix(sigma2_b3_ci_upper_zinb, ncol = 1)
)
sigma2_b3_summary_zinb <- data.frame(
  Parameter = "sigma2_b2",
  True_Value = sigma2_b2_true,
  Bias = round(sigma2_b3_stats_zinb$bias, 4),
  RMSE = round(sigma2_b3_stats_zinb$rmse, 4),
  CI_Lower_Mean = round(sigma2_b3_stats_zinb$ci_lower_mean, 4),
  CI_Upper_Mean = round(sigma2_b3_stats_zinb$ci_upper_mean, 4),
  CP = round(sigma2_b3_stats_zinb$cp, 3)
)
print(sigma2_b3_summary_zinb, row.names = FALSE)

############################################################
# 11. Model comparison summaries
############################################################
cat("\n================ 模型比较汇总 ================\n")
model_compare_summary <- data.frame(
  Model = c("Joint", "Ordinal-only", "ZINB-only"),
  DIC_Mean = c(
    mean(model_compare_joint[, 1], na.rm = TRUE),
    mean(model_compare_ordinal[, 1], na.rm = TRUE),
    mean(model_compare_zinb[, 1], na.rm = TRUE)
  ),
  WAIC_Mean = c(
    mean(model_compare_joint[, 2], na.rm = TRUE),
    mean(model_compare_ordinal[, 2], na.rm = TRUE),
    mean(model_compare_zinb[, 2], na.rm = TRUE)
  ),
  LOOIC_Mean = c(
    mean(model_compare_joint[, 3], na.rm = TRUE),
    mean(model_compare_ordinal[, 3], na.rm = TRUE),
    mean(model_compare_zinb[, 3], na.rm = TRUE)
  )
)
print(model_compare_summary, row.names = FALSE)

joint_better_than_sum_dic <- mean(model_compare_joint[, 1] < (model_compare_ordinal[, 1] + model_compare_zinb[, 1]), na.rm = TRUE) * 100
joint_better_than_sum_waic <- mean(model_compare_joint[, 2] < (model_compare_ordinal[, 2] + model_compare_zinb[, 2]), na.rm = TRUE) * 100
joint_better_than_sum_looic <- mean(model_compare_joint[, 3] < (model_compare_ordinal[, 3] + model_compare_zinb[, 3]), na.rm = TRUE) * 100
criteria_winrate <- data.frame(
  Comparison = "CDPMM_joint_better_than_sum_of_separate",
  DIC_WinPct = round(joint_better_than_sum_dic, 1),
  WAIC_WinPct = round(joint_better_than_sum_waic, 1),
  LOOIC_WinPct = round(joint_better_than_sum_looic, 1),
  DIC_Mean_Joint = mean(model_compare_joint[, 1], na.rm = TRUE),
  DIC_Mean_SepSum = mean(model_compare_ordinal[, 1] + model_compare_zinb[, 1], na.rm = TRUE),
  WAIC_Mean_Joint = mean(model_compare_joint[, 2], na.rm = TRUE),
  WAIC_Mean_SepSum = mean(model_compare_ordinal[, 2] + model_compare_zinb[, 2], na.rm = TRUE),
  LOOIC_Mean_Joint = mean(model_compare_joint[, 3], na.rm = TRUE),
  LOOIC_Mean_SepSum = mean(model_compare_ordinal[, 3] + model_compare_zinb[, 3], na.rm = TRUE)
)
cat("\n================ 准则胜率 (联合模型优于单独模型之和 的百分比) ================\n")
print(criteria_winrate, row.names = FALSE)

############################################################
# 12. Posterior predictive check summaries
############################################################
cat("\n================ 后验预测检验汇总 ================\n")
ppc_summary <- data.frame(
  Model = c("Joint", "Ordinal-only", "ZINB-only"),
  Zero_Prop_Obs = c(
    mean(sim_results[[1]]$dat$y1 == 0),
    NA_real_,
    mean(sim_results[[1]]$dat$y1 == 0)
  ),
  Zero_Prop_Rep_Mean = c(
    mean(ppc_joint_zero, na.rm = TRUE),
    NA_real_,
    mean(ppc_zinb_zero, na.rm = TRUE)
  ),
  Cor_Obs = c(
    safe_cor(sim_results[[1]]$dat$y2, as.numeric(sim_results[[1]]$dat$y1 > 0)),
    safe_cor(sim_results[[1]]$dat$y2, as.numeric(sim_results[[1]]$dat$y1 > 0)),
    safe_cor(sim_results[[1]]$dat$y2, as.numeric(sim_results[[1]]$dat$y1 > 0))
  ),
  Cor_Rep_Mean = c(
    mean(ppc_joint_cor, na.rm = TRUE),
    NA_real_,
    mean(ppc_zinb_cor, na.rm = TRUE)
  )
)
print(ppc_summary, row.names = FALSE)

############################################################
# 13. Variance ratio summaries
############################################################
cat("\n================ 方差比汇总 ================\n")
vr_alpha_mat <- do.call(rbind, vr_alpha_list)
vr_beta_mat <- do.call(rbind, vr_beta_list)
vr_gamma_mat <- do.call(rbind, vr_gamma_list)

vr_alpha_mean <- colMeans(vr_alpha_mat, na.rm = TRUE)
vr_beta_mean <- colMeans(vr_beta_mat, na.rm = TRUE)
vr_gamma_mean <- colMeans(vr_gamma_mat, na.rm = TRUE)
vr_r_mean <- mean(vr_r_list, na.rm = TRUE)
vr_sigma2_ord_mean <- mean(vr_sigma2_ord_list, na.rm = TRUE)
vr_sigma2_b2_mean <- mean(vr_sigma2_b2_list, na.rm = TRUE)
vr_sigma2_b3_mean <- mean(vr_sigma2_b3_list, na.rm = TRUE)

vr_alpha_trim <- apply(vr_alpha_mat, 2, mean, trim = 0.1, na.rm = TRUE)
vr_beta_trim <- apply(vr_beta_mat, 2, mean, trim = 0.1, na.rm = TRUE)
vr_gamma_trim <- apply(vr_gamma_mat, 2, mean, trim = 0.1, na.rm = TRUE)
vr_r_trim <- mean(vr_r_list, na.rm = TRUE, trim = 0.1)
vr_sigma2_ord_trim <- mean(vr_sigma2_ord_list, na.rm = TRUE, trim = 0.1)
vr_sigma2_b2_trim <- mean(vr_sigma2_b2_list, na.rm = TRUE, trim = 0.1)
vr_sigma2_b3_trim <- mean(vr_sigma2_b3_list, na.rm = TRUE, trim = 0.1)

variance_ratio_table <- data.frame(
  Parameter = c(
    paste0("alpha", seq_along(vr_alpha_mean)),
    paste0("beta", seq_along(vr_beta_mean)),
    paste0("gamma", seq_along(vr_gamma_mean)),
    "r",
    "sigma2_b3",
    "sigma2_b1",
    "sigma2_b2"
  ),
  Variance_Ratio_Mean = c(
    vr_alpha_mean,
    vr_beta_mean,
    vr_gamma_mean,
    vr_r_mean,
    vr_sigma2_ord_mean,
    vr_sigma2_b2_mean,
    vr_sigma2_b3_mean
  ),
  Variance_Ratio_TrimMean = c(
    vr_alpha_trim,
    vr_beta_trim,
    vr_gamma_trim,
    vr_r_trim,
    vr_sigma2_ord_trim,
    vr_sigma2_b2_trim,
    vr_sigma2_b3_trim
  )
)
print(variance_ratio_table, row.names = FALSE)

############################################################
# Joint model average CP and RMSE (Sigma diagonal only)
############################################################
joint_cp_values <- c(
  alpha_stats_joint$cp,
  delta_stats_joint$cp,
  beta_stats_joint$cp,
  gamma_stats_joint$cp,
  r_stats_joint$cp,
  diag(Sigma_cp),
  rho_stats_joint$cp
)

joint_rmse_values <- c(
  alpha_stats_joint$rmse,
  delta_stats_joint$rmse,
  beta_stats_joint$rmse,
  gamma_stats_joint$rmse,
  r_stats_joint$rmse,
  diag(Sigma_rmse),
  rho_stats_joint$rmse
)

joint_avg_cp <- mean(joint_cp_values, na.rm = TRUE)
joint_avg_rmse <- mean(joint_rmse_values, na.rm = TRUE)

cat("\n================ 联合模型平均指标（Sigma只取对角元） ================\n")
cat("Average CP  =", round(joint_avg_cp, 3), "\n")
cat("Average RMSE =", round(joint_avg_rmse, 3), "\n")
cat("Average EFF (Mean) =", round(mean(variance_ratio_table$Variance_Ratio_Mean), 3), "\n")
cat("Average EFF (TrimMean) =", round(mean(variance_ratio_table$Variance_Ratio_TrimMean), 3), "\n")

############################################################
# 14. Save summaries
############################################################
tag <- paste0("scen", scenario, "_", re_dist, "_joint_vs_sep_n", n, "_chain", chain)
joint_param_summary <- rbind(
  data.frame(Parameter = paste0("alpha", seq_along(alpha_true)),
             True_Value = alpha_true,
             Bias = round(alpha_stats_joint$bias, 4),
             RMSE = round(alpha_stats_joint$rmse, 4),
             CP = round(alpha_stats_joint$cp, 3)),
  data.frame(Parameter = paste0("delta", seq_along(delta_true)),
             True_Value = delta_true,
             Bias = round(delta_stats_joint$bias, 4),
             RMSE = round(delta_stats_joint$rmse, 4),
             CP = round(delta_stats_joint$cp, 3)),
  data.frame(Parameter = paste0("beta", seq_along(beta_true)),
             True_Value = beta_true,
             Bias = round(beta_stats_joint$bias, 4),
             RMSE = round(beta_stats_joint$rmse, 4),
             CP = round(beta_stats_joint$cp, 3)),
  data.frame(Parameter = paste0("gamma", seq_along(gamma_true)),
             True_Value = gamma_true,
             Bias = round(gamma_stats_joint$bias, 4),
             RMSE = round(gamma_stats_joint$rmse, 4),
             CP = round(gamma_stats_joint$cp, 3)),
  data.frame(Parameter = "r", True_Value = r_true,
             Bias = round(r_stats_joint$bias, 4),
             RMSE = round(r_stats_joint$rmse, 4),
             CP = round(r_stats_joint$cp, 3)),
  data.frame(Parameter = c("Sigma11", "Sigma22", "Sigma33"),
             True_Value = diag(Sigma_true),
             Bias = round(diag(Sigma_bias), 4),
             RMSE = round(diag(Sigma_rmse), 4),
             CP = round(diag(Sigma_cp), 3)),
  data.frame(Parameter = paste0("rho", 1:3),
             True_Value = rho_true,
             Bias = round(rho_stats_joint$bias, 4),
             RMSE = round(rho_stats_joint$rmse, 4),
             CP = round(rho_stats_joint$cp, 3))
)

write.csv(joint_param_summary,
          file.path(OUT_DIR, paste0("joint_param_summary_", tag, ".csv")),
          row.names = FALSE)
write.csv(model_compare_summary,
          file.path(OUT_DIR, paste0("model_compare_", tag, ".csv")),
          row.names = FALSE)
write.csv(criteria_winrate,
          file.path(OUT_DIR, paste0("criteria_winrate_", tag, ".csv")),
          row.names = FALSE)
write.csv(variance_ratio_table,
          file.path(OUT_DIR, paste0("variance_ratio_", tag, ".csv")),
          row.names = FALSE)
write.csv(ppc_summary,
          file.path(OUT_DIR, paste0("ppc_summary_", tag, ".csv")),
          row.names = FALSE)
saveRDS(
  list(
    settings = list(
      scenario = scenario, re_dist = re_dist, re_prior = re_prior,
      n = n, n_sim = n_sim,
      chain = chain, burn = burn, thin = thin, n_success = n_success
    ),
    joint_param_summary = joint_param_summary,
    model_compare_summary = model_compare_summary,
    criteria_winrate = criteria_winrate,
    variance_ratio_table = variance_ratio_table,
    ppc_summary = ppc_summary,
    joint_avg_cp = joint_avg_cp,
    joint_avg_rmse = joint_avg_rmse
  ),
  file.path(OUT_DIR, paste0("sim_summary_", tag, ".rds"))
)
# raw fits are large; keep compact success flag + settings only unless SAVE_RAW=1
if (identical(Sys.getenv("SAVE_RAW", unset = "0"), "1")) {
  saveRDS(sim_results, file.path(OUT_DIR, paste0("sim_raw_", tag, ".rds")))
}
cat("Wrote summaries to: ", OUT_DIR, "\n", sep = "")


} # end SKIP_MAIN_SIM guard
