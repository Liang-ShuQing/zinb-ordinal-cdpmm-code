############################################################
# Program 1 (real data, standalone): CDPMM joint vs separate
# HRS: y1=ZINB(HSPNIT), y2=ordinal(SHLT); RE order (b1,b2,b3)=(zero,count,ord)
# Joint RE prior = multivariate CDPMM only; separate = univariate CDPMM
# Aligned with 模拟研究/prog1_cdpmm_joint_vs_separate.R
# Results -> sibling folder 实际数据分析结果/
############################################################

set.seed(2025)

# ---- defaults (override via env / CLI) ----
n_chains <- 3L
chain_length <- 5000L
burn <- 2000L
thin <- 5L
# Prog1 is CDPMM-only (no Gaussian joint prior)
re_prior <- "cdpmm"
fit_separate <- 1L

G_mix <- 8L
tau_a1 <- 2
tau_a2 <- 4
zeta0_sd2 <- 10
kappa0 <- 1
iw_nu0 <- 6
iw_S0_scale <- 1

delta_min <- -10
delta_max <- 10

covariate_vars <- c("SMOKEV", "HIBP", "DIAB", "LUNG", "BMI")

.env_int <- function(key, default) {
  v <- suppressWarnings(as.integer(Sys.getenv(key, unset = "")))
  if (length(v) == 1L && is.finite(v)) v else default
}
.env_chr <- function(key, default) {
  v <- Sys.getenv(key, unset = "")
  if (nzchar(v)) v else default
}

n_chains <- .env_int("N_CHAINS", n_chains)
chain_length <- .env_int("CHAIN", chain_length)
burn <- .env_int("BURN", burn)
thin <- .env_int("THIN", thin)
# Ignore RE_PRIOR if set to gaussian — this program is CDPMM-only
re_prior <- "cdpmm"
G_mix <- .env_int("G_MIX", G_mix)
fit_separate <- .env_int("FIT_SEPARATE", 1L)
program_tag <- "prog1_joint_vs_sep"

if (!exists("OUT_DIR") || is.null(OUT_DIR) || !nzchar(as.character(OUT_DIR)[1])) {
  args_cli <- commandArgs(trailingOnly = TRUE)
  OUT_DIR <- if (length(args_cli) >= 1L) {
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
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

library(BayesLogit)
library(mvtnorm)
library(MCMCpack)
library(truncnorm)
library(coda)
library(parallel)
library(loo)

use_cdpmm <- identical(re_prior, "cdpmm")

############################################################
# Utility functions (aligned with simulation code.R)
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
  C <- length(delta) + 1L
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
    if (k == start_k && !fix_first) next
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
  if (fix_first) delta_new[1] <- 0
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

safe_col_var <- function(x) {
  apply(x, 2, function(z) stats::var(z, na.rm = TRUE))
}

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

cdpmm_implied_var_1d <- function(pi, mu, omega) {
  v <- sum(pi * (omega + mu^2))
  max(v, 1e-8)
}

init_cdpmm_1d <- function(n, G = G_mix) {
  G <- max(2L, as.integer(G))
  q_joint <- 3L
  nu0 <- max(iw_nu0 - q_joint + 1L, 3L)
  S0 <- iw_S0_scale
  tau <- 1
  nu <- c(rbeta(G - 1L, 1, tau), 1)
  pi_w <- stick_break_weights(nu)
  zeta <- rnorm(1, 0, sqrt(zeta0_sd2))
  omega <- vapply(seq_len(G), function(g) {
    max(as.numeric(safe_riwish(nu0, matrix(S0, 1, 1))), 1e-8)
  }, numeric(1))
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
  st$sigma2 <- cdpmm_implied_var_1d(st$pi_w, st$mu, st$omega)
  st$nclust <- sum(n_g > 0)
  st
}

# Likelihood-invariant location PE: center RE sample mean, shift intercept(s).
# Keeps eta = x'beta + b unchanged; pins intercept vs RE-location ridge.
shift_cdpmm_1d_by <- function(st, b_bar) {
  st$mu_star <- st$mu_star - b_bar
  st$zeta <- st$zeta - b_bar
  mu_bar <- sum(st$pi_w * st$mu_star)
  st$mu <- st$mu_star - mu_bar
  st$sigma2 <- cdpmm_implied_var_1d(st$pi_w, st$mu, st$omega)
  st
}

summarize_posterior <- function(samples, param_names, probs = c(0.025, 0.5, 0.975)) {
  if (!is.matrix(samples)) samples <- matrix(samples, ncol = 1)
  quants <- apply(samples, 2, quantile, probs = probs, na.rm = TRUE)
  data.frame(
    Parameter = param_names,
    Mean = round(colMeans(samples, na.rm = TRUE), 4),
    SD = round(apply(samples, 2, sd, na.rm = TRUE), 4),
    Lower = round(quants[1, ], 4),
    Median = round(quants[2, ], 4),
    Upper = round(quants[3, ], 4),
    row.names = NULL
  )
}

############################################################
# 1. Data preparation (complete-case; paper covariates)
############################################################
cat("=== 数据准备阶段 ===\n")
cat("OUT_DIR =", OUT_DIR, "\n")
cat("re_prior =", re_prior, "  FIT_SEPARATE =", fit_separate, "\n")

data_path <- file.path(getwd(), "Data_model_complete.csv")
if (!file.exists(data_path)) {
  alt <- "D:/Users/lenovo/Desktop/new/实际数据分析/Data_model_complete.csv"
  if (file.exists(alt)) data_path <- alt else stop("找不到 Data_model_complete.csv")
}
dat_raw <- read.csv(data_path)
cat("原始数据维度:", dim(dat_raw), "\n")

if ("wave" %in% colnames(dat_raw)) {
  dat_raw <- dat_raw[order(dat_raw$id, dat_raw$wave), ]
} else {
  dat_raw <- dat_raw[order(dat_raw$id), ]
}

# SHLT recode: 1-2->1, 3->2, 4-5->3
dat_raw$SHLT_num <- as.numeric(dat_raw$SHLT)
dat_raw$SHLT_recoded <- NA_integer_
dat_raw$SHLT_recoded[dat_raw$SHLT_num %in% c(1, 2)] <- 1L
dat_raw$SHLT_recoded[dat_raw$SHLT_num == 3] <- 2L
dat_raw$SHLT_recoded[dat_raw$SHLT_num %in% c(4, 5)] <- 3L

# Drop illegal binary codes (>1) for HIBP / DIAB / LUNG
binary_code_vars <- c("HIBP", "DIAB", "LUNG")
binary_code_vars <- binary_code_vars[binary_code_vars %in% colnames(dat_raw)]
if (length(binary_code_vars) > 0) {
  keep_bin <- rep(TRUE, nrow(dat_raw))
  for (var in binary_code_vars) {
    keep_bin <- keep_bin & !is.na(dat_raw[[var]]) & dat_raw[[var]] <= 1
  }
  n_drop_bin <- sum(!keep_bin)
  if (n_drop_bin > 0) cat("去除 HIBP/DIAB/LUNG > 1 的行数:", n_drop_bin, "\n")
  dat_raw <- dat_raw[keep_bin, ]
}

need_cols <- c("id", "HSPNIT", "SHLT_recoded", covariate_vars)
missing_cols <- setdiff(need_cols, colnames(dat_raw))
if (length(missing_cols) > 0) {
  stop("缺少列: ", paste(missing_cols, collapse = ", "))
}

# Complete-case: no imputation
cc_vars <- c("HSPNIT", "SHLT_recoded", covariate_vars)
complete <- stats::complete.cases(dat_raw[, cc_vars, drop = FALSE])
n_drop_cc <- sum(!complete)
if (n_drop_cc > 0) {
  cat("完整案例分析：去除含 NA 的行数:", n_drop_cc, "\n")
}
dat <- dat_raw[complete, ]

# Standardize BMI on the analysis sample (mean 0, sd 1)
bmi_scale <- NULL
if ("BMI" %in% covariate_vars) {
  bmi_raw <- as.numeric(dat$BMI)
  bmi_mu <- mean(bmi_raw, na.rm = TRUE)
  bmi_sd <- stats::sd(bmi_raw, na.rm = TRUE)
  if (!is.finite(bmi_sd) || bmi_sd < 1e-12) stop("BMI 标准差无效，无法标准化")
  dat$BMI <- (bmi_raw - bmi_mu) / bmi_sd
  bmi_scale <- c(mean = bmi_mu, sd = bmi_sd)
  cat(sprintf("BMI 已标准化: mean = %.4f, sd = %.4f\n", bmi_mu, bmi_sd))
}

y1 <- dat$HSPNIT
y2 <- dat$SHLT_recoded
ID <- dat$id

C <- length(unique(y2[!is.na(y2)]))
cat("SHLT 合并后类别数 C =", C, "\n")
print(table(y2, useNA = "ifany"))
cat("HSPNIT 零值比例:", round(mean(y1 == 0, na.rm = TRUE) * 100, 2), "%\n")

cov_df <- dat[, covariate_vars, drop = FALSE]
X_ordinal <- cbind(1, as.matrix(cov_df))
X_zero <- cbind(1, as.matrix(cov_df))
X_count <- cbind(1, as.matrix(cov_df))
colnames(X_ordinal)[1] <- "Intercept"
colnames(X_zero)[1] <- "Intercept"
colnames(X_count)[1] <- "Intercept"

id <- as.numeric(factor(ID))
ord_rows <- order(id)
id <- id[ord_rows]
y1 <- y1[ord_rows]
y2 <- y2[ord_rows]
X_ordinal <- X_ordinal[ord_rows, , drop = FALSE]
X_zero <- X_zero[ord_rows, , drop = FALSE]
X_count <- X_count[ord_rows, , drop = FALSE]

n <- length(unique(id))
N <- length(id)
nis <- as.numeric(table(id))
id_index <- split(seq_len(N), id)

p_ordinal <- ncol(X_ordinal)
p_zero <- ncol(X_zero)
p_count <- ncol(X_count)

cat("\n========== 协变量（论文对齐）==========\n")
cat("X_ordinal / X_zero / X_count 维:", p_ordinal,
    "  变量:", paste(colnames(X_ordinal), collapse = ", "), "\n")
cat("个体数 n =", n, "  总观测 N =", N,
    "  平均访视 =", round(mean(nis), 2), "\n")

dat_fit <- list(
  id = id, N = N, n = n, nis = nis, C = C,
  X_ordinal = X_ordinal, X_zero = X_zero, X_count = X_count,
  y1 = y1, y2 = y2
)

alpha_names <- paste0("alpha_zero_", colnames(X_zero))
beta_names <- paste0("beta_count_", colnames(X_count))
gamma_names <- paste0("gamma_ord_", colnames(X_ordinal))
delta_names <- paste0("delta", seq_len(C - 1))
rho_names <- c("rho12", "rho13", "rho23")
sigma_diag_names <- c("Sigma11", "Sigma22", "Sigma33")

############################################################
# 2. Single-chain samplers
############################################################

run_joint_chain <- function(chain_id,
                            chain_length = NULL, burn = NULL, thin = NULL,
                            initial_seed = 2025, G = NULL,
                            re_prior_arg = NULL) {
  if (is.null(chain_length)) chain_length <- get("chain_length", envir = .GlobalEnv)
  if (is.null(burn)) burn <- get("burn", envir = .GlobalEnv)
  if (is.null(thin)) thin <- get("thin", envir = .GlobalEnv)
  if (is.null(G)) G <- get("G_mix", envir = .GlobalEnv)
  re_prior_use <- if (is.null(re_prior_arg)) {
    get("re_prior", envir = .GlobalEnv)
  } else re_prior_arg
  re_prior_use <- tolower(as.character(re_prior_use)[1])
  use_cdpmm_local <- identical(re_prior_use, "cdpmm")
  cat(paste0("\n=== 联合模型 链 ", chain_id,
             " (seed=", initial_seed + chain_id,
             ", re_prior=", re_prior_use, ", G=", G, ") ===\n"))
  set.seed(initial_seed + chain_id)

  q_re <- 3L
  alpha0 <- rep(0, p_zero)
  beta0 <- rep(0, p_count)
  gamma0 <- rep(0, p_ordinal)
  T0a <- diag(0.001, p_zero)
  T0b <- diag(0.001, p_count)
  T0g <- diag(0.001, p_ordinal)
  shape_r <- 0.01
  rate_r <- 0.01

  alpha <- rnorm(p_zero, 0, 0.5)
  beta <- rnorm(p_count, 0, 0.5)
  gamma_ord <- rnorm(p_ordinal, 0, 0.5)
  delta <- seq(0, 3, length.out = C - 1)
  b <- matrix(rnorm(n * q_re, 0, 0.5), n, q_re)
  r <- 1.0
  l <- rep(0, N)
  omega_ord <- rep(1, N)
  u_est <- as.numeric(y1 > 0)

  S0_iw <- diag(iw_S0_scale, q_re)
  nu0_iw <- max(iw_nu0, q_re + 2)

  G <- max(2L, as.integer(G))
  zeta0 <- rep(0, q_re)
  Psi0_inv <- diag(1 / zeta0_sd2, q_re)

  if (use_cdpmm_local) {
    tau <- 1
    nu <- c(rbeta(G - 1L, 1, tau), 1)
    pi_w <- stick_break_weights(nu)
    zeta <- as.numeric(rmvnorm(1, zeta0, diag(zeta0_sd2, q_re)))
    mu_star <- matrix(0, G, q_re)
    Omega_list <- vector("list", G)
    Omega_inv_list <- vector("list", G)
    for (g in seq_len(G)) {
      Omega_list[[g]] <- safe_riwish(nu0_iw, S0_iw)
      Omega_inv_list[[g]] <- safe_solve(Omega_list[[g]])
      mu_star[g, ] <- as.numeric(rmvnorm(1, zeta, Omega_list[[g]] / kappa0))
    }
    mu_bar <- as.numeric(crossprod(pi_w, mu_star))
    mu <- sweep(mu_star, 2, mu_bar, "-")
    L <- sample.int(G, n, replace = TRUE, prob = pi_w)
    Sigma <- cdpmm_implied_Sigma(pi_w, mu, Omega_list)
    Sigma_inv <- NULL
  } else {
    tau <- NA_real_
    nu <- pi_w <- mu_star <- mu <- Omega_list <- Omega_inv_list <- L <- NULL
    n_g <- 1L
    Sigma <- safe_riwish(nu0_iw, S0_iw)
    Sigma_inv <- safe_solve(Sigma)
  }

  b1 <- b[, 1]
  b2 <- b[, 2]
  b3 <- b[, 3]

  save_every <- floor((chain_length - burn) / thin)
  Alpha_store <- matrix(NA, save_every, p_zero)
  Beta_store <- matrix(NA, save_every, p_count)
  Gamma_store <- matrix(NA, save_every, p_ordinal)
  Delta_store <- matrix(NA, save_every, C - 1)
  Sigma_store <- matrix(NA, save_every, 9)
  Rho_store <- matrix(NA, save_every, 3)
  R_store <- rep(NA, save_every)
  Tau_store <- rep(NA, save_every)
  Nclust_store <- rep(NA, save_every)
  loglik_y1 <- matrix(NA, save_every, N)
  loglik_y2 <- matrix(NA, save_every, N)
  # Subject-level RE b_i = (b1,b2,b3): running mean/sd; optional thinned draws
  save_b_draws <- identical(Sys.getenv("SAVE_B_DRAWS", unset = "0"), "1")
  B_sum <- matrix(0, n, q_re)
  B_ssq <- matrix(0, n, q_re)
  B_store <- if (save_b_draws) array(NA_real_, dim = c(save_every, n, q_re)) else NULL

  start_time <- proc.time()
  for (iter in seq_len(chain_length)) {
    # ---- 1. Ordinal block (y2): gamma + b3 ----
    eta_ord <- as.numeric(X_ordinal %*% gamma_ord + rep(b3, times = nis))
    l <- update_latent_ordinal(y2, eta_ord, delta, omega_ord, C)
    psi_pg <- pmin(pmax(l - eta_ord, -50), 50)
    omega_ord <- rpg(N, 2, psi_pg)

    V_gamma <- solve(crossprod(sqrt(omega_ord) * X_ordinal) + T0g)
    m_gamma <- V_gamma %*% (T0g %*% gamma0 +
                              crossprod(X_ordinal, omega_ord * (l - rep(b3, times = nis))))
    gamma_ord <- as.numeric(rmvnorm(1, m_gamma, V_gamma))
    delta <- update_thresholds(y2, l, delta, C, delta_min, delta_max, fix_first = TRUE)

    # ---- 2. ZINB block (y1): alpha + b1, beta + b2 ----
    eta_zero <- as.numeric(X_zero %*% alpha + rep(b1, times = nis))
    eta_zero <- pmin(pmax(eta_zero, -10), 10)
    pi_at_risk <- inv_logit(eta_zero)
    pi_at_risk <- pmin(pmax(pi_at_risk, 1e-10), 1 - 1e-10)

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
      r <- rgamma(1, shape = shape_r + sum(k_crp), rate = max(rate_r - log_term, 1e-10))
    }

    # ---- 3. Random effects b = (b1, b2, b3) ----
    for (j in seq_len(n)) {
      idx <- id_index[[as.character(j)]]
      if (is.null(idx)) idx <- id_index[[j]]
      if (use_cdpmm_local) {
        g <- L[j]
        prior_prec <- Omega_inv_list[[g]]
        prior_mean <- mu[g, ]
      } else {
        prior_prec <- Sigma_inv
        prior_mean <- rep(0, q_re)
      }

      Xo_j <- X_ordinal[idx, , drop = FALSE]
      o_ord <- omega_ord[idx]
      d_ord <- sum(o_ord * (l[idx] - as.numeric(Xo_j %*% gamma_ord)))
      w_ord_sum <- sum(o_ord)

      Xz_j <- X_zero[idx, , drop = FALSE]
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

    # ---- 4. CDPMM or Gaussian Sigma update ----
    if (use_cdpmm_local) {
      log_pi <- log(pmax(pi_w, 1e-300))
      log_dens <- matrix(log_pi, n, G, byrow = TRUE)
      for (g in seq_len(G)) {
        Rchol <- tryCatch(chol(Omega_list[[g]]), error = function(e) NULL)
        if (is.null(Rchol)) {
          Omega_list[[g]] <- Omega_list[[g]] + diag(1e-4, q_re)
          Rchol <- chol(Omega_list[[g]])
          Omega_inv_list[[g]] <- safe_solve(Omega_list[[g]])
        }
        resid <- sweep(b, 2, mu[g, ], "-")
        z_std <- t(backsolve(Rchol, t(resid), transpose = TRUE))
        log_dens[, g] <- log_dens[, g] -
          sum(log(diag(Rchol))) - 0.5 * q_re * log(2 * base::pi) -
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
      tau <- rgamma(1, shape = tau_a1 + (G - 1), rate = max(tau_a2 - sum_log, 1e-8))

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
      if (any(!is.finite(eig)) || min(eig) < 1e-8) Sigma <- Sigma + diag(1e-4, q_re)
    } else {
      Scat <- crossprod(b)
      Sigma <- safe_riwish(nu0_iw + n, S0_iw + Scat)
      eig <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
      if (any(!is.finite(eig)) || min(eig) < 1e-8) Sigma <- Sigma + diag(1e-4, q_re)
      Sigma_inv <- safe_solve(Sigma)
      tau <- NA_real_
      n_g <- 1L
    }

    # Location PE: center b, shift intercepts (and CDPMM atoms if used)
    b_bar <- colMeans(b)
    if (any(abs(b_bar) > 0)) {
      b <- sweep(b, 2, b_bar, "-")
      alpha[1] <- alpha[1] + b_bar[1]
      beta[1] <- beta[1] + b_bar[2]
      gamma_ord[1] <- gamma_ord[1] + b_bar[3]
      b1 <- b[, 1]
      b2 <- b[, 2]
      b3 <- b[, 3]
      if (use_cdpmm_local) {
        mu_star <- sweep(mu_star, 2, b_bar, "-")
        zeta <- zeta - b_bar
        mu_bar <- as.numeric(crossprod(pi_w, mu_star))
        mu <- sweep(mu_star, 2, mu_bar, "-")
        Sigma <- cdpmm_implied_Sigma(pi_w, mu, Omega_list)
        eig <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
        if (any(!is.finite(eig)) || min(eig) < 1e-8) Sigma <- Sigma + diag(1e-4, q_re)
      } else {
        Scat <- crossprod(b)
        Sigma <- safe_riwish(nu0_iw + n, S0_iw + Scat)
        eig <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
        if (any(!is.finite(eig)) || min(eig) < 1e-8) Sigma <- Sigma + diag(1e-4, q_re)
        Sigma_inv <- safe_solve(Sigma)
      }
    }

    if (iter > burn && ((iter - burn) %% thin == 0)) {
      s <- (iter - burn) / thin
      Alpha_store[s, ] <- alpha
      Beta_store[s, ] <- beta
      Gamma_store[s, ] <- gamma_ord
      Delta_store[s, ] <- delta
      Sigma_store[s, ] <- c(Sigma)
      Rho_store[s, ] <- c(
        Sigma[1, 2] / sqrt(Sigma[1, 1] * Sigma[2, 2]),
        Sigma[1, 3] / sqrt(Sigma[1, 1] * Sigma[3, 3]),
        Sigma[2, 3] / sqrt(Sigma[2, 2] * Sigma[3, 3])
      )
      R_store[s] <- r
      Tau_store[s] <- if (use_cdpmm_local) tau else NA_real_
      Nclust_store[s] <- if (use_cdpmm_local) sum(n_g > 0) else 1L
      B_sum <- B_sum + b
      B_ssq <- B_ssq + b * b
      if (!is.null(B_store)) B_store[s, , ] <- b

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
    }

    if (iter %% 500 == 0) {
      cat(paste0("  链 ", chain_id, " iter ", iter, "/", chain_length,
                 " r=", round(r, 3), "\n"))
    }
  }

  loglik_total <- loglik_y1 + loglik_y2
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

  dic_res <- compute_dic_from_loglik(loglik_total, theta_bar_loglik_fun)
  waic_res <- tryCatch(compute_waic(loglik_total), error = function(e) list(value = NA_real_))
  loo_res <- tryCatch(compute_looic(loglik_total), error = function(e) list(value = NA_real_))

  elapsed_min <- (proc.time() - start_time)[3] / 60
  cat(sprintf("  链 %d 完成 (%.1f min) DIC=%.1f WAIC=%.1f LOOIC=%.1f\n",
              chain_id, elapsed_min, dic_res$DIC, waic_res$value, loo_res$value))

  B_mean <- B_sum / save_every
  B_var <- pmax(B_ssq / save_every - B_mean * B_mean, 0)
  B_sd <- sqrt(B_var)
  colnames(B_mean) <- colnames(B_sd) <- c("b1_zero", "b2_count", "b3_ord")

  out_chain <- list(
    Alpha = Alpha_store, Beta = Beta_store, Gamma = Gamma_store,
    Delta = Delta_store, Sigma = Sigma_store, Rho = Rho_store,
    R = R_store, Tau = Tau_store, Nclust = Nclust_store,
    B_mean = B_mean, B_sd = B_sd,
    loglik = loglik_total,
    dic = dic_res$DIC, waic = waic_res$value, looic = loo_res$value,
    chain_id = chain_id, elapsed_min = elapsed_min, re_prior = re_prior_use
  )
  if (!is.null(B_store)) out_chain$B <- B_store
  out_chain
}

run_ordinal_chain <- function(chain_id,
                              chain_length = NULL, burn = NULL, thin = NULL,
                              initial_seed = 3025, G = NULL) {
  if (is.null(chain_length)) chain_length <- get("chain_length", envir = .GlobalEnv)
  if (is.null(burn)) burn <- get("burn", envir = .GlobalEnv)
  if (is.null(thin)) thin <- get("thin", envir = .GlobalEnv)
  if (is.null(G)) G <- get("G_mix", envir = .GlobalEnv)
  cat(paste0("\n=== 有序单独模型 链 ", chain_id, " ===\n"))
  set.seed(initial_seed + chain_id)

  X <- X_ordinal
  y <- y2
  p <- p_ordinal
  gamma0 <- rep(0, p)
  T0 <- diag(0.001, p)
  gamma_ord <- rnorm(p, 0, 0.5)
  delta <- seq(0, 3, length.out = C - 1)
  b <- rnorm(n, 0, 0.5)
  cdp <- init_cdpmm_1d(n, G = G)
  l <- rep(0, N)
  omega <- rep(1, N)

  save_every <- floor((chain_length - burn) / thin)
  Gamma_store <- matrix(NA, save_every, p)
  Delta_store <- matrix(NA, save_every, C - 1)
  Sigma_store <- rep(NA, save_every)
  loglik <- matrix(NA, save_every, N)

  for (iter in seq_len(chain_length)) {
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

    # Location PE (ordinal separate): center b, shift gamma intercept
    b_bar <- mean(b)
    if (is.finite(b_bar) && abs(b_bar) > 0) {
      b <- b - b_bar
      gamma_ord[1] <- gamma_ord[1] + b_bar
      cdp <- shift_cdpmm_1d_by(cdp, b_bar)
    }

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
  waic_res <- tryCatch(compute_waic(loglik), error = function(e) list(value = NA_real_))
  loo_res <- tryCatch(compute_looic(loglik), error = function(e) list(value = NA_real_))

  list(
    Gamma = Gamma_store, Delta = Delta_store, Sigma2 = Sigma_store,
    loglik = loglik, dic = dic_res$DIC, waic = waic_res$value, looic = loo_res$value,
    chain_id = chain_id
  )
}

run_zinb_chain <- function(chain_id,
                           chain_length = NULL, burn = NULL, thin = NULL,
                           initial_seed = 4025, G = NULL) {
  if (is.null(chain_length)) chain_length <- get("chain_length", envir = .GlobalEnv)
  if (is.null(burn)) burn <- get("burn", envir = .GlobalEnv)
  if (is.null(thin)) thin <- get("thin", envir = .GlobalEnv)
  if (is.null(G)) G <- get("G_mix", envir = .GlobalEnv)
  cat(paste0("\n=== ZINB 单独模型 链 ", chain_id, " ===\n"))
  set.seed(initial_seed + chain_id)

  Xz <- X_zero
  Xc <- X_count
  pz <- p_zero
  pc <- p_count
  alpha0 <- rep(0, pz)
  beta0 <- rep(0, pc)
  T0a <- diag(0.001, pz)
  T0b <- diag(0.001, pc)

  alpha <- rnorm(pz, 0, 0.5)
  beta <- rnorm(pc, 0, 0.5)
  r <- 1.0
  b1 <- rnorm(n, 0, 0.5)
  b2 <- rnorm(n, 0, 0.5)
  cdp1 <- init_cdpmm_1d(n, G = G)
  cdp2 <- init_cdpmm_1d(n, G = G)
  u_est <- as.numeric(y1 > 0)

  save_every <- floor((chain_length - burn) / thin)
  Alpha_store <- matrix(NA, save_every, pz)
  Beta_store <- matrix(NA, save_every, pc)
  R_store <- rep(NA, save_every)
  Sigma2_b1_store <- rep(NA, save_every)
  Sigma2_b2_store <- rep(NA, save_every)
  loglik <- matrix(NA, save_every, N)

  for (iter in seq_len(chain_length)) {
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

    for (i in seq_len(N)) {
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
      for (j in seq_len(n_pos)) {
        if (y1_pos[j] > 0) {
          probs <- r / (r + 0:(y1_pos[j] - 1))
          probs <- pmin(pmax(probs, 1e-10), 1 - 1e-10)
          k_crp[j] <- sum(rbinom(y1_pos[j], 1, probs))
        }
      }
      eta_count_current <- as.numeric(Xc_pos %*% beta + b2_pos)
      psi_current <- inv_logit(eta_count_current)
      psi_current <- pmin(pmax(psi_current, 1e-10), 1 - 1e-10)
      log_term <- sum(log(1 - psi_current))
      if (!is.finite(log_term)) log_term <- 0
      r <- rgamma(1, shape = 0.01 + sum(k_crp), rate = max(0.01 - log_term, 1e-10))
    }

    for (j in seq_len(n)) {
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

    for (j in seq_len(n)) {
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

    # Location PE (ZINB separate): center b1/b2, shift alpha/beta intercepts
    b1_bar <- mean(b1)
    if (is.finite(b1_bar) && abs(b1_bar) > 0) {
      b1 <- b1 - b1_bar
      alpha[1] <- alpha[1] + b1_bar
      cdp1 <- shift_cdpmm_1d_by(cdp1, b1_bar)
    }
    b2_bar <- mean(b2)
    if (is.finite(b2_bar) && abs(b2_bar) > 0) {
      b2 <- b2 - b2_bar
      beta[1] <- beta[1] + b2_bar
      cdp2 <- shift_cdpmm_1d_by(cdp2, b2_bar)
    }

    if (iter > burn && ((iter - burn) %% thin == 0)) {
      s <- (iter - burn) / thin
      Alpha_store[s, ] <- alpha
      Beta_store[s, ] <- beta
      R_store[s] <- r
      Sigma2_b1_store[s] <- cdp1$sigma2
      Sigma2_b2_store[s] <- cdp2$sigma2
      eta_zero_s <- as.numeric(Xz %*% alpha + rep(b1, times = nis))
      pi_sv <- inv_logit(eta_zero_s)
      eta_count_s <- as.numeric(Xc %*% beta + rep(b2, times = nis))
      phi_sv <- inv_logit(eta_count_s)
      mu_sv <- r * phi_sv / (1 - phi_sv)
      mu_sv <- pmax(mu_sv, 1e-10)
      loglik[s, ] <- ifelse(
        y1 == 0,
        log((1 - pi_sv) + pi_sv * dnbinom(y1, size = r, mu = mu_sv)),
        log(pi_sv) + dnbinom(y1, size = r, mu = mu_sv, log = TRUE)
      )
    }
  }

  alpha_mean <- colMeans(Alpha_store, na.rm = TRUE)
  beta_mean <- colMeans(Beta_store, na.rm = TRUE)
  r_mean <- mean(R_store, na.rm = TRUE)
  theta_bar_loglik_fun <- function() {
    eta_zero <- as.numeric(Xz %*% alpha_mean)
    pi_sv <- inv_logit(eta_zero)
    eta_count <- as.numeric(Xc %*% beta_mean)
    phi_sv <- inv_logit(eta_count)
    mu_sv <- r_mean * phi_sv / (1 - phi_sv)
    mu_sv <- pmax(mu_sv, 1e-10)
    ll <- ifelse(
      y1 == 0,
      log((1 - pi_sv) + pi_sv * dnbinom(y1, size = r_mean, mu = mu_sv)),
      log(pi_sv) + dnbinom(y1, size = r_mean, mu = mu_sv, log = TRUE)
    )
    -2 * sum(ll)
  }
  dic_res <- compute_dic_from_loglik(loglik, theta_bar_loglik_fun)
  waic_res <- tryCatch(compute_waic(loglik), error = function(e) list(value = NA_real_))
  loo_res <- tryCatch(compute_looic(loglik), error = function(e) list(value = NA_real_))

  list(
    Alpha = Alpha_store, Beta = Beta_store, R = R_store,
    Sigma2_b1 = Sigma2_b1_store, Sigma2_b2 = Sigma2_b2_store,
    loglik = loglik, dic = dic_res$DIC, waic = waic_res$value, looic = loo_res$value,
    chain_id = chain_id
  )
}

############################################################
# 3. Parallel multi-chain driver
############################################################

run_chains_parallel <- function(chain_fun, n_chains, export_names, label = "model") {
  n_cores_env <- suppressWarnings(as.integer(Sys.getenv("N_CORES", unset = "")))
  n_cores_env <- n_cores_env[is.finite(n_cores_env) & n_cores_env >= 1L]
  n_cores <- if (length(n_cores_env)) {
    n_cores_env[[1]]
  } else {
    max(1L, parallel::detectCores(logical = TRUE))
  }
  n_cores <- min(n_cores, n_chains)
  cat(sprintf("\n=== %s: %d 链, %d worker ===\n", label, n_chains, n_cores))

  task <- function(chain_id) chain_fun(chain_id)

  if (n_cores <= 1L || n_chains <= 1L) {
    return(lapply(seq_len(n_chains), task))
  }

  # Linux/HPC: fork (mclapply). Windows: PSOCK.
  if (.Platform$OS.type == "unix") {
    return(parallel::mclapply(seq_len(n_chains), task, mc.cores = n_cores))
  }

  cl <- parallel::makeCluster(n_cores, type = "PSOCK")
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterEvalQ(cl, {
    library(BayesLogit)
    library(mvtnorm)
    library(MCMCpack)
    library(truncnorm)
    library(loo)
  })
  parallel::clusterExport(cl, c("chain_fun", "task", export_names), envir = environment())
  parallel::parLapply(cl, seq_len(n_chains), task)
}

joint_export <- c(
  "run_joint_chain",
  "dnbinom_zero", "log_sum_exp", "inv_logit", "ordinal_loglik_vec",
  "update_latent_ordinal", "enforce_delta_order", "update_thresholds",
  "compute_dic_from_loglik", "compute_waic", "compute_looic",
  "stick_break_weights", "cdpmm_implied_Sigma", "safe_solve", "safe_riwish",
  "shift_cdpmm_1d_by",
  "y1", "y2", "id", "n", "N", "nis", "id_index", "C",
  "X_ordinal", "X_zero", "X_count", "p_ordinal", "p_zero", "p_count",
  "delta_min", "delta_max",
  "G_mix", "tau_a1", "tau_a2", "zeta0_sd2", "kappa0", "iw_nu0", "iw_S0_scale",
  "chain_length", "burn", "thin", "re_prior"
)

ordinal_export <- c(
  "run_ordinal_chain",
  "update_latent_ordinal", "enforce_delta_order", "update_thresholds",
  "ordinal_loglik_vec", "compute_dic_from_loglik", "compute_waic", "compute_looic",
  "init_cdpmm_1d", "update_cdpmm_1d", "cdpmm_implied_var_1d", "shift_cdpmm_1d_by",
  "stick_break_weights", "safe_riwish", "safe_solve",
  "y2", "id", "n", "N", "nis", "id_index", "C",
  "X_ordinal", "p_ordinal", "delta_min", "delta_max",
  "G_mix", "tau_a1", "tau_a2", "zeta0_sd2", "kappa0", "iw_nu0", "iw_S0_scale",
  "chain_length", "burn", "thin"
)

zinb_export <- c(
  "run_zinb_chain",
  "dnbinom_zero", "log_sum_exp", "inv_logit",
  "compute_dic_from_loglik", "compute_waic", "compute_looic",
  "init_cdpmm_1d", "update_cdpmm_1d", "cdpmm_implied_var_1d", "shift_cdpmm_1d_by",
  "stick_break_weights", "safe_riwish", "safe_solve",
  "y1", "id", "n", "N", "nis", "id_index",
  "X_zero", "X_count", "p_zero", "p_count",
  "G_mix", "tau_a1", "tau_a2", "zeta0_sd2", "kappa0", "iw_nu0", "iw_S0_scale",
  "chain_length", "burn", "thin"
)

cat("\n=== 联合模型 MCMC ===\n")
chains_results <- run_chains_parallel(
  run_joint_chain, n_chains, joint_export, label = paste0("Joint (", re_prior, ")")
)
saveRDS(
  list(chains_results = chains_results, settings = list(
    n_chains = n_chains, chain_length = chain_length, burn = burn, thin = thin,
    re_prior = re_prior, program_tag = program_tag, stage = "joint_only"
  )),
  file.path(OUT_DIR, paste0("checkpoint_joint_", program_tag, ".rds"))
)
cat("Saved joint checkpoint\n")

ordinal_chains <- NULL
zinb_chains <- NULL
if (isTRUE(as.integer(fit_separate) == 1L)) {
  ordinal_chains <- run_chains_parallel(
    run_ordinal_chain, n_chains, ordinal_export, label = "Ordinal-only"
  )
  zinb_chains <- run_chains_parallel(
    run_zinb_chain, n_chains, zinb_export, label = "ZINB-only"
  )
}

############################################################
# 4. Combine chains & summaries
############################################################
cat("\n=== 合并联合模型链 ===\n")

combine_field <- function(lst, field) {
  do.call(rbind, lapply(lst, `[[`, field))
}

combined_Alpha <- combine_field(chains_results, "Alpha")
combined_Beta <- combine_field(chains_results, "Beta")
combined_Gamma <- combine_field(chains_results, "Gamma")
combined_Delta <- combine_field(chains_results, "Delta")
combined_Sigma <- combine_field(chains_results, "Sigma")
combined_Rho <- combine_field(chains_results, "Rho")
combined_R <- unlist(lapply(chains_results, `[[`, "R"))
combined_Tau <- unlist(lapply(chains_results, `[[`, "Tau"))
combined_Nclust <- unlist(lapply(chains_results, `[[`, "Nclust"))

cat("\n1. 零膨胀固定效应 alpha:\n")
alpha_results <- summarize_posterior(combined_Alpha, alpha_names)
print(alpha_results, row.names = FALSE)

cat("\n2. 计数固定效应 beta:\n")
beta_results <- summarize_posterior(combined_Beta, beta_names)
print(beta_results, row.names = FALSE)

cat("\n3. 有序固定效应 gamma:\n")
gamma_results <- summarize_posterior(combined_Gamma, gamma_names)
print(gamma_results, row.names = FALSE)

cat("\n4. 阈值 delta:\n")
delta_results <- summarize_posterior(combined_Delta, delta_names)
print(delta_results, row.names = FALSE)

cat("\n5. 离散参数 r:\n")
r_results <- summarize_posterior(combined_R, "r")
print(r_results, row.names = FALSE)

sigma_diag <- cbind(combined_Sigma[, 1], combined_Sigma[, 5], combined_Sigma[, 9])
cat("\n6. 随机效应 Sigma 对角 (b1,b2,b3):\n")
sigma_results <- summarize_posterior(sigma_diag, sigma_diag_names)
print(sigma_results, row.names = FALSE)

cat("\n7. 相关系数 rho:\n")
rho_results <- summarize_posterior(combined_Rho, rho_names)
print(rho_results, row.names = FALSE)

if (use_cdpmm) {
  cat("\n8. CDPMM tau / 占用簇数:\n")
  print(summarize_posterior(combined_Tau, "tau"), row.names = FALSE)
  print(summarize_posterior(combined_Nclust, "nclust"), row.names = FALSE)
}

############################################################
# 5. Convergence diagnostics (R-hat / Geweke)
############################################################
cat("\n=== 收敛诊断 (联合模型) ===\n")

rhat_one <- function(chain_list) {
  mcmc_list <- coda::as.mcmc.list(lapply(chain_list, coda::as.mcmc))
  as.numeric(coda::gelman.diag(mcmc_list, autoburnin = FALSE, multivariate = FALSE)$psrf[1, 1])
}

pick_col <- function(res_list, field, j = 1) {
  lapply(res_list, function(x) {
    m <- x[[field]]
    if (is.matrix(m)) m[, j] else m
  })
}

diag_rows <- list()
add_diag <- function(nm, rhat_val, geweke_val) {
  diag_rows[[length(diag_rows) + 1L]] <<- data.frame(
    Parameter = nm,
    Rhat = if (is.finite(rhat_val)) round(rhat_val, 3) else NA_real_,
    Geweke_z = round(geweke_val, 3)
  )
}

for (j in seq_len(p_zero)) {
  add_diag(
    alpha_names[j],
    if (n_chains >= 2) rhat_one(pick_col(chains_results, "Alpha", j)) else NA_real_,
    as.numeric(coda::geweke.diag(coda::as.mcmc(chains_results[[1]]$Alpha[, j]))$z)
  )
}
for (j in seq_len(p_count)) {
  add_diag(
    beta_names[j],
    if (n_chains >= 2) rhat_one(pick_col(chains_results, "Beta", j)) else NA_real_,
    as.numeric(coda::geweke.diag(coda::as.mcmc(chains_results[[1]]$Beta[, j]))$z)
  )
}
for (j in seq_len(p_ordinal)) {
  add_diag(
    gamma_names[j],
    if (n_chains >= 2) rhat_one(pick_col(chains_results, "Gamma", j)) else NA_real_,
    as.numeric(coda::geweke.diag(coda::as.mcmc(chains_results[[1]]$Gamma[, j]))$z)
  )
}
add_diag(
  "r",
  if (n_chains >= 2) rhat_one(pick_col(chains_results, "R", 1)) else NA_real_,
  as.numeric(coda::geweke.diag(coda::as.mcmc(chains_results[[1]]$R))$z)
)
for (j in 1:3) {
  add_diag(
    rho_names[j],
    if (n_chains >= 2) rhat_one(pick_col(chains_results, "Rho", j)) else NA_real_,
    as.numeric(coda::geweke.diag(coda::as.mcmc(chains_results[[1]]$Rho[, j]))$z)
  )
}
diag_table <- do.call(rbind, diag_rows)
print(diag_table, row.names = FALSE)

############################################################
# 6. Separate models: variance ratios & model comparison
############################################################
variance_ratio_table <- NULL
model_compare <- NULL

if (!is.null(ordinal_chains) && !is.null(zinb_chains)) {
  ord_gamma <- combine_field(ordinal_chains, "Gamma")
  ord_sigma2 <- unlist(lapply(ordinal_chains, function(x) x$Sigma2))
  zinb_alpha <- combine_field(zinb_chains, "Alpha")
  zinb_beta <- combine_field(zinb_chains, "Beta")
  zinb_r <- unlist(lapply(zinb_chains, function(x) x$R))
  zinb_s2_b1 <- unlist(lapply(zinb_chains, function(x) x$Sigma2_b1))
  zinb_s2_b2 <- unlist(lapply(zinb_chains, function(x) x$Sigma2_b2))

  vr_alpha <- safe_col_var(zinb_alpha) / safe_col_var(combined_Alpha)
  vr_beta <- safe_col_var(zinb_beta) / safe_col_var(combined_Beta)
  vr_gamma <- safe_col_var(ord_gamma) / safe_col_var(combined_Gamma)
  vr_r <- stats::var(zinb_r, na.rm = TRUE) / stats::var(combined_R, na.rm = TRUE)
  vr_sigma2_b3 <- stats::var(ord_sigma2, na.rm = TRUE) / stats::var(combined_Sigma[, 9], na.rm = TRUE)
  vr_sigma2_b1 <- stats::var(zinb_s2_b1, na.rm = TRUE) / stats::var(combined_Sigma[, 1], na.rm = TRUE)
  vr_sigma2_b2 <- stats::var(zinb_s2_b2, na.rm = TRUE) / stats::var(combined_Sigma[, 5], na.rm = TRUE)

  variance_ratio_table <- data.frame(
    Parameter = c(
      paste0("alpha", seq_along(vr_alpha)),
      paste0("beta", seq_along(vr_beta)),
      paste0("gamma", seq_along(vr_gamma)),
      "r", "sigma2_b3", "sigma2_b1", "sigma2_b2"
    ),
    Variance_Ratio = c(vr_alpha, vr_beta, vr_gamma, vr_r,
                       vr_sigma2_b3, vr_sigma2_b1, vr_sigma2_b2)
  )
  variance_ratio_table$Variance_Ratio <- round(variance_ratio_table$Variance_Ratio, 4)

  cat("\n================ 方差比 (单独/联合 后验方差) ================\n")
  print(variance_ratio_table, row.names = FALSE)

  mean_metric <- function(lst, field) mean(vapply(lst, `[[`, numeric(1), field), na.rm = TRUE)
  model_compare <- data.frame(
    Model = c("Joint", "Ordinal-only", "ZINB-only"),
    DIC_Mean = c(
      mean_metric(chains_results, "dic"),
      mean_metric(ordinal_chains, "dic"),
      mean_metric(zinb_chains, "dic")
    ),
    WAIC_Mean = c(
      mean_metric(chains_results, "waic"),
      mean_metric(ordinal_chains, "waic"),
      mean_metric(zinb_chains, "waic")
    ),
    LOOIC_Mean = c(
      mean_metric(chains_results, "looic"),
      mean_metric(ordinal_chains, "looic"),
      mean_metric(zinb_chains, "looic")
    )
  )
  model_compare[, -1] <- round(model_compare[, -1], 2)

  cat("\n================ 模型比较 (DIC / WAIC / LOOIC 链均值) ================\n")
  print(model_compare, row.names = FALSE)
}

############################################################
# 7. Save results (prog1: CDPMM joint vs separate)
############################################################
rds_tag <- paste0("MCMC_Analysis_Results_", program_tag, "_cdpmm.rds")
rds_path <- file.path(OUT_DIR, rds_tag)

out <- list(
  settings = list(
    program = "prog1_cdpmm_joint_vs_separate",
    re_prior = re_prior,
    covariates = covariate_vars,
    bmi_scale = bmi_scale,
    n = n,
    N = N,
    C = C,
    n_chains = n_chains,
    chain_length = chain_length,
    burn = burn,
    thin = thin,
    G_mix = G_mix,
    fit_separate = as.integer(fit_separate),
    OUT_DIR = OUT_DIR
  ),
  chains_results = chains_results,
  combined_samples = list(
    Alpha = combined_Alpha, Beta = combined_Beta, Gamma = combined_Gamma,
    Delta = combined_Delta, Sigma = combined_Sigma, Rho = combined_Rho,
    R = combined_R, Tau = combined_Tau, Nclust = combined_Nclust
  ),
  summary = list(
    alpha = alpha_results, beta = beta_results, gamma = gamma_results,
    delta = delta_results, r = r_results, sigma_diag = sigma_results,
    rho = rho_results, diagnostics = diag_table
  )
)
# Pooled posterior mean of subject RE b_i across chains
if (!is.null(chains_results[[1]]$B_mean)) {
  B_mean_pool <- Reduce(`+`, lapply(chains_results, `[[`, "B_mean")) / length(chains_results)
  B_sd_pool <- sqrt(Reduce(`+`, lapply(chains_results, function(ch) {
    ch$B_sd^2 + (ch$B_mean - B_mean_pool)^2
  })) / length(chains_results))
  colnames(B_mean_pool) <- colnames(B_sd_pool) <- c("b1_zero", "b2_count", "b3_ord")
  out$B_mean <- B_mean_pool
  out$B_sd <- B_sd_pool
  write.csv(
    data.frame(subject = seq_len(n), B_mean_pool),
    file.path(OUT_DIR, paste0("Posterior_B_mean_", program_tag, ".csv")),
    row.names = FALSE
  )
}
if (!is.null(ordinal_chains)) out$ordinal_chains <- ordinal_chains
if (!is.null(zinb_chains)) out$zinb_chains <- zinb_chains
if (!is.null(variance_ratio_table)) out$variance_ratio <- variance_ratio_table
if (!is.null(model_compare)) out$model_compare <- model_compare

saveRDS(out, rds_path)
write.csv(diag_table, file.path(OUT_DIR, paste0("diagnose_rhat_", program_tag, ".csv")), row.names = FALSE)
write.csv(alpha_results, file.path(OUT_DIR, paste0("Posterior_Joint_Alpha_", program_tag, ".csv")), row.names = FALSE)
write.csv(beta_results, file.path(OUT_DIR, paste0("Posterior_Joint_Beta_", program_tag, ".csv")), row.names = FALSE)
write.csv(gamma_results, file.path(OUT_DIR, paste0("Posterior_Joint_Gamma_", program_tag, ".csv")), row.names = FALSE)
write.csv(delta_results, file.path(OUT_DIR, paste0("Posterior_Joint_Delta_", program_tag, ".csv")), row.names = FALSE)
write.csv(r_results, file.path(OUT_DIR, paste0("Posterior_Joint_r_", program_tag, ".csv")), row.names = FALSE)
write.csv(sigma_results, file.path(OUT_DIR, paste0("Posterior_Joint_SigmaDiag_", program_tag, ".csv")), row.names = FALSE)
write.csv(rho_results, file.path(OUT_DIR, paste0("Posterior_Joint_Rho_", program_tag, ".csv")), row.names = FALSE)
if (!is.null(variance_ratio_table)) {
  write.csv(variance_ratio_table, file.path(OUT_DIR, paste0("variance_ratio_", program_tag, ".csv")), row.names = FALSE)
}
if (!is.null(model_compare)) {
  write.csv(model_compare, file.path(OUT_DIR, paste0("model_compare_", program_tag, ".csv")), row.names = FALSE)
  # Win indicator: joint criteria < sum of separate (like sim prog1)
  if (nrow(model_compare) >= 3L) {
    j <- model_compare[model_compare$Model == "Joint", , drop = FALSE]
    o <- model_compare[model_compare$Model == "Ordinal-only", , drop = FALSE]
    z <- model_compare[model_compare$Model == "ZINB-only", , drop = FALSE]
    if (nrow(j) && nrow(o) && nrow(z)) {
      criteria_win <- data.frame(
        Comparison = "CDPMM_joint_better_than_sum_of_separate",
        DIC_joint_lt_sep = as.integer(j$DIC_Mean < (o$DIC_Mean + z$DIC_Mean)),
        WAIC_joint_lt_sep = as.integer(j$WAIC_Mean < (o$WAIC_Mean + z$WAIC_Mean)),
        LOOIC_joint_lt_sep = as.integer(j$LOOIC_Mean < (o$LOOIC_Mean + z$LOOIC_Mean)),
        DIC_Joint = j$DIC_Mean,
        DIC_SepSum = o$DIC_Mean + z$DIC_Mean,
        WAIC_Joint = j$WAIC_Mean,
        WAIC_SepSum = o$WAIC_Mean + z$WAIC_Mean,
        LOOIC_Joint = j$LOOIC_Mean,
        LOOIC_SepSum = o$LOOIC_Mean + z$LOOIC_Mean
      )
      write.csv(criteria_win, file.path(OUT_DIR, paste0("criteria_compare_", program_tag, ".csv")),
                row.names = FALSE)
      cat("\n================ 准则比较 (联合 vs 单独之和) ================\n")
      print(criteria_win, row.names = FALSE)
    }
  }
}

cat("\n已保存 RDS:", rds_path, "\n")
cat("CSV 与诊断表已写入:", OUT_DIR, "\n")
cat("绘图: 运行 plot_diagnostics.R（可读 prog1 RDS）\n")
cat("=== Prog1 完成 (n=", n, ", N=", N, ") ===\n", sep = "")
