# Utility helpers
# Source: 模拟研究/prog1_cdpmm_joint_vs_separate.R (algorithm bodies preserved).

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
      l[i] <- truncnorm::rtruncnorm(1, a = -Inf, b = delta[1], mean = eta[i], sd = sd_i)
    } else if (k == C) {
      l[i] <- truncnorm::rtruncnorm(1, a = delta[C - 1], b = Inf, mean = eta[i], sd = sd_i)
    } else {
      l[i] <- truncnorm::rtruncnorm(1, a = delta[k - 1], b = delta[k], mean = eta[i], sd = sd_i)
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

