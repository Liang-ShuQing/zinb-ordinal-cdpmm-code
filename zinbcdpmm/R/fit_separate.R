# Separate ordinal / ZINB CDPMM models
# Source: 模拟研究/prog1_cdpmm_joint_vs_separate.R (algorithm bodies preserved).

fit_ordinal_model <- function(dat, chain = 5000, burn = 2000, thin = 5,
                              delta_min = -10, delta_max = 10,
                              G = 8L,
                              iw_nu0 = 6, iw_S0_scale = 1,
                              zeta0_sd2 = 10, kappa0 = 1,
                              tau_a1 = 2, tau_a2 = 4) {
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
  cdp <- init_cdpmm_1d(n, G = G, iw_nu0 = iw_nu0, iw_S0_scale = iw_S0_scale,
                      zeta0_sd2 = zeta0_sd2, kappa0 = kappa0,
                      tau_a1 = tau_a1, tau_a2 = tau_a2)
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
    omega <- BayesLogit::rpg(N, 2, psi_pg)

    V_g <- solve(crossprod(sqrt(omega) * X) + T0)
    m_g <- V_g %*% (T0 %*% gamma0 + crossprod(X, omega * (l - rep(b, times = nis))))
    gamma_ord <- as.numeric(mvtnorm::rmvnorm(1, m_g, V_g))

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
                           G = 8L,
                           iw_nu0 = 6, iw_S0_scale = 1,
                           zeta0_sd2 = 10, kappa0 = 1,
                           tau_a1 = 2, tau_a2 = 4) {
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
  cdp1 <- init_cdpmm_1d(n, G = G, iw_nu0 = iw_nu0, iw_S0_scale = iw_S0_scale,
                       zeta0_sd2 = zeta0_sd2, kappa0 = kappa0,
                       tau_a1 = tau_a1, tau_a2 = tau_a2)
  cdp2 <- init_cdpmm_1d(n, G = G, iw_nu0 = iw_nu0, iw_S0_scale = iw_S0_scale,
                       zeta0_sd2 = zeta0_sd2, kappa0 = kappa0,
                       tau_a1 = tau_a1, tau_a2 = tau_a2)
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
    
    omega_zero <- BayesLogit::rpg(N, 1, eta_zero)
    z_a <- (u_est - 0.5) / pmax(omega_zero, 1e-10)
    
    if (n_pos > 0) {
      y1_pos <- y1[pos_idx]
      Xc_pos <- Xc[pos_idx, , drop = FALSE]
      b2_pos <- rep(b2, times = nis)[pos_idx]
      eta_count_pos <- as.numeric(Xc_pos %*% beta + b2_pos)
      eta_count_pos <- pmin(pmax(eta_count_pos, -10), 10)
      w_count <- BayesLogit::rpg(n_pos, y1_pos + r, eta_count_pos)
      z_b <- (y1_pos - r) / (2 * pmax(w_count, 1e-10))
    }
    
    V_a <- solve(crossprod(sqrt(omega_zero) * Xz) + T0a)
    m_a <- V_a %*% (T0a %*% alpha0 + crossprod(Xz, omega_zero * (z_a - rep(b1, times = nis))))
    alpha <- as.numeric(mvtnorm::rmvnorm(1, m_a, V_a))
    
    if (n_pos > 0) {
      V_b <- solve(crossprod(sqrt(w_count) * Xc_pos) + T0b)
      m_b <- V_b %*% (T0b %*% beta0 + crossprod(Xc_pos, w_count * (z_b - b2_pos)))
      beta <- as.numeric(mvtnorm::rmvnorm(1, m_b, V_b))
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

