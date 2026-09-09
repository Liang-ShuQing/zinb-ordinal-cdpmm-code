# Joint ZINB + ordinal CDPMM Gibbs
# Source: 模拟研究/prog1_cdpmm_joint_vs_separate.R (algorithm bodies preserved).

fit_joint_model <- function(dat, chain = 5000, burn = 2000, thin = 5,
                            delta_min = -10, delta_max = 10,
                            G = 8L, init = NULL,
                            alpha0 = NULL, beta0 = NULL, gamma0 = NULL,
                            prior_var = 1000,
                            r_update = c("crt", "mh"),
                            mh_r_sd = 0.20,
                            iw_nu0 = 6, iw_S0_scale = 1,
                            zeta0_sd2 = 10, kappa0 = 1,
                            tau_a1 = 2, tau_a2 = 4) {
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
  zeta <- as.numeric(mvtnorm::rmvnorm(1, zeta0, diag(zeta0_sd2, q_re)))
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
    mu_star[g, ] <- as.numeric(mvtnorm::rmvnorm(1, zeta, Omega_list[[g]] / kappa0))
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
          b[i, ] <- as.numeric(mvtnorm::rmvnorm(1, mu[g, ], Omega_list[[g]]))
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
    omega_ord <- BayesLogit::rpg(N, 2, psi_pg)

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

    omega_zero <- BayesLogit::rpg(N, 1, eta_zero)
    z_alpha <- (u_est - 0.5) / pmax(omega_zero, 1e-10)

    if (n_pos > 0) {
      y1_pos <- y1[pos_idx]
      X_pos <- X_count[pos_idx, , drop = FALSE]
      b2_pos <- rep(b2, times = nis)[pos_idx]

      eta_count_pos <- as.numeric(X_pos %*% beta + b2_pos)
      eta_count_pos <- pmin(pmax(eta_count_pos, -10), 10)
      w_count <- BayesLogit::rpg(n_pos, y1_pos + r, eta_count_pos)
      z_beta <- (y1_pos - r) / (2 * pmax(w_count, 1e-10))
    }

    V_alpha <- solve(crossprod(sqrt(omega_zero) * X_zero) + T0a)
    m_alpha <- V_alpha %*% (T0a %*% alpha0 +
                              crossprod(X_zero, omega_zero * (z_alpha - rep(b1, times = nis))))
    alpha <- as.numeric(mvtnorm::rmvnorm(1, m_alpha, V_alpha))

    if (n_pos > 0) {
      V_beta <- solve(crossprod(sqrt(w_count) * X_pos) + T0b)
      m_beta <- V_beta %*% (T0b %*% beta0 +
                              crossprod(X_pos, w_count * (z_beta - b2_pos)))
      beta <- as.numeric(mvtnorm::rmvnorm(1, m_beta, V_beta))
    }

    V_gamma <- solve(crossprod(sqrt(omega_ord) * X_ordinal) + T0g)
    m_gamma <- V_gamma %*% (T0g %*% gamma0 +
                              crossprod(X_ordinal, omega_ord * (l - rep(b3, times = nis))))
    gamma_ord <- as.numeric(mvtnorm::rmvnorm(1, m_gamma, V_gamma))

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
      b[j, ] <- as.numeric(mvtnorm::rmvnorm(1, post_mean, post_var))
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
          mu_star[g, ] <- as.numeric(mvtnorm::rmvnorm(1, zeta, Omega_list[[g]] / kappa0))
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
          mu_star[g, ] <- as.numeric(mvtnorm::rmvnorm(1, m_n, Omega_list[[g]] / kn))
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
      zeta <- as.numeric(mvtnorm::rmvnorm(1, V_zeta %*% num_zeta, V_zeta))

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

