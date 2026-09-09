# CDPMM stick-breaking / 1D NIW updates
# Source: 模拟研究/prog1_cdpmm_joint_vs_separate.R (algorithm bodies preserved).

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

init_cdpmm_1d <- function(n, G = 8L,
                          iw_nu0 = 6, iw_S0_scale = 1,
                          zeta0_sd2 = 10, kappa0 = 1,
                          tau_a1 = 2, tau_a2 = 4) {
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
    mu_star = mu_star, mu = mu, omega = omega, L = L,
    zeta0_sd2 = zeta0_sd2, kappa0 = kappa0,
    tau_a1 = tau_a1, tau_a2 = tau_a2
  )
}

# Strict 1D NIW updates (Tang stick-breaking/centering retained)
update_cdpmm_1d <- function(st, b,
                            tau_a1 = NULL, tau_a2 = NULL,
                            zeta0_sd2 = NULL, kappa0 = NULL) {
  if (is.null(tau_a1)) tau_a1 <- if (!is.null(st$tau_a1)) st$tau_a1 else 2
  if (is.null(tau_a2)) tau_a2 <- if (!is.null(st$tau_a2)) st$tau_a2 else 4
  if (is.null(zeta0_sd2)) zeta0_sd2 <- if (!is.null(st$zeta0_sd2)) st$zeta0_sd2 else 10
  if (is.null(kappa0)) kappa0 <- if (!is.null(st$kappa0)) st$kappa0 else 1
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

