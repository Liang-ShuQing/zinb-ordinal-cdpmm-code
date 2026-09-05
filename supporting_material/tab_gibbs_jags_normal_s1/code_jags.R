############################################################
# JAGS-only Monte Carlo (scheme 1: ~1000 retained draws)
# Env: SCENARIO, N / N_SUBJ, RE_DIST, NSIM, OUT_DIR
# Seeds match main Gibbs table: generate_data(..., seed = 2026 + s)
# Retained: JAGS_CHAINS * floor(JAGS_ITER / JAGS_THIN) ≈ 1000
############################################################

.env_int <- function(key, default) {
  v <- suppressWarnings(as.integer(Sys.getenv(key, unset = "")))
  if (length(v) == 1L && is.finite(v)) v else as.integer(default)
}
.env_chr <- function(key, default) {
  v <- Sys.getenv(key, unset = "")
  if (nzchar(v)) v else as.character(default)
}

args <- commandArgs(trailingOnly = TRUE)
OUT_DIR <- if (length(args) >= 1L && nzchar(args[[1]])) args[[1]] else .env_chr("OUT_DIR", ".")
if (!nzchar(as.character(OUT_DIR)[1])) OUT_DIR <- "."
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_DIR <- normalizePath(OUT_DIR, winslash = "/", mustWork = TRUE)

n_sim <- .env_int("NSIM", 500L)
n <- .env_int("N_SUBJ", .env_int("N", 200L))
G_mix <- .env_int("G_MIX", 8L)
scenario <- .env_int("SCENARIO", 2L)
re_dist <- tolower(.env_chr("RE_DIST", "mixture"))
C <- 5L

jags_chains <- .env_int("JAGS_CHAINS", 2L)
jags_adapt <- .env_int("JAGS_ADAPT", 500L)
jags_burn <- .env_int("JAGS_BURN", 500L)
jags_iter <- .env_int("JAGS_ITER", 500L)
jags_thin <- .env_int("JAGS_THIN", 1L)

cat("METHOD=JAGS OUT_DIR=", OUT_DIR, "\n")
cat("NSIM=", n_sim, " N=", n, " G=", G_mix,
    " SCENARIO=", scenario, " RE_DIST=", re_dist, "\n", sep = "")
cat("JAGS:", jags_chains, "chains, adapt", jags_adapt,
    " burn", jags_burn, " iter", jags_iter, " thin", jags_thin, "\n")

user_lib <- path.expand("~/R_libs")
dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(user_lib, .libPaths()))

suppressPackageStartupMessages({
  library(coda)
  library(foreach)
  library(doParallel)
})

SKIP_MAIN_SIM <- TRUE
Sys.setenv(
  SCENARIO = as.character(scenario),
  RE_DIST = re_dist,
  G_MIX = as.character(G_mix),
  N = as.character(n),
  N_SIM = "1",
  FIT_SEPARATE = "0"
)
source("prog1_cdpmm_joint_vs_separate.R", local = FALSE)

n_sim <- .env_int("NSIM", 500L)
n <- .env_int("N_SUBJ", .env_int("N", 200L))
G_mix <- .env_int("G_MIX", 8L)
scenario <- .env_int("SCENARIO", 2L)
re_dist <- tolower(.env_chr("RE_DIST", "mixture"))
C <- 5L
random_nis <- TRUE
nis_range <- 6:9
nis_fixed <- 7L
mvt_df <- 4

vech3 <- function(S) {
  S <- as.matrix(S)
  as.numeric(S[lower.tri(S, diag = TRUE)])
}
rho_from_Sigma <- function(S) {
  S <- as.matrix(S)
  c(
    S[1, 2] / sqrt(S[1, 1] * S[2, 2]),
    S[1, 3] / sqrt(S[1, 1] * S[3, 3]),
    S[2, 3] / sqrt(S[2, 2] * S[3, 3])
  )
}

sigma_from_jags_atoms <- function(pi_w, mu_star, Tau_arr) {
  G <- length(pi_w)
  mu_bar <- as.numeric(crossprod(pi_w, mu_star))
  mu <- sweep(mu_star, 2, mu_bar, "-")
  Sigma <- matrix(0, 3, 3)
  for (g in seq_len(G)) {
    Om <- tryCatch(solve(Tau_arr[g, , ]), error = function(e) diag(3))
    Sigma <- Sigma + pi_w[g] * (Om + tcrossprod(mu[g, ]))
  }
  Sigma
}

param_ess_mean <- function(mat) {
  if (is.null(dim(mat)) || ncol(mat) < 1) return(NA_real_)
  ess <- apply(mat, 2, function(x) {
    tryCatch(as.numeric(coda::effectiveSize(coda::as.mcmc(x))), error = function(e) NA_real_)
  })
  mean(ess, na.rm = TRUE)
}

true_vec <- function() {
  fe <- true_fixed_effects(scenario)
  c(fe$alpha, fe$beta, fe$gamma, 2)
}
fe_names <- function() {
  fe <- true_fixed_effects(scenario)
  c(
    paste0("alpha", seq_along(fe$alpha)),
    paste0("beta", seq_along(fe$beta)),
    paste0("gamma", seq_along(fe$gamma)),
    "r"
  )
}
sigma_vech_names <- c("Sigma11", "Sigma21", "Sigma22", "Sigma31", "Sigma32", "Sigma33")
rho_names <- c("rho12", "rho13", "rho23")

calculate_stats <- function(estimates, true_value, ci_lower, ci_upper) {
  n_params <- length(true_value)
  bias <- rmse <- cp <- rep(NA_real_, n_params)
  for (j in seq_len(n_params)) {
    bias[j] <- mean(estimates[, j] - true_value[j], na.rm = TRUE)
    rmse[j] <- sqrt(mean((estimates[, j] - true_value[j])^2, na.rm = TRUE))
    cp[j] <- mean(ci_lower[, j] <= true_value[j] & true_value[j] <= ci_upper[, j], na.rm = TRUE)
  }
  list(bias = bias, rmse = rmse, cp = cp)
}

calculate_stats_varying_true <- function(estimates, true_mat, ci_lower, ci_upper) {
  n_params <- ncol(estimates)
  bias <- rmse <- cp <- rep(NA_real_, n_params)
  for (j in seq_len(n_params)) {
    bias[j] <- mean(estimates[, j] - true_mat[, j], na.rm = TRUE)
    rmse[j] <- sqrt(mean((estimates[, j] - true_mat[, j])^2, na.rm = TRUE))
    cp[j] <- mean(ci_lower[, j] <= true_mat[, j] & true_mat[, j] <= ci_upper[, j], na.rm = TRUE)
  }
  list(bias = bias, rmse = rmse, cp = cp)
}

fit_jags_cdpmm <- function(dat, jags_file, G,
                           n_chains = 2L, n_adapt = 500L,
                           n_burn = 500L, n_iter = 500L, n_thin = 1L,
                           seed = 1L) {
  suppressPackageStartupMessages({
    library(rjags)
    library(coda)
  })
  set.seed(seed)
  C_loc <- as.integer(dat$C)
  jags_data <- list(
    N = as.integer(dat$n),
    Ntot = as.integer(dat$N),
    C = C_loc,
    G = as.integer(G),
    pz = ncol(dat$X_zero),
    pc = ncol(dat$X_count),
    po = ncol(dat$X_ordinal),
    id = as.integer(dat$id),
    y1 = as.integer(dat$y1),
    y2 = as.integer(dat$y2),
    Xz = dat$X_zero,
    Xc = dat$X_count,
    Xo = dat$X_ordinal,
    R = diag(3),
    wish_df = 5,
    tau_a1 = 1,
    tau_a2 = 1
  )
  inits <- function(chain_id) {
    list(
      alpha = rnorm(ncol(dat$X_zero), 0, 0.1),
      beta = rnorm(ncol(dat$X_count), 0, 0.1),
      gamma = rnorm(ncol(dat$X_ordinal), 0, 0.1),
      r = 2,
      tau = 1,
      nu = rep(0.5, G - 1L),
      .RNG.name = "base::Wichmann-Hill",
      .RNG.seed = as.integer(seed + 100L * chain_id)
    )
  }
  init_list <- lapply(seq_len(n_chains), inits)
  mod_j <- rjags::jags.model(
    file = jags_file,
    data = jags_data,
    inits = init_list,
    n.chains = n_chains,
    n.adapt = n_adapt,
    quiet = TRUE
  )
  update(mod_j, n_burn, progress.bar = "none")
  samp <- rjags::coda.samples(
    mod_j,
    variable.names = c("alpha", "beta", "gamma", "r", "pi", "mu_star", "Tau"),
    n.iter = n_iter,
    thin = n_thin,
    progress.bar = "none"
  )
  mat <- as.matrix(samp)

  pick_cols <- function(pattern) {
    nm <- grep(pattern, colnames(mat), value = TRUE)
    nm[order(nm)]
  }
  fe_nm <- c(
    pick_cols("^alpha\\["),
    pick_cols("^beta\\["),
    pick_cols("^gamma\\["),
    pick_cols("^r$")
  )
  fe_mat <- mat[, fe_nm, drop = FALSE]
  fe_est <- colMeans(fe_mat)
  fe_lo <- apply(fe_mat, 2, quantile, probs = 0.025, na.rm = TRUE)
  fe_hi <- apply(fe_mat, 2, quantile, probs = 0.975, na.rm = TRUE)
  ess_fe <- param_ess_mean(fe_mat)

  n_draw <- nrow(mat)
  # Cap Sigma reconstruction cost while keeping enough draws for CI
  take <- if (n_draw > 1000L) {
    as.integer(round(seq(1, n_draw, length.out = 1000)))
  } else {
    seq_len(n_draw)
  }
  sig_draws <- matrix(NA_real_, length(take), 6L)
  rho_draws <- matrix(NA_real_, length(take), 3L)
  pi_nm <- pick_cols("^pi\\[")
  for (ii in seq_along(take)) {
    t <- take[[ii]]
    pi_w <- as.numeric(mat[t, pi_nm])
    pi_w <- pmax(pi_w, 0)
    if (sum(pi_w) <= 0) next
    pi_w <- pi_w / sum(pi_w)
    mu_star <- matrix(0, G, 3)
    for (g in seq_len(G)) for (k in 1:3) {
      nm <- sprintf("mu_star[%d,%d]", g, k)
      if (nm %in% colnames(mat)) mu_star[g, k] <- mat[t, nm]
    }
    Tau_by_g <- array(0, dim = c(G, 3, 3))
    for (g in seq_len(G)) for (i in 1:3) for (j in 1:3) {
      nm <- sprintf("Tau[%d,%d,%d]", g, i, j)
      if (nm %in% colnames(mat)) Tau_by_g[g, i, j] <- mat[t, nm]
    }
    Sig <- sigma_from_jags_atoms(pi_w, mu_star, Tau_by_g)
    sig_draws[ii, ] <- vech3(Sig)
    rho_draws[ii, ] <- rho_from_Sigma(Sig)
  }
  ok_row <- rowSums(is.finite(sig_draws)) == 6L
  sig_draws <- sig_draws[ok_row, , drop = FALSE]
  rho_draws <- rho_draws[ok_row, , drop = FALSE]
  if (nrow(sig_draws) < 1L) {
    sig_est <- sig_lo <- sig_hi <- rep(NA_real_, 6L)
    rho_est <- rho_lo <- rho_hi <- rep(NA_real_, 3L)
  } else {
    sig_est <- colMeans(sig_draws)
    sig_lo <- apply(sig_draws, 2, quantile, probs = 0.025, na.rm = TRUE)
    sig_hi <- apply(sig_draws, 2, quantile, probs = 0.975, na.rm = TRUE)
    rho_est <- colMeans(rho_draws)
    rho_lo <- apply(rho_draws, 2, quantile, probs = 0.025, na.rm = TRUE)
    rho_hi <- apply(rho_draws, 2, quantile, probs = 0.975, na.rm = TRUE)
  }

  list(
    fe_est = as.numeric(fe_est),
    fe_lo = as.numeric(fe_lo),
    fe_hi = as.numeric(fe_hi),
    sig_est = as.numeric(sig_est),
    sig_lo = as.numeric(sig_lo),
    sig_hi = as.numeric(sig_hi),
    rho_est = as.numeric(rho_est),
    rho_lo = as.numeric(rho_lo),
    rho_hi = as.numeric(rho_hi),
    ess = ess_fe
  )
}

n_cores <- as.integer(Sys.getenv(
  c("SLURM_NTASKS_PER_NODE", "SLURM_CPUS_PER_TASK", "N_CORES"),
  unset = NA_character_
))
n_cores <- n_cores[is.finite(n_cores) & n_cores >= 1L]
n_cores <- if (length(n_cores)) n_cores[[1]] else parallel::detectCores()
n_cores <- max(1L, as.integer(n_cores))
n_workers <- max(1L, min(as.integer(n_sim), n_cores))
cat("Parallel sims: n_workers=", n_workers, " of n_cores=", n_cores, "\n", sep = "")

run_one_sim <- function(s) {
  seed_s <- 2026L + as.integer(s)
  cat("\n===== JAGS replicate", s, "/", n_sim, " pid=", Sys.getpid(), " =====\n", sep = "")
  flush.console()

  dat <- generate_data(
    n = n, nis = nis_fixed, seed = seed_s,
    random_nis = random_nis, nis_range = nis_range,
    scenario = scenario, C = C,
    re_dist = re_dist, mvt_df = mvt_df
  )
  truth <- true_vec()
  Sigma_true <- dat$Sigma_true
  rho_true <- rho_from_Sigma(Sigma_true)
  sigma_true_vech <- vech3(Sigma_true)

  jags_file <- file.path(OUT_DIR, "joint_cdpmm.jags")
  if (!file.exists(jags_file)) jags_file <- "joint_cdpmm.jags"

  n_fe <- length(truth)
  fe_est <- fe_lo <- fe_hi <- rep(NA_real_, n_fe)
  sig_est <- sig_lo <- sig_hi <- rep(NA_real_, 6L)
  rho_est <- rho_lo <- rho_hi <- rep(NA_real_, 3L)
  aess <- NA_real_

  t0 <- proc.time()[["elapsed"]]
  jfit <- tryCatch(
    fit_jags_cdpmm(
      dat, jags_file = jags_file, G = G_mix,
      n_chains = jags_chains, n_adapt = jags_adapt,
      n_burn = jags_burn, n_iter = jags_iter, n_thin = jags_thin,
      seed = seed_s
    ),
    error = function(e) {
      cat("JAGS failed (sim=", s, "):", conditionMessage(e), "\n")
      NULL
    }
  )
  sec <- proc.time()[["elapsed"]] - t0

  if (!is.null(jfit)) {
    fe_est <- jfit$fe_est
    fe_lo <- jfit$fe_lo
    fe_hi <- jfit$fe_hi
    if (length(fe_est) < n_fe) {
      fe_est <- c(fe_est, rep(NA_real_, n_fe - length(fe_est)))
      fe_lo <- c(fe_lo, rep(NA_real_, n_fe - length(fe_lo)))
      fe_hi <- c(fe_hi, rep(NA_real_, n_fe - length(fe_hi)))
    } else if (length(fe_est) > n_fe) {
      fe_est <- fe_est[seq_len(n_fe)]
      fe_lo <- fe_lo[seq_len(n_fe)]
      fe_hi <- fe_hi[seq_len(n_fe)]
    }
    sig_est <- jfit$sig_est
    sig_lo <- jfit$sig_lo
    sig_hi <- jfit$sig_hi
    rho_est <- jfit$rho_est
    rho_lo <- jfit$rho_lo
    rho_hi <- jfit$rho_hi
    aess <- jfit$ess
  }

  list(
    sim = s,
    sec = sec,
    aess = aess,
    truth_fe = truth,
    truth_sigma = sigma_true_vech,
    truth_rho = rho_true,
    fe_est = fe_est, fe_lo = fe_lo, fe_hi = fe_hi,
    sig_est = sig_est, sig_lo = sig_lo, sig_hi = sig_hi,
    rho_est = rho_est, rho_lo = rho_lo, rho_hi = rho_hi
  )
}

registerDoParallel(cores = n_workers)
sim_out <- tryCatch(
  foreach(s = seq_len(n_sim), .errorhandling = "pass") %dopar% {
    run_one_sim(s)
  },
  finally = {
    stopImplicitCluster()
  }
)

ok <- list()
for (i in seq_along(sim_out)) {
  oi <- sim_out[[i]]
  if (inherits(oi, "error") || inherits(oi, "simpleError")) {
    cat("Worker error index", i, ":", conditionMessage(oi), "\n")
    next
  }
  if (!is.list(oi) || is.null(oi$fe_est)) next
  ok[[length(ok) + 1L]] <- oi
}
if (!length(ok)) stop("All JAGS sims failed.")

n_ok <- length(ok)
nm_fe <- fe_names()
fe_est_m <- do.call(rbind, lapply(ok, function(z) z$fe_est))
fe_lo_m <- do.call(rbind, lapply(ok, function(z) z$fe_lo))
fe_hi_m <- do.call(rbind, lapply(ok, function(z) z$fe_hi))
sig_est_m <- do.call(rbind, lapply(ok, function(z) z$sig_est))
sig_lo_m <- do.call(rbind, lapply(ok, function(z) z$sig_lo))
sig_hi_m <- do.call(rbind, lapply(ok, function(z) z$sig_hi))
sig_true_m <- do.call(rbind, lapply(ok, function(z) z$truth_sigma))
rho_est_m <- do.call(rbind, lapply(ok, function(z) z$rho_est))
rho_lo_m <- do.call(rbind, lapply(ok, function(z) z$rho_lo))
rho_hi_m <- do.call(rbind, lapply(ok, function(z) z$rho_hi))
rho_true_m <- do.call(rbind, lapply(ok, function(z) z$truth_rho))
truth_fe <- ok[[1]]$truth_fe

fe_stats <- calculate_stats(fe_est_m, truth_fe, fe_lo_m, fe_hi_m)
sig_stats <- calculate_stats_varying_true(sig_est_m, sig_true_m, sig_lo_m, sig_hi_m)
rho_stats <- calculate_stats_varying_true(rho_est_m, rho_true_m, rho_lo_m, rho_hi_m)

rep_tab <- data.frame(
  sim = vapply(ok, function(z) z$sim, 1L),
  N = n,
  Sec = vapply(ok, function(z) z$sec, 1.0),
  AESS = vapply(ok, function(z) z$aess, 1.0),
  stringsAsFactors = FALSE
)
write.csv(rep_tab, file.path(OUT_DIR, "jags_replicates.csv"), row.names = FALSE)

param_summary <- rbind(
  data.frame(
    group = "fixed", param = nm_fe, True = truth_fe,
    Bias = round(fe_stats$bias, 4), RMSE = round(fe_stats$rmse, 4), CP = round(fe_stats$cp, 3),
    stringsAsFactors = FALSE
  ),
  data.frame(
    group = "Sigma", param = sigma_vech_names, True = colMeans(sig_true_m),
    Bias = round(sig_stats$bias, 4), RMSE = round(sig_stats$rmse, 4), CP = round(sig_stats$cp, 3),
    stringsAsFactors = FALSE
  ),
  data.frame(
    group = "rho", param = rho_names, True = colMeans(rho_true_m),
    Bias = round(rho_stats$bias, 4), RMSE = round(rho_stats$rmse, 4), CP = round(rho_stats$cp, 3),
    stringsAsFactors = FALSE
  )
)
write.csv(param_summary, file.path(OUT_DIR, "jags_param_summary.csv"), row.names = FALSE)

timing <- data.frame(
  Method = "JAGS",
  N = n,
  Scenario = scenario,
  RE_DIST = re_dist,
  N_ok = n_ok,
  N_sim = n_sim,
  Mean_Sec = mean(rep_tab$Sec, na.rm = TRUE),
  Mean_AESS = mean(rep_tab$AESS, na.rm = TRUE),
  ESS_per_sec = mean(rep_tab$AESS / rep_tab$Sec, na.rm = TRUE),
  Adapt = jags_adapt,
  Burn = jags_burn,
  Iter = jags_iter,
  Thin = jags_thin,
  Chains = jags_chains,
  stringsAsFactors = FALSE
)
write.csv(timing, file.path(OUT_DIR, "jags_timing_summary.csv"), row.names = FALSE)

est_long <- do.call(rbind, lapply(ok, function(z) {
  data.frame(
    sim = z$sim,
    group = c(rep("fixed", length(nm_fe)), rep("Sigma", 6L), rep("rho", 3L)),
    param = c(nm_fe, sigma_vech_names, rho_names),
    truth = c(z$truth_fe, z$truth_sigma, z$truth_rho),
    est = c(z$fe_est, z$sig_est, z$rho_est),
    lo = c(z$fe_lo, z$sig_lo, z$rho_lo),
    hi = c(z$fe_hi, z$sig_hi, z$rho_hi),
    stringsAsFactors = FALSE
  )
}))
write.csv(est_long, file.path(OUT_DIR, "jags_param_estimates.csv"), row.names = FALSE)

cat("DONE JAGS n_ok=", n_ok, "/", n_sim, "\n", sep = "")
print(timing)
print(param_summary)
