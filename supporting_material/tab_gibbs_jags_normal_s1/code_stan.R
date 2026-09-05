############################################################
# Stan-only Monte Carlo (scheme 1: ~1000 retained draws)
# Env: SCENARIO, N / N_SUBJ, RE_DIST, NSIM, OUT_DIR
# Seeds match main Gibbs table: generate_data(..., seed = 2026 + s)
############################################################

.env_int <- function(key, default) {
  v <- suppressWarnings(as.integer(Sys.getenv(key, unset = "")))
  if (length(v) == 1L && is.finite(v)) v else as.integer(default)
}
.env_num <- function(key, default) {
  v <- suppressWarnings(as.numeric(Sys.getenv(key, unset = "")))
  if (length(v) == 1L && is.finite(v)) v else as.numeric(default)
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

stan_warmup <- .env_int("STAN_WARMUP", 500L)
stan_sampling <- .env_int("STAN_SAMPLING", 500L)
stan_chains <- .env_int("STAN_CHAINS", 2L)
stan_refresh <- .env_int("STAN_REFRESH", 0L)
# Force chain-serial when sims are forked (avoid oversubscribe)
stan_par_chains_force <- .env_int("STAN_PARALLEL_CHAINS", 1L)

cmdstan_root <- .env_chr("CMDSTAN", path.expand("~/cmdstan"))

cat("METHOD=Stan OUT_DIR=", OUT_DIR, "\n")
cat("NSIM=", n_sim, " N=", n, " G=", G_mix,
    " SCENARIO=", scenario, " RE_DIST=", re_dist, "\n", sep = "")
cat("Stan:", stan_chains, "chains,", stan_warmup, "+", stan_sampling,
    " parallel_chains=", stan_par_chains_force, "\n", sep = "")

user_lib <- path.expand("~/R_libs")
dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(user_lib, .libPaths()))

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(coda)
  library(foreach)
  library(doParallel)
})

cands <- list.files(cmdstan_root, pattern = "^cmdstan", full.names = TRUE)
if (!length(cands)) {
  stop("CmdStan not found under ", cmdstan_root, ". Install CmdStan first.")
}
dts_bin <- "/opt/rh/devtoolset-7/root/usr/bin"
cxx <- file.path(dts_bin, "g++")
cc <- file.path(dts_bin, "gcc")
if (!file.exists(cxx)) {
  cxx <- "g++"
  cc <- "gcc"
}
tbb_dir <- file.path(cands[[1]], "stan", "lib", "stan_math", "lib", "tbb")
lp <- Sys.getenv("LIBRARY_PATH", unset = "")
ldp <- Sys.getenv("LD_LIBRARY_PATH", unset = "")
Sys.setenv(
  TBB_CXX_TYPE = "gcc",
  CXX = cxx,
  CC = cc,
  PATH = paste0(dts_bin, ":", Sys.getenv("PATH")),
  LIBRARY_PATH = if (nzchar(lp)) paste0(tbb_dir, ":", lp) else tbb_dir,
  LD_LIBRARY_PATH = if (nzchar(ldp)) paste0(tbb_dir, ":", ldp) else tbb_dir
)
set_cmdstan_path(cands[[1]])
cat("cmdstan_path=", cmdstan_path(), " version=", as.character(cmdstan_version()), "\n")

local_mk <- file.path(cmdstan_path(), "make", "local")
writeLines(c(
  paste0("CXX = ", cxx),
  paste0("CC = ", cc),
  "TBB_CXX_TYPE = gcc",
  paste0("LDFLAGS += -L", tbb_dir, " -Wl,-rpath,", tbb_dir)
), local_mk)

stan_file <- file.path(OUT_DIR, "joint_cdpmm.stan")
if (!file.exists(stan_file)) stan_file <- "joint_cdpmm.stan"
mod <- cmdstan_model(
  stan_file,
  force_recompile = FALSE,
  cpp_options = list(CXX = cxx, CC = cc, TBB_CXX_TYPE = "gcc")
)
cat("Stan model ready\n")

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
rho1 <- 0.5

as_stan_data <- function(dat) {
  list(
    N = dat$n,
    Ntot = dat$N,
    C = dat$C,
    G = G_mix,
    pz = ncol(dat$X_zero),
    pc = ncol(dat$X_count),
    po = ncol(dat$X_ordinal),
    nis = as.integer(dat$nis),
    id = as.integer(dat$id),
    y1 = as.integer(dat$y1),
    y2 = as.integer(dat$y2),
    Xz = dat$X_zero,
    Xc = dat$X_count,
    Xo = dat$X_ordinal,
    prior_sd = sqrt(1000),
    tau_a1 = 1,
    tau_a2 = 1
  )
}

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

n_cores <- as.integer(Sys.getenv(
  c("SLURM_NTASKS_PER_NODE", "SLURM_CPUS_PER_TASK", "N_CORES"),
  unset = NA_character_
))
n_cores <- n_cores[is.finite(n_cores) & n_cores >= 1L]
n_cores <- if (length(n_cores)) n_cores[[1]] else parallel::detectCores()
n_cores <- max(1L, as.integer(n_cores))
n_workers <- max(1L, min(as.integer(n_sim), n_cores))
stan_par_chains <- max(1L, min(as.integer(stan_par_chains_force), as.integer(stan_chains)))
cat("Parallel sims: n_workers=", n_workers, " of n_cores=", n_cores,
    " Stan parallel_chains=", stan_par_chains, "\n", sep = "")

run_one_sim <- function(s) {
  seed_s <- 2026L + as.integer(s)  # match main Gibbs Monte Carlo
  cat("\n===== Stan replicate", s, "/", n_sim, " pid=", Sys.getpid(), " =====\n", sep = "")
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

  sd <- as_stan_data(dat)
  t0 <- proc.time()[["elapsed"]]
  sfit <- tryCatch(
    mod$sample(
      data = sd,
      seed = seed_s,
      chains = stan_chains,
      parallel_chains = stan_par_chains,
      iter_warmup = stan_warmup,
      iter_sampling = stan_sampling,
      refresh = stan_refresh,
      adapt_delta = 0.9,
      max_treedepth = 12,
      show_messages = FALSE
    ),
    error = function(e) {
      cat("Stan failed (sim=", s, "):", conditionMessage(e), "\n")
      NULL
    }
  )
  sec <- proc.time()[["elapsed"]] - t0

  n_fe <- length(truth)
  fe_est <- fe_lo <- fe_hi <- rep(NA_real_, n_fe)
  sig_est <- sig_lo <- sig_hi <- rep(NA_real_, 6L)
  rho_est <- rho_lo <- rho_hi <- rep(NA_real_, 3L)
  aess <- n_div <- NA_real_

  if (!is.null(sfit)) {
    dr_fe <- sfit$draws(variables = c("alpha", "beta", "gamma", "r"), format = "draws_matrix")
    fe_est <- as.numeric(colMeans(dr_fe))
    fe_lo <- as.numeric(apply(dr_fe, 2, quantile, probs = 0.025, na.rm = TRUE))
    fe_hi <- as.numeric(apply(dr_fe, 2, quantile, probs = 0.975, na.rm = TRUE))
    aess <- param_ess_mean(dr_fe)

    S_dr <- sfit$draws(variables = "Sigma", format = "draws_matrix")
    # CmdStan may flatten Sigma as Sigma[1,1],...
    sig_draws <- matrix(NA_real_, nrow(S_dr), 6L)
    for (t in seq_len(nrow(S_dr))) {
      Sm <- matrix(as.numeric(S_dr[t, ]), 3, 3)
      sig_draws[t, ] <- vech3(Sm)
    }
    sig_est <- colMeans(sig_draws)
    sig_lo <- apply(sig_draws, 2, quantile, probs = 0.025, na.rm = TRUE)
    sig_hi <- apply(sig_draws, 2, quantile, probs = 0.975, na.rm = TRUE)

    r_dr <- sfit$draws(variables = "rho", format = "draws_matrix")
    rho_est <- as.numeric(colMeans(r_dr))
    rho_lo <- as.numeric(apply(r_dr, 2, quantile, probs = 0.025, na.rm = TRUE))
    rho_hi <- as.numeric(apply(r_dr, 2, quantile, probs = 0.975, na.rm = TRUE))

    diag_sum <- tryCatch(sfit$diagnostic_summary(), error = function(e) NULL)
    n_div <- if (!is.null(diag_sum) && !is.null(diag_sum$num_divergent))
      sum(diag_sum$num_divergent, na.rm = TRUE) else NA_real_
  }

  list(
    sim = s,
    sec = sec,
    aess = aess,
    n_div = n_div,
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
if (!length(ok)) stop("All Stan sims failed.")

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
  Divergent = vapply(ok, function(z) z$n_div, 1.0),
  stringsAsFactors = FALSE
)
write.csv(rep_tab, file.path(OUT_DIR, "stan_replicates.csv"), row.names = FALSE)

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
write.csv(param_summary, file.path(OUT_DIR, "stan_param_summary.csv"), row.names = FALSE)

timing <- data.frame(
  Method = "Stan",
  N = n,
  Scenario = scenario,
  RE_DIST = re_dist,
  N_ok = n_ok,
  N_sim = n_sim,
  Mean_Sec = mean(rep_tab$Sec, na.rm = TRUE),
  Mean_AESS = mean(rep_tab$AESS, na.rm = TRUE),
  ESS_per_sec = mean(rep_tab$AESS / rep_tab$Sec, na.rm = TRUE),
  Mean_Divergent = mean(rep_tab$Divergent, na.rm = TRUE),
  Warmup = stan_warmup,
  Sampling = stan_sampling,
  Chains = stan_chains,
  stringsAsFactors = FALSE
)
write.csv(timing, file.path(OUT_DIR, "stan_timing_summary.csv"), row.names = FALSE)

# Per-sim point estimates (for optional pooling later)
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
write.csv(est_long, file.path(OUT_DIR, "stan_param_estimates.csv"), row.names = FALSE)

cat("DONE Stan n_ok=", n_ok, "/", n_sim, "\n", sep = "")
print(timing)
print(param_summary)
