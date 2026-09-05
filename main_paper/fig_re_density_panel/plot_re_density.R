############################################################
# Tang-style RE density panel figure
# (a) bi truth = normal: Scenario 1/2 x (b1,b2,b3)
# (b) bi truth = mixture: Scenario 1/2 x (b1,b2,b3)
# Defaults: N=100, N_SIM_PLOT=128, true density + KDE of posterior means
# Override via env: CHAIN, BURN, THIN, N_CORES, OUT_DIR / CLI OUT_DIR
############################################################

args <- commandArgs(trailingOnly = TRUE)
OUT_DIR <- if (length(args) >= 1) args[[1]] else {
  code_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)
  parent <- dirname(code_dir)
  base <- basename(code_dir)
  if (grepl("[\u4e00-\u9fff]", base)) {
    file.path(parent, paste0(base, "结果"))
  } else {
    file.path(parent, paste0(base, "_results"))
  }
}
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SKIP_MAIN_SIM <- TRUE
sim_candidates <- c("prog1_cdpmm_joint_vs_separate.R", "code.R", "simulation_code_hpc.R")
sim_file <- sim_candidates[file.exists(sim_candidates)][1]
if (is.na(sim_file)) stop("Cannot find prog1_cdpmm_joint_vs_separate.R")
source(sim_file, encoding = "UTF-8", local = FALSE)

library(parallel)
library(doParallel)
library(foreach)

# Paper-oriented defaults (override on smoke tests via env)
n_sim_plot <- as.integer(Sys.getenv("N_SIM_PLOT", unset = "128"))
chain_plot <- as.integer(Sys.getenv("CHAIN", unset = "5000"))
burn_plot <- as.integer(Sys.getenv("BURN", unset = "2000"))
thin_plot <- as.integer(Sys.getenv("THIN", unset = "5"))
n_plot <- as.integer(Sys.getenv("N", unset = "100"))
if (!is.finite(n_plot) || n_plot < 1L) n_plot <- 100L
L_grid <- 50L

############################################################
# True marginal densities for each component of b
# Ordering: (b1,b2,b3) = (zero, count, ordinal)
############################################################
true_marginal_density <- function(x, m, re_dist, Sigma_true, re_info = NULL,
                                  mvt_df = 4) {
  if (identical(re_dist, "mixture")) {
    mix <- generate_mixture_re(n = 2L)
    dens <- numeric(length(x))
    for (g in seq_along(mix$pi)) {
      dens <- dens + mix$pi[g] * stats::dnorm(
        x,
        mean = mix$mu[g, m],
        sd = sqrt(max(mix$Omega_list[[g]][m, m], 1e-10))
      )
    }
    return(dens)
  }
  if (identical(re_dist, "normal")) {
    return(stats::dnorm(x, mean = 0, sd = sqrt(max(Sigma_true[m, m], 1e-10))))
  }
  if (identical(re_dist, "mvt")) {
    sc <- if (!is.null(re_info$scale_Sigma)) {
      re_info$scale_Sigma[m, m]
    } else {
      ((mvt_df - 2) / mvt_df) * Sigma_true[m, m]
    }
    sc <- max(sc, 1e-10)
    df_use <- if (!is.null(re_info$df)) re_info$df else mvt_df
    return(stats::dt(x / sqrt(sc), df = df_use) / sqrt(sc))
  }
  stop("Unknown re_dist")
}

rmse_density_b <- function(b_hat, re_dist, Sigma_true, re_info = NULL,
                           mvt_df = 4, L = L_grid) {
  q <- ncol(b_hat)
  acc <- 0
  for (m in seq_len(q)) {
    bh <- b_hat[, m]
    bh <- bh[is.finite(bh)]
    if (length(bh) < 5L) return(NA_real_)
    probs <- seq(0.01, 0.99, length.out = L)
    h <- as.numeric(stats::quantile(bh, probs = probs, names = FALSE, na.rm = TRUE))
    kd <- stats::density(bh, n = 512)
    p_hat <- stats::approx(kd$x, kd$y, xout = h, rule = 2)$y
    p_true <- true_marginal_density(h, m, re_dist, Sigma_true, re_info, mvt_df)
    acc <- acc + sum((p_true - p_hat)^2)
  }
  sqrt(acc / (q * L))
}

run_one_setting <- function(scenario_use, re_dist_use, seed_base) {
  cat(sprintf(
    "\n=== Setting: scenario=%d, re_dist=%s | %d reps, n=%d, chain=%d/%d thin=%d ===\n",
    scenario_use, re_dist_use, n_sim_plot, n_plot, chain_plot, burn_plot, thin_plot
  ))

  plot_reps <- foreach::foreach(
    i = seq_len(n_sim_plot),
    .packages = c("BayesLogit", "mvtnorm", "MCMCpack", "truncnorm", "loo", "coda")
  ) %dopar% {
    tryCatch({
      dat <- generate_data(
        n = n_plot,
        nis = nis_fixed,
        seed = seed_base + i,
        random_nis = random_nis,
        nis_range = nis_range,
        scenario = scenario_use,
        C = C,
        re_dist = re_dist_use,
        mvt_df = mvt_df
      )
      fit <- fit_joint_model(
        dat,
        chain = chain_plot,
        burn = burn_plot,
        thin = thin_plot,
        delta_min = delta_min,
        delta_max = delta_max
      )
      rmse <- rmse_density_b(
        fit$b_est, dat$re_dist, dat$Sigma_true, dat$re_info, mvt_df
      )
      list(
        i = i,
        b_est = fit$b_est,
        Sigma_true = dat$Sigma_true,
        re_dist = dat$re_dist,
        re_info = dat$re_info,
        rmse = rmse
      )
    }, error = function(e) {
      list(i = i, error = conditionMessage(e))
    })
  }

  ok <- vapply(plot_reps, function(z) is.null(z$error) && is.finite(z$rmse), logical(1))
  if (any(!ok)) {
    errs <- vapply(plot_reps[!ok], function(z) {
      if (!is.null(z$error)) z$error else "non-finite rmse"
    }, character(1))
    cat("Failed reps:\n"); print(unique(errs))
  }
  plot_reps <- plot_reps[ok]
  if (!length(plot_reps)) {
    stop("All replications failed for scenario=", scenario_use, " re_dist=", re_dist_use)
  }
  rmse_vec <- vapply(plot_reps, `[[`, numeric(1), "rmse")
  pick <- which.min(rmse_vec)
  rep_plot <- plot_reps[[pick]]
  cat(sprintf(
    "Successful: %d/%d; min dens-RMSE=%.4f; plotting rep #%d\n",
    length(plot_reps), n_sim_plot, rep_plot$rmse, rep_plot$i
  ))
  list(
    scenario = scenario_use,
    re_dist = re_dist_use,
    n_success = length(plot_reps),
    rmse_vec = rmse_vec,
    pick_rep_id = rep_plot$i,
    min_rmse = rep_plot$rmse,
    rep = rep_plot
  )
}

############################################################
# Parallel cluster once; four settings sequentially
############################################################
n_cores <- suppressWarnings(as.integer(Sys.getenv(
  c("N_CORES", "SLURM_NTASKS_PER_NODE", "SLURM_CPUS_PER_TASK")
)))
n_cores <- n_cores[is.finite(n_cores) & n_cores >= 1L]
n_cores <- if (length(n_cores)) n_cores[[1]] else max(1L, parallel::detectCores() - 1L)
n_cores <- min(n_cores, n_sim_plot)
doParallel::registerDoParallel(cores = n_cores)
cat(sprintf("Registered %d cores for RE density panel (N_SIM_PLOT=%d)\n", n_cores, n_sim_plot))

settings <- list(
  list(scenario = 1L, re_dist = "normal",  seed_base = 91000L),
  list(scenario = 2L, re_dist = "normal",  seed_base = 92000L),
  list(scenario = 1L, re_dist = "mixture", seed_base = 93000L),
  list(scenario = 2L, re_dist = "mixture", seed_base = 94000L)
)

results <- vector("list", length(settings))
for (k in seq_along(settings)) {
  s <- settings[[k]]
  results[[k]] <- run_one_setting(s$scenario, s$re_dist, s$seed_base)
  names(results)[k] <- paste0("s", s$scenario, "_", s$re_dist)
}

doParallel::stopImplicitCluster()

saveRDS(
  list(
    settings = list(
      n = n_plot, n_sim_plot = n_sim_plot,
      chain = chain_plot, burn = burn_plot, thin = thin_plot,
      designs = c("normal x S1/S2", "mixture x S1/S2")
    ),
    results = results
  ),
  file.path(OUT_DIR, paste0(
    "re_density_panel_meta_n", n_plot, "_nsim", n_sim_plot,
    "_chain", chain_plot, ".rds"
  ))
)

############################################################
# Draw combined PDF: (a) normal 2x3, (b) mixture 2x3
############################################################
panel_labs <- list(
  expression(b[1]~"(zero)"),
  expression(b[2]~"(count)"),
  expression(b[3]~"(ordinal)")
)

draw_one_panel <- function(m, rep_plot, show_legend = FALSE, ylab = "Density") {
  bh <- rep_plot$b_est[, m]
  kd <- stats::density(bh, n = 512)
  xg <- kd$x
  yt <- true_marginal_density(
    xg, m, rep_plot$re_dist, rep_plot$Sigma_true, rep_plot$re_info, mvt_df
  )
  ylim <- range(c(kd$y, yt), finite = TRUE)
  ylim[1] <- 0
  ylim[2] <- ylim[2] * 1.08
  plot(
    xg, yt,
    type = "l", lwd = 2, col = "black",
    xlab = panel_labs[[m]], ylab = ylab,
    ylim = ylim, main = ""
  )
  lines(kd$x, kd$y, col = "red", lty = 2, lwd = 2)
  if (isTRUE(show_legend)) {
    legend(
      "topright",
      legend = c("True", "Estimated (KDE)"),
      col = c("black", "red"), lty = c(1, 2), lwd = 2,
      bty = "n", cex = 0.8
    )
  }
}

get_res <- function(scen, redist) {
  results[[paste0("s", scen, "_", redist)]]$rep
}

pdf_file <- file.path(
  OUT_DIR,
  paste0("re_density_true_vs_est_panel_ab_n", n_plot, "_nsim", n_sim_plot, ".pdf")
)

# layout: title (a), normal 2x3, title (b), mixture 2x3
grDevices::pdf(pdf_file, width = 10.5, height = 9.5)
lay <- matrix(c(
  1, 1, 1,
  2, 3, 4,
  5, 6, 7,
  8, 8, 8,
  9, 10, 11,
  12, 13, 14
), nrow = 6, ncol = 3, byrow = TRUE)
graphics::layout(lay, heights = c(0.22, 1, 1, 0.22, 1, 1))

graphics::par(mar = c(0, 0, 0, 0))
plot.new()
graphics::text(0.02, 0.5, "(a) Normal random-effects truth",
               adj = 0, cex = 1.25, font = 2, xpd = NA)

graphics::par(mar = c(4.0, 4.0, 1.0, 0.6))
for (scen in c(1L, 2L)) {
  rep_plot <- get_res(scen, "normal")
  for (m in 1:3) {
    draw_one_panel(
      m, rep_plot,
      show_legend = (scen == 1L && m == 1L),
      ylab = if (m == 1L) paste0("Scenario ", scen) else "Density"
    )
  }
}

graphics::par(mar = c(0, 0, 0, 0))
plot.new()
graphics::text(0.02, 0.5, "(b) Mixture random-effects truth",
               adj = 0, cex = 1.25, font = 2, xpd = NA)

graphics::par(mar = c(4.0, 4.0, 1.0, 0.6))
for (scen in c(1L, 2L)) {
  rep_plot <- get_res(scen, "mixture")
  for (m in 1:3) {
    draw_one_panel(
      m, rep_plot,
      show_legend = FALSE,
      ylab = if (m == 1L) paste0("Scenario ", scen) else "Density"
    )
  }
}

grDevices::dev.off()

cat("Wrote:\n  ", pdf_file, "\n", sep = "")
