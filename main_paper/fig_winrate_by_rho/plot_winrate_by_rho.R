############################################################
# Model-selection win rates by common RE correlation rho
# Joint CDPMM vs sum of separate (WAIC / LOOIC only)
# Defaults: Scenario 1, N=100, re_dist=normal, S=128
# Env: N, N_SIM, CHAIN, BURN, THIN, SCENARIO, RHO_GRID, OUT_DIR
############################################################

args <- commandArgs(trailingOnly = TRUE)
OUT_DIR <- if (length(args) >= 1L && nzchar(args[[1]])) {
  args[[1]]
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
if (!nzchar(as.character(OUT_DIR)[1])) {
  OUT_DIR <- normalizePath(".", winslash = "/", mustWork = TRUE)
}
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SKIP_MAIN_SIM <- TRUE
sim_candidates <- c("prog1_cdpmm_joint_vs_separate.R", "code.R")
sim_file <- sim_candidates[file.exists(sim_candidates)][1]
if (is.na(sim_file)) stop("Cannot find prog1_cdpmm_joint_vs_separate.R")
source(sim_file, encoding = "UTF-8", local = FALSE)

library(parallel)
library(doParallel)
library(foreach)

n_sim_wr <- as.integer(Sys.getenv("N_SIM", unset = "128"))
chain_wr <- as.integer(Sys.getenv("CHAIN", unset = "5000"))
burn_wr <- as.integer(Sys.getenv("BURN", unset = "2000"))
thin_wr <- as.integer(Sys.getenv("THIN", unset = "5"))
n_wr <- as.integer(Sys.getenv("N", unset = "100"))
scenario_wr <- as.integer(Sys.getenv("SCENARIO", unset = "1"))
re_dist_wr <- "normal"
if (!is.finite(n_wr) || n_wr < 1L) n_wr <- 100L
if (!is.finite(n_sim_wr) || n_sim_wr < 1L) n_sim_wr <- 128L

rho_grid_env <- Sys.getenv("RHO_GRID", unset = "")
if (nzchar(rho_grid_env)) {
  rho_grid <- as.numeric(strsplit(rho_grid_env, "[,;[:space:]]+")[[1]])
  rho_grid <- rho_grid[is.finite(rho_grid)]
} else {
  rho_grid <- seq(-1, 1, by = 0.1)
}
if (!length(rho_grid)) stop("Empty RHO_GRID")

n_cores <- suppressWarnings(as.integer(Sys.getenv(
  c("N_CORES", "SLURM_NTASKS_PER_NODE", "SLURM_CPUS_PER_TASK")
)))
n_cores <- n_cores[is.finite(n_cores) & n_cores >= 1L]
n_cores <- if (length(n_cores)) n_cores[[1]] else max(1L, parallel::detectCores() - 1L)
n_cores <- min(n_cores, n_sim_wr)
doParallel::registerDoParallel(cores = n_cores)

cat(sprintf(
  "Win-rate by rho | scenario=%d N=%d re_dist=%s S=%d chain=%d/%d thin=%d cores=%d\n",
  scenario_wr, n_wr, re_dist_wr, n_sim_wr, chain_wr, burn_wr, thin_wr, n_cores
))
cat("rho grid:", paste(sprintf("%.2f", rho_grid), collapse = ", "), "\n")
cat("OUT_DIR=", OUT_DIR, "\n", sep = "")

csv_file <- file.path(
  OUT_DIR,
  sprintf(
    "winrate_by_rho_scen%d_%s_n%d_nsim%d.csv",
    scenario_wr, re_dist_wr, n_wr, n_sim_wr
  )
)
meta_file <- file.path(
  OUT_DIR,
  sprintf(
    "winrate_by_rho_scen%d_%s_n%d_nsim%d_meta.rds",
    scenario_wr, re_dist_wr, n_wr, n_sim_wr
  )
)
pdf_file <- file.path(
  OUT_DIR,
  sprintf(
    "win_rate_waic_looic_by_rho_scen%d_%s_n%d_nsim%d.pdf",
    scenario_wr, re_dist_wr, n_wr, n_sim_wr
  )
)

run_one_rho <- function(rho_use, seed_base) {
  cat(sprintf("\n=== rho = %.2f | %d reps ===\n", rho_use, n_sim_wr))
  reps <- foreach::foreach(
    i = seq_len(n_sim_wr),
    .packages = c("BayesLogit", "mvtnorm", "MCMCpack", "truncnorm", "loo", "coda")
  ) %dopar% {
    tryCatch({
      dat <- generate_data(
        n = n_wr, nis = nis_fixed, seed = seed_base + i,
        random_nis = random_nis, nis_range = nis_range,
        scenario = scenario_wr, C = C,
        re_dist = re_dist_wr, mvt_df = mvt_df,
        rho = rho_use
      )
      joint_fit <- fit_joint_model(
        dat, chain = chain_wr, burn = burn_wr, thin = thin_wr,
        delta_min = delta_min, delta_max = delta_max
      )
      ordinal_fit <- fit_ordinal_model(
        dat, chain = chain_wr, burn = burn_wr, thin = thin_wr,
        delta_min = delta_min, delta_max = delta_max
      )
      zinb_fit <- fit_zinb_model(
        dat, chain = chain_wr, burn = burn_wr, thin = thin_wr
      )
      list(
        i = i,
        waic_joint = joint_fit$waic,
        waic_sep = ordinal_fit$waic + zinb_fit$waic,
        looic_joint = joint_fit$looic,
        looic_sep = ordinal_fit$looic + zinb_fit$looic
      )
    }, error = function(e) {
      list(i = i, error = conditionMessage(e))
    })
  }

  ok <- vapply(reps, function(z) {
    is.null(z$error) &&
      is.finite(z$waic_joint) && is.finite(z$waic_sep) &&
      is.finite(z$looic_joint) && is.finite(z$looic_sep)
  }, logical(1))
  if (any(!ok)) {
    errs <- vapply(reps[!ok], function(z) {
      if (!is.null(z$error)) z$error else "non-finite criteria"
    }, character(1))
    cat("Failed reps for rho=", rho_use, ":\n", sep = "")
    print(utils::head(unique(errs), 5L))
  }
  reps <- reps[ok]
  n_ok <- length(reps)
  if (n_ok < 1L) {
    return(data.frame(
      rho = rho_use, n_success = 0L, n_sim = n_sim_wr,
      WAIC_WinPct = NA_real_, LOOIC_WinPct = NA_real_,
      WAIC_Mean_Joint = NA_real_, WAIC_Mean_SepSum = NA_real_,
      LOOIC_Mean_Joint = NA_real_, LOOIC_Mean_SepSum = NA_real_
    ))
  }

  waic_j <- vapply(reps, `[[`, numeric(1), "waic_joint")
  waic_s <- vapply(reps, `[[`, numeric(1), "waic_sep")
  loo_j <- vapply(reps, `[[`, numeric(1), "looic_joint")
  loo_s <- vapply(reps, `[[`, numeric(1), "looic_sep")

  out <- data.frame(
    rho = rho_use,
    n_success = n_ok,
    n_sim = n_sim_wr,
    WAIC_WinPct = mean(waic_j < waic_s) * 100,
    LOOIC_WinPct = mean(loo_j < loo_s) * 100,
    WAIC_Mean_Joint = mean(waic_j),
    WAIC_Mean_SepSum = mean(waic_s),
    LOOIC_Mean_Joint = mean(loo_j),
    LOOIC_Mean_SepSum = mean(loo_s)
  )
  cat(sprintf(
    "rho=%.2f success=%d/%d WAIC=%.1f%% LOOIC=%.1f%%\n",
    rho_use, n_ok, n_sim_wr, out$WAIC_WinPct, out$LOOIC_WinPct
  ))
  out
}

rows <- vector("list", length(rho_grid))
for (k in seq_along(rho_grid)) {
  rho_k <- rho_grid[[k]]
  # Distinct seed bases per rho so data sets differ across the grid
  seed_base <- 810000L + as.integer(round(rho_k * 1000)) * 1000L
  rows[[k]] <- run_one_rho(rho_k, seed_base)
  win_tab <- do.call(rbind, rows[seq_len(k)])
  utils::write.csv(win_tab, csv_file, row.names = FALSE)
  saveRDS(
    list(
      settings = list(
        scenario = scenario_wr, re_dist = re_dist_wr, n = n_wr,
        n_sim = n_sim_wr, chain = chain_wr, burn = burn_wr, thin = thin_wr,
        rho_grid = rho_grid
      ),
      winrate = win_tab
    ),
    meta_file
  )
}

doParallel::stopImplicitCluster()
win_tab <- do.call(rbind, rows)
utils::write.csv(win_tab, csv_file, row.names = FALSE)

############################################################
# Tang-style win-rate plot (WAIC orange, LOOIC blue)
############################################################
draw_winrate_plot <- function(tab, file) {
  tab <- tab[is.finite(tab$rho) & is.finite(tab$WAIC_WinPct) & is.finite(tab$LOOIC_WinPct), ]
  tab <- tab[order(tab$rho), ]
  if (!nrow(tab)) stop("No finite win-rate rows to plot")

  grDevices::pdf(file, width = 8.2, height = 5.2)
  op <- graphics::par(mar = c(4.2, 4.2, 1.2, 1.2), mgp = c(2.4, 0.7, 0))
  on.exit({
    graphics::par(op)
    grDevices::dev.off()
  }, add = FALSE)

  xlim <- range(tab$rho)
  ylim <- c(0, 100)
  plot(
    tab$rho, tab$WAIC_WinPct,
    type = "n",
    xlim = xlim, ylim = ylim,
    xlab = expression(Correlation ~ rho),
    ylab = "Win rate (%)",
    xaxt = "n", yaxt = "n",
    xaxs = "i", yaxs = "i"
  )
  graphics::axis(1, at = seq(-1, 1, by = 0.2))
  graphics::axis(2, at = seq(0, 100, by = 10), las = 1)
  graphics::abline(h = seq(0, 100, by = 10), col = "gray85", lwd = 1)
  graphics::abline(v = 0, col = "gray50", lty = 2, lwd = 1.2)

  col_waic <- "#E69F00"
  col_loo <- "#0072B2"
  graphics::lines(tab$rho, tab$WAIC_WinPct, col = col_waic, lwd = 2)
  graphics::points(tab$rho, tab$WAIC_WinPct, pch = 16, col = col_waic, cex = 1.05)
  graphics::lines(tab$rho, tab$LOOIC_WinPct, col = col_loo, lwd = 2)
  graphics::points(tab$rho, tab$LOOIC_WinPct, pch = 16, col = col_loo, cex = 1.05)

  graphics::legend(
    "top",
    legend = c("WAIC Win (%)", "LOOIC Win (%)"),
    col = c(col_waic, col_loo),
    lty = 1, lwd = 2, pch = 16, pt.cex = 1.05,
    bty = "n", horiz = TRUE, cex = 0.95
  )
  graphics::box()
}

draw_winrate_plot(win_tab, pdf_file)
cat("Wrote:\n  ", csv_file, "\n  ", pdf_file, "\n  ", meta_file, "\n", sep = "")
