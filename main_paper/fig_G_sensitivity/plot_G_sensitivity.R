############################################################
# Truncation G sensitivity: ARMSE / ACP vs G
# Defaults: Scenario 2, N=100, re_dist via RE_DIST (normal|mixture)
# G grid: 5,8,12,15,20,25,30,50,80 ; S=100 (joint CDPMM only)
# ARMSE/ACP average over alpha,beta,gamma,r, Sigma diag, rho (exclude delta)
# Env: N, N_SIM, CHAIN, BURN, THIN, SCENARIO, RE_DIST, G_GRID, OUT_DIR
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

n_sim_g <- as.integer(Sys.getenv("N_SIM", unset = "200"))
chain_g <- as.integer(Sys.getenv("CHAIN", unset = "5000"))
burn_g <- as.integer(Sys.getenv("BURN", unset = "2000"))
thin_g <- as.integer(Sys.getenv("THIN", unset = "5"))
n_g <- as.integer(Sys.getenv("N", unset = "100"))
scenario_g <- as.integer(Sys.getenv("SCENARIO", unset = "2"))
re_dist_g <- Sys.getenv("RE_DIST", unset = "normal")
if (!nzchar(re_dist_g)) re_dist_g <- "normal"
re_dist_g <- tolower(re_dist_g)
if (!re_dist_g %in% c("normal", "mixture")) {
  stop("RE_DIST must be 'normal' or 'mixture', got: ", re_dist_g)
}
if (!is.finite(n_g) || n_g < 1L) n_g <- 100L
if (!is.finite(n_sim_g) || n_sim_g < 1L) n_sim_g <- 200L

g_env <- Sys.getenv("G_GRID", unset = "")
G_grid <- if (nzchar(g_env)) {
  as.integer(as.numeric(strsplit(g_env, "[,;[:space:]]+")[[1]]))
} else {
  c(5L, 8L, 12L, 15L, 20L, 25L, 30L, 50L, 80L)
}
G_grid <- unique(G_grid[is.finite(G_grid) & G_grid >= 2L])
if (!length(G_grid)) stop("Empty G_GRID")

n_cores <- suppressWarnings(as.integer(Sys.getenv(
  c("N_CORES", "SLURM_NTASKS_PER_NODE", "SLURM_CPUS_PER_TASK")
)))
n_cores <- n_cores[is.finite(n_cores) & n_cores >= 1L]
n_cores <- if (length(n_cores)) n_cores[[1]] else max(1L, parallel::detectCores() - 1L)
n_cores <- min(n_cores, n_sim_g)
doParallel::registerDoParallel(cores = n_cores)

cat(sprintf(
  "G-sensitivity | scenario=%d N=%d re_dist=%s S=%d chain=%d/%d thin=%d cores=%d\n",
  scenario_g, n_g, re_dist_g, n_sim_g, chain_g, burn_g, thin_g, n_cores
))
cat("G grid:", paste(G_grid, collapse = ", "), "\n")
cat("OUT_DIR=", OUT_DIR, "\n", sep = "")

tag <- sprintf("scen%d_%s_n%d_nsim%d", scenario_g, re_dist_g, n_g, n_sim_g)
csv_file <- file.path(OUT_DIR, paste0("G_sensitivity_", tag, ".csv"))
meta_file <- file.path(OUT_DIR, paste0("G_sensitivity_", tag, "_meta.rds"))
pdf_file <- file.path(OUT_DIR, paste0("G_sensitivity_armse_acp_", tag, ".pdf"))

# True FE: Scenario 1 is 5/6/4; Scenario 2 (no I(U>0)) is 4/5/3
fe_true <- true_fixed_effects(scenario_g)
alpha_true <- fe_true$alpha
beta_true <- fe_true$beta
gamma_true <- fe_true$gamma
r_true <- 2
cat(sprintf(
  "FE dims: alpha=%d beta=%d gamma=%d (scenario %d)\n",
  length(alpha_true), length(beta_true), length(gamma_true), scenario_g
))

extract_armse_acp <- function(fits, Sigma_true_ref, rho_true_ref) {
  n_ok <- length(fits)
  p_a <- length(alpha_true)
  p_b <- length(beta_true)
  p_g <- length(gamma_true)

  a_est <- matrix(NA_real_, n_ok, p_a)
  a_lo <- matrix(NA_real_, n_ok, p_a)
  a_hi <- matrix(NA_real_, n_ok, p_a)
  b_est <- matrix(NA_real_, n_ok, p_b)
  b_lo <- matrix(NA_real_, n_ok, p_b)
  b_hi <- matrix(NA_real_, n_ok, p_b)
  g_est <- matrix(NA_real_, n_ok, p_g)
  g_lo <- matrix(NA_real_, n_ok, p_g)
  g_hi <- matrix(NA_real_, n_ok, p_g)
  r_est <- rep(NA_real_, n_ok)
  r_lo <- rep(NA_real_, n_ok)
  r_hi <- rep(NA_real_, n_ok)
  s_est <- matrix(NA_real_, n_ok, 3L)
  s_lo <- matrix(NA_real_, n_ok, 3L)
  s_hi <- matrix(NA_real_, n_ok, 3L)
  rho_est <- matrix(NA_real_, n_ok, 3L)
  rho_lo <- matrix(NA_real_, n_ok, 3L)
  rho_hi <- matrix(NA_real_, n_ok, 3L)
  s_true <- matrix(NA_real_, n_ok, 3L)
  rho_true_mat <- matrix(NA_real_, n_ok, 3L)

  for (i in seq_len(n_ok)) {
    fit <- fits[[i]]$fit
    Sig_t <- fits[[i]]$Sigma_true
    a_est[i, ] <- fit$alpha_est
    a_lo[i, ] <- fit$alpha_ci_lower
    a_hi[i, ] <- fit$alpha_ci_upper
    b_est[i, ] <- fit$beta_est
    b_lo[i, ] <- fit$beta_ci_lower
    b_hi[i, ] <- fit$beta_ci_upper
    g_est[i, ] <- fit$gamma_est
    g_lo[i, ] <- fit$gamma_ci_lower
    g_hi[i, ] <- fit$gamma_ci_upper
    r_est[i] <- fit$r_est
    r_lo[i] <- fit$r_ci_lower
    r_hi[i] <- fit$r_ci_upper
    s_est[i, ] <- diag(fit$Sigma_est)
    s_lo[i, ] <- diag(fit$Sigma_ci_lower)
    s_hi[i, ] <- diag(fit$Sigma_ci_upper)
    rho_est[i, ] <- fit$Rho_est
    rho_lo[i, ] <- fit$Rho_ci_lower
    rho_hi[i, ] <- fit$Rho_ci_upper
    s_true[i, ] <- diag(Sig_t)
    rho_true_mat[i, ] <- c(
      Sig_t[1, 2] / sqrt(Sig_t[1, 1] * Sig_t[2, 2]),
      Sig_t[1, 3] / sqrt(Sig_t[1, 1] * Sig_t[3, 3]),
      Sig_t[2, 3] / sqrt(Sig_t[2, 2] * Sig_t[3, 3])
    )
  }

  st_a <- calculate_stats(a_est, alpha_true, a_lo, a_hi)
  st_b <- calculate_stats(b_est, beta_true, b_lo, b_hi)
  st_g <- calculate_stats(g_est, gamma_true, g_lo, g_hi)
  st_r <- calculate_stats(matrix(r_est, ncol = 1), r_true, matrix(r_lo, ncol = 1), matrix(r_hi, ncol = 1))
  st_s <- calculate_stats_varying_true(s_est, s_true, s_lo, s_hi)
  st_rho <- calculate_stats_varying_true(rho_est, rho_true_mat, rho_lo, rho_hi)

  rmse_vec <- c(st_a$rmse, st_b$rmse, st_g$rmse, st_r$rmse, st_s$rmse, st_rho$rmse)
  cp_vec <- c(st_a$cp, st_b$cp, st_g$cp, st_r$cp, st_s$cp, st_rho$cp)
  list(
    ARMSE = mean(rmse_vec, na.rm = TRUE),
    ACP = mean(cp_vec, na.rm = TRUE),
    n_param = length(rmse_vec),
    Sigma_true_ref = Sigma_true_ref,
    rho_true_ref = rho_true_ref
  )
}

run_one_G <- function(G_use, seed_base) {
  cat(sprintf("\n=== G=%d | %d reps | re_dist=%s ===\n", G_use, n_sim_g, re_dist_g))
  reps <- foreach::foreach(
    i = seq_len(n_sim_g),
    .packages = c("BayesLogit", "mvtnorm", "MCMCpack", "truncnorm", "loo", "coda")
  ) %dopar% {
    tryCatch({
      dat <- generate_data(
        n = n_g, nis = nis_fixed, seed = seed_base + i,
        random_nis = random_nis, nis_range = nis_range,
        scenario = scenario_g, C = C,
        re_dist = re_dist_g, mvt_df = mvt_df
      )
      fit <- fit_joint_model(
        dat, chain = chain_g, burn = burn_g, thin = thin_g,
        delta_min = delta_min, delta_max = delta_max,
        G = G_use
      )
      list(i = i, fit = fit, Sigma_true = dat$Sigma_true)
    }, error = function(e) {
      list(i = i, error = conditionMessage(e))
    })
  }

  ok <- vapply(reps, function(z) is.null(z$error) && !is.null(z$fit), logical(1))
  if (any(!ok)) {
    errs <- vapply(reps[!ok], function(z) {
      if (!is.null(z$error)) z$error else "missing fit"
    }, character(1))
    cat("Failed reps:\n")
    print(utils::head(unique(errs), 5L))
  }
  reps <- reps[ok]
  n_ok <- length(reps)
  if (n_ok < 1L) {
    return(data.frame(
      G = G_use, re_dist = re_dist_g, n_success = 0L, n_sim = n_sim_g,
      ARMSE = NA_real_, ACP = NA_real_, n_param = NA_integer_
    ))
  }

  Sigma_ref <- reps[[1]]$Sigma_true
  rho_ref <- c(
    Sigma_ref[1, 2] / sqrt(Sigma_ref[1, 1] * Sigma_ref[2, 2]),
    Sigma_ref[1, 3] / sqrt(Sigma_ref[1, 1] * Sigma_ref[3, 3]),
    Sigma_ref[2, 3] / sqrt(Sigma_ref[2, 2] * Sigma_ref[3, 3])
  )
  summ <- extract_armse_acp(reps, Sigma_ref, rho_ref)
  out <- data.frame(
    G = G_use,
    re_dist = re_dist_g,
    n_success = n_ok,
    n_sim = n_sim_g,
    ARMSE = summ$ARMSE,
    ACP = summ$ACP,
    n_param = summ$n_param
  )
  cat(sprintf(
    "G=%d success=%d/%d ARMSE=%.4f ACP=%.4f\n",
    G_use, n_ok, n_sim_g, out$ARMSE, out$ACP
  ))
  out
}

rows <- vector("list", length(G_grid))
for (k in seq_along(G_grid)) {
  G_k <- G_grid[[k]]
  seed_base <- 720000L + as.integer(G_k) * 1000L +
    if (identical(re_dist_g, "mixture")) 500000L else 0L
  rows[[k]] <- run_one_G(G_k, seed_base)
  tab <- do.call(rbind, rows[seq_len(k)])
  utils::write.csv(tab, csv_file, row.names = FALSE)
  saveRDS(
    list(
      settings = list(
        scenario = scenario_g, re_dist = re_dist_g, n = n_g,
        n_sim = n_sim_g, chain = chain_g, burn = burn_g, thin = thin_g,
        G_grid = G_grid, exclude_delta = TRUE
      ),
      summary = tab
    ),
    meta_file
  )
}

doParallel::stopImplicitCluster()
tab <- do.call(rbind, rows)
utils::write.csv(tab, csv_file, row.names = FALSE)

############################################################
# Single-series plot for this RE_DIST (combine offline for 2-line figure)
############################################################
draw_one <- function(tab, file, re_lab) {
  tab <- tab[is.finite(tab$G) & is.finite(tab$ARMSE) & is.finite(tab$ACP), ]
  tab <- tab[order(tab$G), ]
  if (!nrow(tab)) {
    cat("Skip plot: no finite rows\n")
    return(invisible(NULL))
  }
  df <- rbind(
    data.frame(G = tab$G, value = tab$ARMSE, Metric = "ARMSE"),
    data.frame(G = tab$G, value = tab$ACP, Metric = "ACP")
  )
  df$Metric <- factor(df$Metric, levels = c("ARMSE", "ACP"))
  col <- if (identical(re_lab, "mixture")) "#7B3FA0" else "#1B9E91"

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    ref95 <- data.frame(Metric = factor("ACP", levels = c("ARMSE", "ACP")), yint = 0.95)
    p <- ggplot2::ggplot(df, ggplot2::aes(x = G, y = value)) +
      ggplot2::geom_hline(
        data = ref95, ggplot2::aes(yintercept = yint),
        color = "red", linewidth = 0.5
      ) +
      ggplot2::geom_line(color = col, linewidth = 0.7, linetype = if (identical(re_lab, "mixture")) "dashed" else "solid") +
      ggplot2::geom_point(color = col, size = 2.4,
                          shape = if (identical(re_lab, "mixture")) 17 else 16) +
      ggplot2::facet_wrap(~Metric, ncol = 1, scales = "free_y") +
      ggplot2::scale_x_continuous(breaks = tab$G, labels = tab$G) +
      ggplot2::labs(
        x = expression(Truncation~level~G), y = NULL,
        title = sprintf("G sensitivity (%s) — Scenario 2, N=%d", re_lab, n_g)
      ) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(
        strip.background = ggplot2::element_rect(fill = "grey85", color = "grey40"),
        strip.text = ggplot2::element_text(face = "bold", hjust = 0.02),
        panel.grid.minor = ggplot2::element_blank()
      ) +
      ggplot2::geom_blank(data = data.frame(
        Metric = factor("ACP", levels = c("ARMSE", "ACP")),
        G = min(tab$G), value = 1.0
      ))
    ggplot2::ggsave(file, p, width = 7.2, height = 6.5, device = grDevices::pdf)
  } else {
    grDevices::pdf(file, width = 7.2, height = 6.5)
    graphics::par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))
    plot(tab$G, tab$ARMSE, type = "b", xlab = "G", ylab = "ARMSE", main = re_lab)
    plot(tab$G, tab$ACP, type = "b", xlab = "G", ylab = "ACP", ylim = c(0.9, 1))
    graphics::abline(h = 0.95, col = "red")
    grDevices::dev.off()
  }
  invisible(NULL)
}

draw_one(tab, pdf_file, re_dist_g)
cat("Wrote:\n  ", csv_file, "\n  ", pdf_file, "\n  ", meta_file, "\n", sep = "")
