############################################################
# Fixed-effect prior sensitivity (Tang-style TYPE I/II/III x kappa)
# Defaults: Scenario 2, N=100, re_dist=mixture, S=200, G=8
# Prior: alpha,beta,gamma ~ N(mean, kappa * I); kappa = prior variance
# ARMSE/ACP over alpha,beta,gamma,r, Sigma diag, rho (exclude delta)
# Env: N, N_SIM, CHAIN, BURN, THIN, SCENARIO, RE_DIST, PRIOR_TYPE,
#      KAPPA_GRID, G_MIX, OUT_DIR
############################################################

args <- commandArgs(trailingOnly = TRUE)
OUT_DIR <- if (length(args) >= 1L && nzchar(args[[1]])) {
  args[[1]]
} else {
  env_out <- Sys.getenv("OUT_DIR", unset = "")
  if (nzchar(env_out)) env_out else {
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

n_sim_ps <- as.integer(Sys.getenv("N_SIM", unset = "200"))
chain_ps <- as.integer(Sys.getenv("CHAIN", unset = "5000"))
burn_ps <- as.integer(Sys.getenv("BURN", unset = "2000"))
thin_ps <- as.integer(Sys.getenv("THIN", unset = "5"))
n_ps <- as.integer(Sys.getenv("N", unset = "100"))
scenario_ps <- as.integer(Sys.getenv("SCENARIO", unset = "2"))
re_dist_ps <- tolower(Sys.getenv("RE_DIST", unset = "mixture"))
G_ps <- as.integer(Sys.getenv("G_MIX", unset = "8"))
if (!nzchar(re_dist_ps)) re_dist_ps <- "mixture"
if (!re_dist_ps %in% c("normal", "mixture")) {
  stop("RE_DIST must be normal or mixture")
}
if (!is.finite(n_ps) || n_ps < 1L) n_ps <- 100L
if (!is.finite(n_sim_ps) || n_sim_ps < 1L) n_sim_ps <- 200L
if (!is.finite(G_ps) || G_ps < 2L) G_ps <- 8L

# True FE: Scenario 1 is 5/6/4; Scenario 2 (no I(U>0)) is 4/5/3
fe_true <- true_fixed_effects(scenario_ps)
alpha_true <- fe_true$alpha
beta_true <- fe_true$beta
gamma_true <- fe_true$gamma
r_true <- 2
cat(sprintf(
  "FE dims: alpha=%d beta=%d gamma=%d (scenario %d)\n",
  length(alpha_true), length(beta_true), length(gamma_true), scenario_ps
))

prior_type_env <- toupper(Sys.getenv("PRIOR_TYPE", unset = "ALL"))
# Accept I/II/III, 1/2/3, TYPE I, ALL
normalize_type <- function(x) {
  x <- gsub("^TYPE[[:space:]]*", "", x)
  x <- gsub("[[:space:]]+", "", x)
  if (x %in% c("1", "I")) return("TYPE I")
  if (x %in% c("2", "II")) return("TYPE II")
  if (x %in% c("3", "III")) return("TYPE III")
  if (x %in% c("ALL", "")) return("ALL")
  stop("PRIOR_TYPE must be I/II/III or ALL, got: ", x)
}
prior_type_req <- normalize_type(prior_type_env)
type_grid <- if (identical(prior_type_req, "ALL")) {
  c("TYPE I", "TYPE II", "TYPE III")
} else {
  prior_type_req
}

k_env <- Sys.getenv("KAPPA_GRID", unset = "")
kappa_grid <- if (nzchar(k_env)) {
  as.numeric(strsplit(k_env, "[,;[:space:]]+")[[1]])
} else {
  c(0.1, 0.5, 1, 2, 5, 10, 50, 100, 1000)
}
kappa_grid <- kappa_grid[is.finite(kappa_grid) & kappa_grid > 0]
if (!length(kappa_grid)) stop("Empty KAPPA_GRID")

prior_means_for <- function(type_lab) {
  if (identical(type_lab, "TYPE I")) {
    list(a = alpha_true, b = beta_true, g = gamma_true)
  } else if (identical(type_lab, "TYPE II")) {
    list(a = rep(0, length(alpha_true)),
         b = rep(0, length(beta_true)),
         g = rep(0, length(gamma_true)))
  } else if (identical(type_lab, "TYPE III")) {
    list(a = -alpha_true, b = -beta_true, g = -gamma_true)
  } else {
    stop("Unknown type: ", type_lab)
  }
}

n_cores <- suppressWarnings(as.integer(Sys.getenv(
  c("N_CORES", "SLURM_NTASKS_PER_NODE", "SLURM_CPUS_PER_TASK")
)))
n_cores <- n_cores[is.finite(n_cores) & n_cores >= 1L]
n_cores <- if (length(n_cores)) n_cores[[1]] else max(1L, parallel::detectCores() - 1L)
n_cores <- min(n_cores, n_sim_ps)
doParallel::registerDoParallel(cores = n_cores)

cat(sprintf(
  "Prior sensitivity | scenario=%d N=%d re_dist=%s S=%d G=%d types=%s\n",
  scenario_ps, n_ps, re_dist_ps, n_sim_ps, G_ps, paste(type_grid, collapse = ",")
))
cat("kappa grid:", paste(kappa_grid, collapse = ", "), "\n")
cat("OUT_DIR=", OUT_DIR, "\n", sep = "")

tag <- sprintf(
  "scen%d_%s_n%d_nsim%d",
  scenario_ps, re_dist_ps, n_ps, n_sim_ps
)
csv_file <- file.path(OUT_DIR, paste0("prior_sensitivity_", tag, ".csv"))
meta_file <- file.path(OUT_DIR, paste0("prior_sensitivity_", tag, "_meta.rds"))
pdf_file <- file.path(OUT_DIR, paste0("prior_sensitivity_", tag, ".pdf"))

extract_armse_acp <- function(fits) {
  n_ok <- length(fits)
  p_a <- length(alpha_true)
  p_b <- length(beta_true)
  p_g <- length(gamma_true)
  a_est <- matrix(NA_real_, n_ok, p_a); a_lo <- a_est; a_hi <- a_est
  b_est <- matrix(NA_real_, n_ok, p_b); b_lo <- b_est; b_hi <- b_est
  g_est <- matrix(NA_real_, n_ok, p_g); g_lo <- g_est; g_hi <- g_est
  r_est <- rep(NA_real_, n_ok); r_lo <- r_est; r_hi <- r_est
  s_est <- matrix(NA_real_, n_ok, 3L); s_lo <- s_est; s_hi <- s_est
  rho_est <- matrix(NA_real_, n_ok, 3L); rho_lo <- rho_est; rho_hi <- rho_est
  s_true <- matrix(NA_real_, n_ok, 3L); rho_true_mat <- matrix(NA_real_, n_ok, 3L)

  for (i in seq_len(n_ok)) {
    fit <- fits[[i]]$fit
    Sig_t <- fits[[i]]$Sigma_true
    a_est[i, ] <- fit$alpha_est; a_lo[i, ] <- fit$alpha_ci_lower; a_hi[i, ] <- fit$alpha_ci_upper
    b_est[i, ] <- fit$beta_est; b_lo[i, ] <- fit$beta_ci_lower; b_hi[i, ] <- fit$beta_ci_upper
    g_est[i, ] <- fit$gamma_est; g_lo[i, ] <- fit$gamma_ci_lower; g_hi[i, ] <- fit$gamma_ci_upper
    r_est[i] <- fit$r_est; r_lo[i] <- fit$r_ci_lower; r_hi[i] <- fit$r_ci_upper
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
  st_r <- calculate_stats(matrix(r_est, ncol = 1), r_true,
                          matrix(r_lo, ncol = 1), matrix(r_hi, ncol = 1))
  st_s <- calculate_stats_varying_true(s_est, s_true, s_lo, s_hi)
  st_rho <- calculate_stats_varying_true(rho_est, rho_true_mat, rho_lo, rho_hi)
  rmse_vec <- c(st_a$rmse, st_b$rmse, st_g$rmse, st_r$rmse, st_s$rmse, st_rho$rmse)
  cp_vec <- c(st_a$cp, st_b$cp, st_g$cp, st_r$cp, st_s$cp, st_rho$cp)
  list(ARMSE = mean(rmse_vec, na.rm = TRUE),
       ACP = mean(cp_vec, na.rm = TRUE),
       n_param = length(rmse_vec))
}

run_one_cell <- function(type_lab, kappa, seed_base) {
  means <- prior_means_for(type_lab)
  cat(sprintf("\n=== %s | kappa=%.4g | %d reps ===\n", type_lab, kappa, n_sim_ps))
  reps <- foreach::foreach(
    i = seq_len(n_sim_ps),
    .packages = c("BayesLogit", "mvtnorm", "MCMCpack", "truncnorm", "loo", "coda")
  ) %dopar% {
    tryCatch({
      dat <- generate_data(
        n = n_ps, nis = nis_fixed, seed = seed_base + i,
        random_nis = random_nis, nis_range = nis_range,
        scenario = scenario_ps, C = C,
        re_dist = re_dist_ps, mvt_df = mvt_df
      )
      fit <- fit_joint_model(
        dat, chain = chain_ps, burn = burn_ps, thin = thin_ps,
        delta_min = delta_min, delta_max = delta_max,
        G = G_ps,
        alpha0 = means$a, beta0 = means$b, gamma0 = means$g,
        prior_var = kappa
      )
      list(fit = fit, Sigma_true = dat$Sigma_true)
    }, error = function(e) list(error = conditionMessage(e)))
  }
  ok <- vapply(reps, function(z) is.null(z$error) && !is.null(z$fit), logical(1))
  if (any(!ok)) {
    errs <- vapply(reps[!ok], function(z) {
      if (!is.null(z$error)) z$error else "missing fit"
    }, character(1))
    cat("Failed:\n"); print(utils::head(unique(errs), 5L))
  }
  reps <- reps[ok]
  n_ok <- length(reps)
  if (n_ok < 1L) {
    return(data.frame(
      Type = type_lab, kappa = kappa, n_success = 0L, n_sim = n_sim_ps,
      ARMSE = NA_real_, ACP = NA_real_, n_param = NA_integer_
    ))
  }
  summ <- extract_armse_acp(reps)
  out <- data.frame(
    Type = type_lab, kappa = kappa,
    n_success = n_ok, n_sim = n_sim_ps,
    ARMSE = summ$ARMSE, ACP = summ$ACP, n_param = summ$n_param
  )
  cat(sprintf("%s kappa=%.4g success=%d/%d ARMSE=%.4f ACP=%.4f\n",
              type_lab, kappa, n_ok, n_sim_ps, out$ARMSE, out$ACP))
  out
}

rows <- list()
idx <- 0L
for (type_lab in type_grid) {
  type_code <- match(type_lab, c("TYPE I", "TYPE II", "TYPE III"))
  for (kappa in kappa_grid) {
    idx <- idx + 1L
    seed_base <- 830000L + type_code * 100000L + as.integer(round(kappa * 1000))
    rows[[idx]] <- run_one_cell(type_lab, kappa, seed_base)
    tab <- do.call(rbind, rows)
    utils::write.csv(tab, csv_file, row.names = FALSE)
    saveRDS(list(
      settings = list(
        scenario = scenario_ps, re_dist = re_dist_ps, n = n_ps,
        n_sim = n_sim_ps, chain = chain_ps, burn = burn_ps, thin = thin_ps,
        G = G_ps, type_grid = type_grid, kappa_grid = kappa_grid,
        exclude_delta = TRUE
      ),
      summary = tab
    ), meta_file)
  }
}

doParallel::stopImplicitCluster()
tab <- do.call(rbind, rows)
utils::write.csv(tab, csv_file, row.names = FALSE)

############################################################
# Tang-style plot (log-x kappa)
############################################################
draw_prior_sens_plot <- function(tab, file) {
  tab <- tab[is.finite(tab$kappa) & is.finite(tab$ARMSE) & is.finite(tab$ACP), ]
  if (!nrow(tab)) {
    cat("Skip plot: no finite rows\n")
    return(invisible(NULL))
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cat("ggplot2 not available; skip PDF\n")
    return(invisible(NULL))
  }
  tab$Type <- factor(tab$Type, levels = c("TYPE I", "TYPE II", "TYPE III"))
  df <- rbind(
    data.frame(hyper = tab$kappa, value = tab$ARMSE, Type = tab$Type, Metric = "ARMSE"),
    data.frame(hyper = tab$kappa, value = tab$ACP, Type = tab$Type, Metric = "ACP")
  )
  df$Metric <- factor(df$Metric, levels = c("ARMSE", "ACP"))
  hyper <- sort(unique(tab$kappa))
  ref95 <- data.frame(Metric = factor("ACP", levels = c("ARMSE", "ACP")), yint = 0.95)
  pal <- c("TYPE I" = "#1B9E91", "TYPE II" = "#7B3FA0", "TYPE III" = "#E69500")
  shapes <- c("TYPE I" = 16, "TYPE II" = 17, "TYPE III" = 15)
  ltypes <- c("TYPE I" = "solid", "TYPE II" = "dashed", "TYPE III" = "dotdash")

  p <- ggplot2::ggplot(df, ggplot2::aes(hyper, value, color = Type, shape = Type, linetype = Type)) +
    ggplot2::geom_hline(data = ref95, ggplot2::aes(yintercept = yint),
                        color = "red", linewidth = 0.5) +
    ggplot2::geom_line(linewidth = 0.6) +
    ggplot2::geom_point(size = 2.2, fill = "white", stroke = 0.5) +
    ggplot2::facet_wrap(~Metric, ncol = 1, scales = "free_y") +
    ggplot2::scale_x_log10(breaks = hyper, labels = hyper) +
    ggplot2::scale_color_manual(values = pal) +
    ggplot2::scale_shape_manual(values = shapes) +
    ggplot2::scale_linetype_manual(values = ltypes) +
    ggplot2::labs(x = "Prior hyperparameter", y = NULL,
                  color = NULL, shape = NULL, linetype = NULL) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      legend.position = "inside",
      legend.position.inside = c(0.82, 0.92),
      legend.background = ggplot2::element_rect(
        fill = scales::alpha("white", 0.7), color = NA
      ),
      legend.key.width = grid::unit(1.2, "cm"),
      strip.background = ggplot2::element_rect(fill = "grey85", color = "grey40"),
      strip.text = ggplot2::element_text(face = "bold", hjust = 0.02),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "grey92", linewidth = 0.3)
    ) +
    ggplot2::geom_blank(data = data.frame(
      Metric = factor("ACP", levels = c("ARMSE", "ACP")),
      hyper = min(hyper), value = 1.0,
      Type = factor("TYPE I", levels = c("TYPE I", "TYPE II", "TYPE III"))
    )) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.05))) +
    ggplot2::geom_text(
      data = data.frame(
        Metric = factor("ACP", levels = c("ARMSE", "ACP")),
        hyper = min(hyper) * 1.2, value = 0.95
      ),
      ggplot2::aes(x = hyper, y = value, label = "0.95"),
      inherit.aes = FALSE, color = "red", vjust = 1.6, hjust = 0, size = 3.2
    )
  ggplot2::ggsave(file, p, width = 6.5, height = 6.5, device = grDevices::pdf)
  invisible(NULL)
}

draw_prior_sens_plot(tab, pdf_file)
cat("Wrote:\n  ", csv_file, "\n  ", pdf_file, "\n  ", meta_file, "\n", sep = "")
