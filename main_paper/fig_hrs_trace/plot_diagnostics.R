#############################################
# ZINB + 有序（CDPMM）实际数据：收敛诊断可视化
# 布局对齐 二值-ZINB项目/实际数据分析对比单独模型/realdataplot.r
# 用法：
#   Rscript plot_diagnostics.R [RDS路径或含RDS的目录] [OUT_DIR]
#   无参数时自动查找 ../实际数据分析结果/ 下所有 MCMC_Analysis_Results_*.rds
#   多个 RDS 时分别出图（文件名带 RDS stem 前缀，避免覆盖）
#############################################

library(ggplot2)
library(coda)

for (pkg in c("patchwork", "cowplot", "gridExtra")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("请安装 ", pkg, ": install.packages('", pkg, "')")
  }
}

cli_args <- commandArgs(trailingOnly = TRUE)

default_sidecar <- function() {
  code_dir <- normalizePath(".", winslash = "/", mustWork = FALSE)
  parent <- dirname(code_dir)
  base <- basename(code_dir)
  if (grepl("[\u4e00-\u9fff]", base)) {
    file.path(parent, paste0(base, "结果"))
  } else {
    file.path(parent, paste0(base, "_results"))
  }
}

collect_rds_paths <- function(root = NULL) {
  sidecar <- default_sidecar()
  search_dirs <- character(0)
  if (!is.null(root) && nzchar(root)) {
    root <- normalizePath(root, winslash = "/", mustWork = FALSE)
    if (file.exists(root) && !dir.exists(root) && grepl("\\.rds$", root, ignore.case = TRUE)) {
      return(root)
    }
    if (dir.exists(root)) search_dirs <- c(search_dirs, root)
  }
  code_dir <- normalizePath(".", winslash = "/", mustWork = FALSE)
  search_dirs <- unique(c(search_dirs, sidecar, code_dir, "."))
  # Also scan sibling hpc_jobs hrs_* folders (downloaded job bundles)
  parent <- dirname(code_dir)
  hpc <- file.path(parent, "hpc_jobs")
  if (dir.exists(hpc)) {
    hrs <- list.dirs(hpc, full.names = TRUE, recursive = FALSE)
    hrs <- hrs[grepl("hrs_prog", basename(hrs))]
    search_dirs <- c(search_dirs, hrs)
  }

  prefer_names <- c(
    "MCMC_Analysis_Results_prog1_joint_vs_sep_cdpmm.rds",
    "MCMC_Analysis_Results_prog2_cdpmm_vs_gauss.rds"
  )
  hits <- character(0)
  for (d in search_dirs) {
    if (!dir.exists(d)) next
    found <- list.files(d, pattern = "^MCMC_Analysis_Results_.*\\.rds$", full.names = TRUE)
    hits <- c(hits, found)
  }
  hits <- unique(normalizePath(hits, winslash = "/", mustWork = FALSE))
  if (!length(hits)) return(character(0))

  # Prefer prog1/prog2 formal outputs; ignore old smoke-test RDS unless nothing else
  ordered <- character(0)
  for (nm in prefer_names) {
    match <- hits[basename(hits) == nm]
    if (length(match)) ordered <- c(ordered, match)
  }
  if (length(ordered)) return(unique(ordered))

  # Fallback: any MCMC_Analysis_Results_*.rds except tiny smoke tags if possible
  rest <- hits[order(file.info(hits)$mtime, decreasing = TRUE)]
  rest
}

rds_paths <- if (length(cli_args) >= 1L) {
  collect_rds_paths(cli_args[[1]])
} else {
  collect_rds_paths(NULL)
}
if (!length(rds_paths)) {
  stop("找不到 MCMC_Analysis_Results_*.rds（请先从超算下载 prog1/prog2 结果到 ../实际数据分析结果/）")
}

out_dir_arg <- if (length(cli_args) >= 2L) cli_args[[2]] else NA_character_

############################################################
# 工具函数（对齐 realdataplot.r）
############################################################

is_intercept_cov <- function(cov_name) {
  is.na(cov_name) || !nzchar(cov_name) ||
    identical(cov_name, "1") ||
    identical(cov_name, "Intercept") ||
    identical(cov_name, "(Intercept)") ||
    grepl("Intercept$", cov_name)
}

# Trace 标题：plotmath，如 alpha[SMOKEV]、beta[Intercept]
make_coef_label <- function(greek, cov_name) {
  if (is_intercept_cov(cov_name)) {
    sprintf("%s[Intercept]", greek)
  } else {
    sprintf("%s[%s]", greek, cov_name)
  }
}

# Geweke / EPSR 横轴：参数名+变量名，如 alpha(SMOKEV)
make_coef_display <- function(greek, cov_name) {
  if (is_intercept_cov(cov_name)) {
    sprintf("%s(Intercept)", greek)
  } else {
    sprintf("%s(%s)", greek, cov_name)
  }
}

parsed_plotmath_title <- function(label) {
  expr <- tryCatch(parse(text = label), error = function(e) NULL)
  if (length(expr) >= 1L) expr[[1L]] else label
}

extract_param_matrix <- function(chains_results, param_type) {
  if (param_type == "Alpha") return(lapply(chains_results, function(x) x$Alpha))
  if (param_type == "Beta")  return(lapply(chains_results, function(x) x$Beta))
  if (param_type == "Gamma") return(lapply(chains_results, function(x) x$Gamma))
  if (param_type == "Delta") return(lapply(chains_results, function(x) x$Delta))
  if (param_type == "R")     return(lapply(chains_results, function(x) matrix(x$R, ncol = 1)))
  if (param_type == "Tau")   return(lapply(chains_results, function(x) matrix(x$Tau, ncol = 1)))
  if (param_type == "Nclust") return(lapply(chains_results, function(x) matrix(x$Nclust, ncol = 1)))
  if (param_type == "Sigma") {
    return(lapply(chains_results, function(x) {
      n_s <- nrow(x$Sigma)
      out <- matrix(0, n_s, 3)
      for (i in seq_len(n_s)) {
        S <- matrix(x$Sigma[i, ], 3, 3)
        out[i, ] <- c(S[1, 1], S[2, 2], S[3, 3])
      }
      out
    }))
  }
  if (param_type == "Rho") {
    return(lapply(chains_results, function(x) {
      if (!is.null(x$Rho)) return(x$Rho)
      n_s <- nrow(x$Sigma)
      out <- matrix(0, n_s, 3)
      for (i in seq_len(n_s)) {
        S <- matrix(x$Sigma[i, ], 3, 3)
        out[i, 1] <- S[1, 2] / sqrt(S[1, 1] * S[2, 2])
        out[i, 2] <- S[1, 3] / sqrt(S[1, 1] * S[3, 3])
        out[i, 3] <- S[2, 3] / sqrt(S[2, 2] * S[3, 3])
      }
      out
    }))
  }
  stop("Unknown param_type: ", param_type)
}

plot_trace_panel <- function(combined_data, ylab = "Value") {
  ggplot(combined_data, aes(x = iteration, y = value, color = chain)) +
    geom_line(linewidth = 0.45, alpha = 0.75) +
    labs(x = NULL, y = ylab, color = "Chain") +
    theme_bw(base_size = 8) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 7),
      axis.text.y = element_text(size = 6),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      plot.margin = margin(2, 3, 2, 2, "pt")
    ) +
    scale_color_manual(values = c("Chain 1" = "blue", "Chain 2" = "green", "Chain 3" = "red"))
}

plot_density_panel <- function(combined_data, ylab = "Density") {
  ggplot(combined_data, aes(x = value, color = chain)) +
    geom_density(linewidth = 0.45, alpha = 0.75, fill = NA) +
    labs(x = NULL, y = ylab, color = "Chain") +
    theme_bw(base_size = 8) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 7),
      axis.text.y = element_text(size = 6),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      plot.margin = margin(2, 3, 2, 2, "pt")
    ) +
    scale_color_manual(values = c("Chain 1" = "blue", "Chain 2" = "green", "Chain 3" = "red"))
}

build_param_combined_data <- function(chains_results, param_type, param_idx) {
  mats <- extract_param_matrix(chains_results, param_type)
  out <- data.frame()
  for (chain in seq_along(mats)) {
    vals <- if (ncol(mats[[chain]]) == 1L) mats[[chain]][, 1] else mats[[chain]][, param_idx]
    out <- rbind(out, data.frame(
      iteration = seq_along(vals),
      value = vals,
      chain = factor(paste0("Chain ", chain))
    ))
  }
  out
}

build_mcmc_grid <- function(chains_results, param_specs, grid_cols = 5L, grid_rows = NULL) {
  n_params <- length(param_specs)
  if (is.null(grid_rows)) grid_rows <- ceiling(n_params / grid_cols)
  if (n_params > grid_cols * grid_rows) {
    stop("参数个数 (", n_params, ") 超过 ", grid_rows, "×", grid_cols, " 网格容量")
  }

  trace_plots <- vector("list", n_params)
  density_plots <- vector("list", n_params)
  for (k in seq_len(n_params)) {
    spec <- param_specs[[k]]
    dat <- build_param_combined_data(chains_results, spec$type, spec$idx)
    trace_plots[[k]] <- plot_trace_panel(dat, ylab = if (k == 1L) "Value" else NULL) +
      labs(title = parsed_plotmath_title(spec$label)) +
      theme(plot.title = element_text(size = 8, hjust = 0.5, face = "plain"))
    density_plots[[k]] <- plot_density_panel(dat, ylab = if (k == 1L) "Density" else NULL)
  }

  row_blocks <- list()
  for (r in seq_len(grid_rows)) {
    idx_start <- (r - 1L) * grid_cols + 1L
    idx_end <- min(r * grid_cols, n_params)
    if (idx_start > n_params) break
    idx <- idx_start:idx_end
    n_in_row <- length(idx)

    if (n_in_row < grid_cols) {
      blank <- ggplot() + theme_void()
      blanks <- rep(list(blank), grid_cols - n_in_row)
      trace_row <- patchwork::wrap_plots(c(trace_plots[idx], blanks), ncol = grid_cols)
      dens_row <- patchwork::wrap_plots(c(density_plots[idx], blanks), ncol = grid_cols)
    } else {
      trace_row <- patchwork::wrap_plots(trace_plots[idx], ncol = grid_cols)
      dens_row <- patchwork::wrap_plots(density_plots[idx], ncol = grid_cols)
    }

    label_trace <- cowplot::ggdraw() +
      cowplot::draw_label("Trace", x = 0.5, y = 0.5, angle = 90, size = 9, fontface = "bold")
    label_dens <- cowplot::ggdraw() +
      cowplot::draw_label("Density", x = 0.5, y = 0.5, angle = 90, size = 9, fontface = "bold")

    row_blocks[[2L * r - 1L]] <- cowplot::plot_grid(
      label_trace, trace_row, ncol = 2, rel_widths = c(0.04, 1)
    )
    row_blocks[[2L * r]] <- cowplot::plot_grid(
      label_dens, dens_row, ncol = 2, rel_widths = c(0.04, 1)
    )
  }
  patchwork::wrap_plots(row_blocks, ncol = 1)
}

rhat_safe <- function(samples_list) {
  if (length(samples_list) < 2L) return(NA_real_)
  mcmc_list <- tryCatch(
    coda::as.mcmc.list(lapply(samples_list, coda::as.mcmc)),
    error = function(e) NULL
  )
  if (is.null(mcmc_list)) return(NA_real_)
  gd <- tryCatch(
    coda::gelman.diag(mcmc_list, autoburnin = FALSE, multivariate = FALSE),
    error = function(e) NULL
  )
  if (is.null(gd)) return(NA_real_)
  as.numeric(gd$psrf[1, 1])
}

cov_from_param <- function(nm, prefixes) {
  for (prefix in prefixes) {
    if (startsWith(nm, prefix)) return(sub(paste0("^", prefix), "", nm))
  }
  nm
}

param_names_from_summary <- function(df, fallback_n, prefix) {
  if (!is.null(df) && "Parameter" %in% names(df)) return(as.character(df$Parameter))
  paste0(prefix, seq_len(fallback_n) - 1L)
}

run_one_diagnostics <- function(rds_path, out_dir) {
  cat("\n========== 读取结果:", rds_path, "==========\n")
  res <- readRDS(rds_path)
  chains_results <- res$chains_results
  if (is.null(chains_results) || !length(chains_results)) {
    stop("RDS 中无 chains_results: ", rds_path)
  }
  n_chains <- length(chains_results)
  stem <- sub("\\.rds$", "", basename(rds_path), ignore.case = TRUE)
  # Short tag for filenames
  tag <- if (grepl("prog1", stem, ignore.case = TRUE)) {
    "prog1"
  } else if (grepl("prog2", stem, ignore.case = TRUE)) {
    "prog2"
  } else {
    stem
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # checkpoint_*.rds 通常无 summary；HRS prog1/2 固定效应为 Intercept + 5 协变量
  hrs_covs <- c("Intercept", "SMOKEV", "HIBP", "DIAB", "LUNG", "BMI")
  fallback_coef_names <- function(prefix, n) {
    if (identical(as.integer(n), 6L)) {
      paste0(prefix, hrs_covs)
    } else {
      paste0(prefix, seq_len(n) - 1L)
    }
  }

  alpha_names <- if (!is.null(res$summary$alpha) && "Parameter" %in% names(res$summary$alpha)) {
    as.character(res$summary$alpha$Parameter)
  } else {
    fallback_coef_names("alpha_zero_", ncol(chains_results[[1]]$Alpha))
  }
  beta_names <- if (!is.null(res$summary$beta) && "Parameter" %in% names(res$summary$beta)) {
    as.character(res$summary$beta$Parameter)
  } else {
    fallback_coef_names("beta_count_", ncol(chains_results[[1]]$Beta))
  }
  gamma_names <- if (!is.null(res$summary$gamma) && "Parameter" %in% names(res$summary$gamma)) {
    as.character(res$summary$gamma$Parameter)
  } else {
    fallback_coef_names("gamma_ord_", ncol(chains_results[[1]]$Gamma))
  }
  delta_names <- if (!is.null(res$summary$delta) && "Parameter" %in% names(res$summary$delta)) {
    as.character(res$summary$delta$Parameter)
  } else {
    paste0("delta", seq_len(ncol(chains_results[[1]]$Delta)))
  }
  rho_names <- if (!is.null(res$summary$rho) && "Parameter" %in% names(res$summary$rho)) {
    as.character(res$summary$rho$Parameter)
  } else {
    paste0("rho", seq_len(ncol(chains_results[[1]]$Rho)))
  }
  sigma_names <- c("Sigma[1,1]", "Sigma[2,2]", "Sigma[3,3]")

  alpha_cov <- vapply(alpha_names, cov_from_param, character(1),
                      prefixes = c("alpha_ord_", "alpha_zero_", "alpha"))
  beta_cov  <- vapply(beta_names,  cov_from_param, character(1),
                      prefixes = c("beta_zero_", "beta_count_", "beta"))
  gamma_cov <- vapply(gamma_names, cov_from_param, character(1),
                      prefixes = c("gamma_count_", "gamma_ord_", "gamma"))

  cat("=== ZINB-有序 CDPMM 诊断图 (", tag, ") ===\n")
  cat("链数:", n_chains, " 每链样本:", nrow(chains_results[[1]]$Alpha), "\n")

  ############################################################
  # 1. Trace / Density（固定效应 α/β/γ + 离散参数 r）
  ############################################################
  cat("\n1. 创建 Trace/Density 诊断图...\n")

  param_specs <- list()
  for (i in seq_along(alpha_names)) {
    param_specs[[length(param_specs) + 1L]] <- list(
      type = "Alpha", idx = i, name = alpha_names[i],
      label = make_coef_label("alpha", alpha_cov[i]),
      display = make_coef_display("alpha", alpha_cov[i]),
      family = "Alpha (zero)"
    )
  }
  for (i in seq_along(beta_names)) {
    param_specs[[length(param_specs) + 1L]] <- list(
      type = "Beta", idx = i, name = beta_names[i],
      label = make_coef_label("beta", beta_cov[i]),
      display = make_coef_display("beta", beta_cov[i]),
      family = "Beta (count)"
    )
  }
  for (i in seq_along(gamma_names)) {
    param_specs[[length(param_specs) + 1L]] <- list(
      type = "Gamma", idx = i, name = gamma_names[i],
      label = make_coef_label("gamma", gamma_cov[i]),
      display = make_coef_display("gamma", gamma_cov[i]),
      family = "Gamma (ordinal)"
    )
  }
  param_specs[[length(param_specs) + 1L]] <- list(
    type = "R", idx = 1L, name = "r", label = "r",
    display = "r", family = "Dispersion"
  )

  grid_cols <- 4L
  grid_rows <- ceiling(length(param_specs) / grid_cols)
  cat(sprintf("  共 %d 个参数 (α/β/γ/r) -> %d 行 × %d 列 (Trace+Density 共 %d 行)\n",
              length(param_specs), grid_rows, grid_cols, 2L * grid_rows))

  mcmc_grid <- build_mcmc_grid(chains_results, param_specs, grid_cols, grid_rows)
  fig_h <- max(8, 2.0 * grid_rows * 2)
  f_trace <- file.path(out_dir, paste0(tag, "_MCMC_Comprehensive_Convergence_Diagnostics.pdf"))
  ggplot2::ggsave(
    f_trace, plot = mcmc_grid, width = 14, height = fig_h,
    device = grDevices::pdf, limitsize = FALSE
  )
  cat("已保存:", f_trace, "\n")

  ############################################################
  # 1b. Trace / Density：随机效应方差对角 + 相关系数
  ############################################################
  cat("\n1b. 创建 Sigma/Rho Trace/Density 诊断图...\n")
  sigma_rho_specs <- list()
  sigma_labels <- c("Sigma[1,1]", "Sigma[2,2]", "Sigma[3,3]")
  sigma_displays <- c("Sigma11", "Sigma22", "Sigma33")
  for (i in seq_along(sigma_labels)) {
    sigma_rho_specs[[length(sigma_rho_specs) + 1L]] <- list(
      type = "Sigma", idx = i, name = sigma_names[i],
      label = sigma_labels[i], display = sigma_displays[i],
      family = "Sigma (diag)"
    )
  }
  rho_labels <- c("rho[12]", "rho[13]", "rho[23]")
  n_rho <- if (!is.null(chains_results[[1]]$Rho)) ncol(chains_results[[1]]$Rho) else 3L
  for (i in seq_len(min(3L, n_rho))) {
    nm <- if (i <= length(rho_names)) rho_names[i] else paste0("rho", i)
    sigma_rho_specs[[length(sigma_rho_specs) + 1L]] <- list(
      type = "Rho", idx = i, name = nm,
      label = rho_labels[i], display = nm,
      family = "Rho"
    )
  }
  sr_cols <- 3L
  sr_rows <- ceiling(length(sigma_rho_specs) / sr_cols)
  cat(sprintf("  共 %d 个参数 (Sigma diag + rho) -> %d 行 × %d 列\n",
              length(sigma_rho_specs), sr_rows, sr_cols))
  sigma_rho_grid <- build_mcmc_grid(chains_results, sigma_rho_specs, sr_cols, sr_rows)
  f_sigma_rho <- file.path(out_dir, paste0(tag, "_MCMC_Sigma_Rho_Convergence_Diagnostics.pdf"))
  ggplot2::ggsave(
    f_sigma_rho, plot = sigma_rho_grid, width = 12, height = max(6, 2.0 * sr_rows * 2),
    device = grDevices::pdf, limitsize = FALSE
  )
  cat("已保存:", f_sigma_rho, "\n")

  ############################################################
  # 2. Geweke
  ############################################################
  cat("\n2. 创建 Geweke 诊断图...\n")
  geweke_df <- data.frame()
  for (chain in seq_len(n_chains)) {
    for (spec in param_specs) {
      mats <- extract_param_matrix(list(chains_results[[chain]]), spec$type)[[1]]
      vals <- if (ncol(mats) == 1L) mats[, 1] else mats[, spec$idx]
      if (length(vals) < 10L) next
      z <- tryCatch(as.numeric(coda::geweke.diag(coda::as.mcmc(vals))$z),
                    error = function(e) NA_real_)
      geweke_df <- rbind(geweke_df, data.frame(
        chain = paste0("Chain ", chain),
        display_name = spec$display,
        z = z,
        family = spec$family,
        stringsAsFactors = FALSE
      ))
    }
  }

  if (nrow(geweke_df) > 0) {
    geweke_df$display_name <- factor(
      geweke_df$display_name,
      levels = unique(vapply(param_specs, `[[`, "", "display"))
    )
    geweke_plot <- ggplot(geweke_df, aes(x = display_name, y = z, color = chain, shape = chain)) +
      geom_hline(yintercept = c(-2, 2), linetype = "dashed", color = "grey50") +
      geom_point(size = 2, position = position_dodge(width = 0.5)) +
      labs(x = NULL, y = "Geweke z", color = "Chain", shape = "Chain",
           title = paste0("Geweke diagnostic (", tag, ")")) +
      theme_bw(base_size = 10) +
      theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7),
            legend.position = "bottom") +
      scale_color_manual(values = c("Chain 1" = "blue", "Chain 2" = "green", "Chain 3" = "red"))
  } else {
    geweke_plot <- ggplot() + theme_void() +
      labs(title = "Geweke skipped (too few posterior draws)")
  }

  ############################################################
  # 3. EPSR / R-hat
  ############################################################
  cat("\n3. 创建 EPSR (R-hat) 图...\n")
  epsr_df <- data.frame()
  if (n_chains >= 2L) {
    for (spec in param_specs) {
      mats <- extract_param_matrix(chains_results, spec$type)
      samples_list <- lapply(mats, function(m) {
        if (ncol(m) == 1L) m[, 1] else m[, spec$idx]
      })
      rh <- rhat_safe(samples_list)
      epsr_df <- rbind(epsr_df, data.frame(
        parameter = spec$display,
        R_hat = rh,
        type = spec$family,
        stringsAsFactors = FALSE
      ))
    }
    epsr_plot <- ggplot(epsr_df, aes(x = parameter, y = R_hat, color = type, shape = type)) +
      geom_hline(yintercept = 1.1, linetype = "dashed", color = "orange") +
      geom_hline(yintercept = 1.2, linetype = "dashed", color = "red") +
      geom_point(size = 2.5) +
      labs(x = NULL, y = "EPSR (R-hat)", color = "Type", shape = "Type",
           title = paste0("Gelman–Rubin EPSR (", tag, ")")) +
      theme_bw(base_size = 10) +
      theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7),
            legend.position = "bottom") +
      ylim(0.8, max(2.0, max(epsr_df$R_hat, na.rm = TRUE) + 0.2, na.rm = TRUE))
  } else {
    epsr_df <- data.frame(parameter = character(), R_hat = numeric(), type = character())
    epsr_plot <- ggplot() + theme_void() +
      labs(title = paste0("EPSR requires >= 2 chains (n_chains = ", n_chains, ")"))
  }

  ############################################################
  # 4. Save Geweke + EPSR
  ############################################################
  diag_combo <- cowplot::plot_grid(geweke_plot, epsr_plot, ncol = 1, rel_heights = c(1, 1))
  f_epsr <- file.path(out_dir, paste0(tag, "_MCMC_Geweke_EPSR_Diagnostics.pdf"))
  ggplot2::ggsave(f_epsr, plot = diag_combo, width = 12, height = 10, device = grDevices::pdf)
  cat("已保存:", f_epsr, "\n")

  if (nrow(epsr_df) > 0) {
    f_csv <- file.path(out_dir, paste0(tag, "_EPSR_Joint_ZINB_Ordinal_RealData.csv"))
    write.csv(epsr_df, f_csv, row.names = FALSE)
    cat("已保存:", f_csv, "\n")
  }
  if (nrow(geweke_df) > 0) {
    f_csv <- file.path(out_dir, paste0(tag, "_Geweke_Joint_ZINB_Ordinal_RealData.csv"))
    write.csv(geweke_df, f_csv, row.names = FALSE)
    cat("已保存:", f_csv, "\n")
  }

  if (nrow(epsr_df) > 0 && any(is.finite(epsr_df$R_hat))) {
    cat("\n=== EPSR 摘要 (", tag, ") ===\n")
    cat("  良好 (<1.1):", sum(epsr_df$R_hat < 1.1, na.rm = TRUE), "\n")
    cat("  可接受 [1.1,1.2):", sum(epsr_df$R_hat >= 1.1 & epsr_df$R_hat < 1.2, na.rm = TRUE), "\n")
    cat("  需关注 (>=1.2):", sum(epsr_df$R_hat >= 1.2, na.rm = TRUE), "\n")
    bad <- epsr_df[!is.na(epsr_df$R_hat) & epsr_df$R_hat >= 1.2, ]
    if (nrow(bad) > 0) {
      cat("  R-hat>=1.2 参数:\n")
      print(bad, row.names = FALSE)
    }
  }

  cat("\n=== 诊断图完成 (", tag, ") ===\n")
  invisible(list(tag = tag, trace = f_trace, epsr = f_epsr))
}

sidecar <- default_sidecar()
cat("将处理", length(rds_paths), "个 RDS:\n")
print(rds_paths)

for (rds_path in rds_paths) {
  out_dir <- if (!is.na(out_dir_arg) && nzchar(out_dir_arg)) {
    out_dir_arg
  } else {
    sidecar
  }
  run_one_diagnostics(rds_path, out_dir)
}

cat("\n全部诊断图已写入:", if (!is.na(out_dir_arg) && nzchar(out_dir_arg)) out_dir_arg else sidecar, "\n")
