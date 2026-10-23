# Combine HRS PPC into one (a)(b)(c) PDF from saved summaries.
# (a) count margin; (b) cor(HSPNIT, SHLT); (c) E(HSPNIT | SHLT).
args <- commandArgs(trailingOnly = TRUE)
OUT_DIR <- if (length(args) >= 1) args[[1]] else {
  code_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)
  parent <- dirname(code_dir)
  base <- basename(code_dir)
  sibling <- if (grepl("[\u4e00-\u9fff]", base)) {
    file.path(parent, paste0(base, "结果"))
  } else {
    file.path(parent, paste0(base, "_results"))
  }
  sibling
}
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

proj <- if (length(args) >= 2) args[[2]] else {
  dirname(normalizePath(".", winslash = "/", mustWork = TRUE))
}
ab_rds <- if (length(args) >= 3) args[[3]] else {
  file.path(proj, "hpc_jobs", "hrs_ppc_g8k23_xh5_20260915", "hrs_ppc_reps_summary.rds")
}
mean_rds <- if (length(args) >= 4) args[[4]] else {
  file.path(OUT_DIR, "hrs_ppc_mean_y1_by_y2.rds")
}
if (!file.exists(ab_rds)) stop("AB RDS not found: ", ab_rds)
if (!file.exists(mean_rds)) stop("Mean RDS not found: ", mean_rds)

x <- readRDS(ab_rds)
m <- readRDS(mean_rds)$summary

ks <- as.integer(x$ks)
obs <- as.numeric(x$obs_prop)
pred <- as.numeric(x$pred_prop)
lo <- as.numeric(x$lo_prop)
hi <- as.numeric(x$hi_prop)
k_cut <- ks[length(ks)]
idx_head <- seq_len(length(ks) - 1L)
obs_plot <- c(obs[idx_head], 1 - sum(obs[idx_head]))
pred_plot <- c(pred[idx_head], 1 - sum(pred[idx_head]))
lo_head <- lo[idx_head]
hi_head <- hi[idx_head]
lab_a <- c(as.character(ks[idx_head]), paste0(k_cut, "+"))

obs_mean <- as.numeric(m$Obs)
pred_mean <- as.numeric(m$Pred)
lo_mean <- as.numeric(m$Lo)
hi_mean <- as.numeric(m$Hi)
lab_c <- as.character(m$SHLT)

col_obs  <- "#0072B2"
col_pred <- "#E69F00"

out_pdf <- file.path(OUT_DIR, "hrs_ppc_abc.pdf")
# Horizontal (a)(b)(c): give (a) a bit more width for the many count bins.
grDevices::pdf(out_pdf, width = 11.5, height = 3.6)
op <- graphics::par(no.readonly = TRUE)
graphics::layout(matrix(1:3, nrow = 1L), widths = c(1.35, 1, 1))
graphics::par(mar = c(4.6, 3.8, 2.0, 0.7), cex.axis = 0.85, cex.lab = 0.95)

ylim_a <- max(c(obs_plot, hi_head), na.rm = TRUE) * 1.08
bp_a <- graphics::barplot(
  obs_plot, names.arg = lab_a, col = col_obs, border = NA,
  ylim = c(0, ylim_a), xlab = "HSPNIT count", ylab = "Proportion", main = "(a)",
  cex.names = 0.55, las = 2, mgp = c(2.4, 0.55, 0)
)
graphics::arrows(
  x0 = bp_a[idx_head], y0 = lo_head, y1 = hi_head,
  angle = 90, code = 3, length = 0.03, col = col_pred
)
graphics::points(bp_a, pred_plot, pch = 16, col = col_pred, cex = 0.7)

graphics::hist(
  x$rep_cor_y1_y2, breaks = 30, col = col_pred, border = "white",
  main = "(b)", xlab = expression(cor(HSPNIT, SHLT)), ylab = "Count",
  mgp = c(2.4, 0.7, 0)
)
graphics::abline(v = x$obs_cor_y1_y2, col = col_obs, lwd = 2)

ylim_c <- max(c(obs_mean, hi_mean), na.rm = TRUE) * 1.12
bp_c <- graphics::barplot(
  obs_mean, names.arg = lab_c, col = col_obs, border = NA,
  ylim = c(0, ylim_c),
  xlab = "SHLT",
  ylab = expression(Mean ~ of ~ HSPNIT),
  main = "(c)",
  mgp = c(2.4, 0.7, 0)
)
graphics::arrows(
  x0 = bp_c, y0 = lo_mean, y1 = hi_mean,
  angle = 90, code = 3, length = 0.05, col = col_pred, lwd = 1.5
)
graphics::points(bp_c, pred_mean, pch = 16, col = col_pred, cex = 1.1)

graphics::par(op)
grDevices::dev.off()
cat("Wrote", out_pdf, "\n")
