# Redraw HRS PPC panels (a) count margin + (b) cor(y1, y2) from saved RDS.
# Colors differ from Paper 1 teal/red; last bar pools counts >= 20.
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

rds_path <- if (length(args) >= 2) args[[2]] else {
  file.path(
    dirname(normalizePath(".", winslash = "/", mustWork = TRUE)),
    "hpc_jobs",
    "hrs_ppc_g8k23_xh5_20260915",
    "hrs_ppc_reps_summary.rds"
  )
}
if (!file.exists(rds_path)) {
  stop("RDS not found: ", rds_path)
}
x <- readRDS(rds_path)

ks <- as.integer(x$ks)
obs <- as.numeric(x$obs_prop)
pred <- as.numeric(x$pred_prop)
lo <- as.numeric(x$lo_prop)
hi <- as.numeric(x$hi_prop)
if (length(ks) < 2L || !identical(ks, seq.int(0L, length(ks) - 1L))) {
  stop("Expected ks = 0:(K) contiguous; got: ", paste(ks, collapse = ","))
}
# Keep exact 0:(K-1); replace last bin with P(y >= K) = 1 - sum_{0:(K-1)}.
# Y1rep was not saved, so pointwise 95% bands are shown only for 0:(K-1);
# the pooled tail shows the predictive mean only.
k_cut <- ks[length(ks)]
idx_head <- seq_len(length(ks) - 1L)
obs_plot <- c(obs[idx_head], 1 - sum(obs[idx_head]))
pred_plot <- c(pred[idx_head], 1 - sum(pred[idx_head]))
lo_head <- lo[idx_head]
hi_head <- hi[idx_head]
lab <- c(as.character(ks[idx_head]), paste0(k_cut, "+"))

# Okabe–Ito blue / amber (vs Paper 1 teal #1B9E91 + red #D7301F)
col_obs  <- "#0072B2"
col_pred <- "#E69F00"
ylim_top <- max(c(obs_plot, hi_head), na.rm = TRUE) * 1.08

out_pdf <- file.path(OUT_DIR, "hrs_ppc_ab.pdf")
grDevices::pdf(out_pdf, width = 8.2, height = 9.2)
op <- graphics::par(mfrow = c(2, 1), mar = c(4.2, 4.2, 2.2, 1.2))
bp <- graphics::barplot(
  obs_plot,
  names.arg = lab,
  col = col_obs,
  border = NA,
  ylim = c(0, ylim_top),
  xlab = "HSPNIT count",
  ylab = "Proportion",
  main = "(a)"
)
graphics::arrows(
  x0 = bp[idx_head], y0 = lo_head, y1 = hi_head,
  angle = 90, code = 3, length = 0.04, col = col_pred
)
graphics::points(bp, pred_plot, pch = 16, col = col_pred, cex = 0.9)
graphics::hist(
  x$rep_cor_y1_y2,
  breaks = 30,
  col = col_pred,
  border = "white",
  main = "(b)",
  xlab = expression(cor(HSPNIT, SHLT)),
  ylab = "Count"
)
graphics::abline(v = x$obs_cor_y1_y2, col = col_obs, lwd = 2)
graphics::par(op)
grDevices::dev.off()

cat("Wrote", out_pdf, "\n")
cat(sprintf(
  "Tail %d+: obs=%.5f pred=%.5f (bands omitted; Y1rep not in RDS)\n",
  k_cut, obs_plot[length(obs_plot)], pred_plot[length(pred_plot)]
))
