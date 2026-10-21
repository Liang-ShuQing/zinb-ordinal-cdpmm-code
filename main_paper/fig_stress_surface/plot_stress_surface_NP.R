#!/usr/bin/env Rscript
# Figure 9 (paper label fig:stress_surface): N x P stress surface.
# Three-panel heatmap: (a) ARMSE, (b) ACP, (c) max |bias| of zero-effect coefs.
# x = P (extra null covariates per submodel), y = N (sample size).
# Row labels are plain N values; ACP is rounded UP to the nearest 0.002.
#
# Usage:
#   Rscript plot_stress_surface_NP.R <input_dir> [output_dir]
#
# <input_dir> must contain stress_surface_long.csv with columns
#   N, HD, ARMSE, ACP, nullMaxBias
# where HD is the number of extra null covariates (P) and the 16 rows cover
# N in {50,100,200,400} x HD in {0,5,20,30}. See aggregate_stress_cells.R.
#
# ARMSE/ACP are the mean over the 19 base parameters
# (alpha1-4, beta1-5, gamma1-3, r, Sigma11, Sigma22, Sigma33, rho1-3);
# delta and the off-diagonal Sigma entries are excluded.
suppressWarnings(suppressMessages({ library(ggplot2); library(gridExtra); library(grid) }))

args <- commandArgs(trailingOnly = TRUE)
in_dir  <- if (length(args) >= 1L) args[[1]] else "."
out_dir <- if (length(args) >= 2L) args[[2]] else in_dir

f_in <- file.path(in_dir, "stress_surface_long.csv")
if (!file.exists(f_in)) stop("missing: ", f_in)
d <- read.csv(f_in, stringsAsFactors = FALSE)

Ns  <- c(50, 100, 200, 400)
HDs <- c(0, 5, 20, 30)
mk <- function(col) {
  m <- matrix(NA_real_, length(Ns), length(HDs), dimnames = list(as.character(Ns), paste0("+", HDs)))
  for (i in seq_len(nrow(d))) if (d$N[i] %in% Ns) m[as.character(d$N[i]), paste0("+", d$HD[i])] <- d[[col]][i]
  m
}
M_armse <- mk("ARMSE"); M_acp <- mk("ACP"); M_null <- mk("nullMaxBias")
M_acp <- ceiling(M_acp / 0.002 - 1e-9) * 0.002
write.csv(M_armse, file.path(out_dir, "stress_surface_NP_ARMSE.csv"))
write.csv(M_acp,   file.path(out_dir, "stress_surface_NP_ACP.csv"))
write.csv(M_null,  file.path(out_dir, "stress_surface_NP_nullMaxBias.csv"))

mk_heat <- function(mat, title, low, high) {
  df <- as.data.frame(as.table(mat)); names(df) <- c("N", "HD", "value")
  df$N  <- factor(df$N,  levels = rownames(mat))
  df$HD <- factor(df$HD, levels = colnames(mat))
  ggplot(df, aes(HD, N, fill = value)) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(aes(label = ifelse(is.na(value), "—", sprintf("%.3f", value))), size = 3.6) +
    scale_fill_gradient(low = low, high = high, na.value = "grey88") +
    labs(title = title, x = "P", y = "N") +
    theme_bw(base_size = 12) + theme(panel.grid = element_blank())
}
p1 <- mk_heat(M_armse, "(a) ARMSE", "white", "#c0392b")
p2 <- mk_heat(M_acp,   "(b) ACP",   "#2166ac", "white")
p3 <- mk_heat(M_null,  "(c) max |bias| of zero-effect coefs", "white", "#762a83")

pdf(file.path(out_dir, "stress_surface_NP_heatmap.pdf"), width = 14.0, height = 4.4)
grid.arrange(p1, p2, p3, ncol = 3)
invisible(dev.off())
cat("Wrote:", file.path(out_dir, "stress_surface_NP_heatmap.pdf"), "\n")
