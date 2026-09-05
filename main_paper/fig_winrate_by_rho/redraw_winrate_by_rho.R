# Restyle win-rate-by-rho figure from saved CSV (no refit).
# Outputs PDF into a sibling results folder under 模拟研究结果/.

csv_default <- file.path(
  dirname(normalizePath(".")),
  "模拟研究结果", "fig_n200_taiyuan_20260817", "winrate_rho_s1_normal",
  "winrate_by_rho_scen1_normal_n200_nsim128.csv"
)
args <- commandArgs(trailingOnly = TRUE)
csv_file <- if (length(args) >= 1L) args[[1]] else csv_default
if (!file.exists(csv_file)) {
  # When cwd is 模拟研究工具/
  csv_file <- file.path(
    "..", "模拟研究结果", "fig_n200_taiyuan_20260817", "winrate_rho_s1_normal",
    "winrate_by_rho_scen1_normal_n200_nsim128.csv"
  )
}
csv_file <- normalizePath(csv_file, winslash = "/", mustWork = TRUE)

out_dir <- Sys.getenv("OUT_DIR", unset = "")
if (!nzchar(out_dir)) {
  out_dir <- file.path(
    dirname(dirname(csv_file)), # fig_n200_...
    "winrate_rho_s1_normal_restyle"
  )
  # Prefer sibling of 模拟研究工具 -> 模拟研究结果/...
  tools_parent <- dirname(normalizePath(".", winslash = "/", mustWork = TRUE))
  cand <- file.path(
    tools_parent, "模拟研究结果", "fig_n200_taiyuan_20260817",
    "winrate_rho_s1_normal_restyle"
  )
  out_dir <- cand
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

tab <- utils::read.csv(csv_file, stringsAsFactors = FALSE)
tab <- tab[is.finite(tab$rho) & is.finite(tab$WAIC_WinPct) & is.finite(tab$LOOIC_WinPct), ]
tab <- tab[order(tab$rho), ]

pdf_file <- file.path(
  out_dir, "win_rate_waic_looic_by_rho_scen1_normal_n200_nsim128.pdf"
)

# Visual fixes vs original:
# - slight headroom above 100% so the |rho| large plateau is not glued to the frame
# - LOOIC open circles, WAIC filled, so overlap is readable
# - tiny labels where win rate is in (99, 100) so a 127/128 dip is not mistaken for a quirk
grDevices::pdf(pdf_file, width = 8.2, height = 5.2)
op <- graphics::par(mar = c(4.2, 4.2, 1.2, 1.2), mgp = c(2.4, 0.7, 0))
on.exit({
  graphics::par(op)
  grDevices::dev.off()
}, add = FALSE)

xlim <- range(tab$rho)
ylim <- c(0, 104)
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
graphics::abline(h = 100, col = "gray70", lty = 3, lwd = 1)
graphics::abline(v = 0, col = "gray50", lty = 2, lwd = 1.2)

col_waic <- "#E69F00"
col_loo <- "#0072B2"
graphics::lines(tab$rho, tab$WAIC_WinPct, col = col_waic, lwd = 2)
graphics::points(tab$rho, tab$WAIC_WinPct, pch = 16, col = col_waic, cex = 1.1)
graphics::lines(tab$rho, tab$LOOIC_WinPct, col = col_loo, lwd = 2)
graphics::points(
  tab$rho, tab$LOOIC_WinPct,
  pch = 21, col = col_loo, bg = "white", cex = 1.15, lwd = 1.6
)

# Annotate near-100% dips (one missed replication out of 128 => 99.21875)
lab_idx <- which(tab$WAIC_WinPct > 99 & tab$WAIC_WinPct < 100 |
  tab$LOOIC_WinPct > 99 & tab$LOOIC_WinPct < 100)
if (length(lab_idx)) {
  for (i in lab_idx) {
    y_lab <- max(tab$WAIC_WinPct[i], tab$LOOIC_WinPct[i])
    graphics::text(
      tab$rho[i], y_lab + 2.2,
      labels = sprintf("%.1f%%", y_lab),
      cex = 0.72, col = "gray25"
    )
  }
}

graphics::legend(
  "top",
  legend = c("WAIC Win (%)", "LOOIC Win (%)"),
  col = c(col_waic, col_loo),
  lty = 1, lwd = 2,
  pch = c(16, 21), pt.bg = c(NA, "white"), pt.cex = 1.1,
  bty = "n", horiz = TRUE, cex = 0.95
)
graphics::box()

cat("Wrote ", pdf_file, "\n", sep = "")
