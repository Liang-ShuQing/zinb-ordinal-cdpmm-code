############################################################
# Combine normal + mixture G-sensitivity CSVs into Tang-style figure
# Usage: Rscript combine_G_sensitivity_plot.R <dir_with_csvs> [OUT_DIR]
############################################################

args <- commandArgs(trailingOnly = TRUE)
in_dir <- if (length(args) >= 1L) args[[1]] else "."
OUT_DIR <- if (length(args) >= 2L) args[[2]] else in_dir
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

library(ggplot2)

csvs <- list.files(in_dir, pattern = "^G_sensitivity_scen2_.*\\.csv$", full.names = TRUE)
# also accept recursive one-level if empty
if (!length(csvs)) {
  csvs <- list.files(in_dir, pattern = "^G_sensitivity_scen2_.*\\.csv$",
                     full.names = TRUE, recursive = TRUE)
}
csvs <- csvs[!grepl("_meta|_combined", csvs)]
if (!length(csvs)) stop("No G_sensitivity_scen2_*.csv found in ", in_dir)

tab <- do.call(rbind, lapply(csvs, utils::read.csv, stringsAsFactors = FALSE))
tab <- tab[is.finite(tab$G) & is.finite(tab$ARMSE) & is.finite(tab$ACP), ]
tab$Type <- factor(tab$re_dist, levels = c("normal", "mixture"))
stopifnot(nrow(tab) > 0)

n_tag <- if (any(grepl("_n200", basename(csvs)))) {
  "n200"
} else if (any(grepl("_n100", basename(csvs)))) {
  "n100"
} else {
  "combined"
}

df <- rbind(
  data.frame(G = tab$G, value = tab$ARMSE, Type = tab$Type, Metric = "ARMSE"),
  data.frame(G = tab$G, value = tab$ACP, Type = tab$Type, Metric = "ACP")
)
df$Metric <- factor(df$Metric, levels = c("ARMSE", "ACP"))
G_breaks <- sort(unique(tab$G))

# Tang-style generous y-ranges; N=200 ARMSE is lower than the N=100 panel
ylim_armse <- if (identical(n_tag, "n200")) c(0.08, 0.20) else c(0.12, 0.30)
ylim_acp <- c(0.85, 1.00)
pad <- rbind(
  data.frame(
    Metric = factor("ARMSE", levels = c("ARMSE", "ACP")),
    G = min(G_breaks), value = ylim_armse,
    Type = factor("normal", levels = c("normal", "mixture"))
  ),
  data.frame(
    Metric = factor("ACP", levels = c("ARMSE", "ACP")),
    G = min(G_breaks), value = ylim_acp,
    Type = factor("normal", levels = c("normal", "mixture"))
  )
)

ref95 <- data.frame(Metric = factor("ACP", levels = c("ARMSE", "ACP")), yint = 0.95)
pal <- c(normal = "#1B9E91", mixture = "#7B3FA0")
shapes <- c(normal = 16, mixture = 17)
ltypes <- c(normal = "solid", mixture = "dashed")

p <- ggplot(df, aes(G, value, color = Type, shape = Type, linetype = Type)) +
  geom_hline(data = ref95, aes(yintercept = yint), color = "red", linewidth = 0.5) +
  geom_blank(data = pad) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.4, fill = "white", stroke = 0.55) +
  facet_wrap(~Metric, ncol = 1, scales = "free_y") +
  scale_x_continuous(breaks = G_breaks, labels = G_breaks) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 5),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_color_manual(values = pal) +
  scale_shape_manual(values = shapes) +
  scale_linetype_manual(values = ltypes) +
  labs(x = expression(Truncation~level~G), y = NULL,
       color = NULL, shape = NULL, linetype = NULL) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.82, 0.88),
    legend.background = element_rect(fill = scales::alpha("white", 0.7), color = NA),
    legend.key.width = grid::unit(1.2, "cm"),
    strip.background = element_rect(fill = "grey85", color = "grey40"),
    strip.text = element_text(face = "bold", hjust = 0.02),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    axis.text.x = element_text(size = 9)
  ) +
  geom_text(
    data = data.frame(
      Metric = factor("ACP", levels = c("ARMSE", "ACP")),
      G = min(G_breaks) + 0.5, value = 0.95
    ),
    aes(x = G, y = value, label = "0.95"),
    inherit.aes = FALSE, color = "red", vjust = -0.6, hjust = 0, size = 3.2
  )

pdf_file <- file.path(OUT_DIR, sprintf("G_sensitivity_armse_acp_scen2_%s_combined.pdf", n_tag))
ggsave(pdf_file, p, width = 7.2, height = 7.2, device = grDevices::pdf)
utils::write.csv(tab, file.path(OUT_DIR, sprintf("G_sensitivity_scen2_%s_combined.csv", n_tag)), row.names = FALSE)
cat("Wrote ", pdf_file, "\n", sep = "")
