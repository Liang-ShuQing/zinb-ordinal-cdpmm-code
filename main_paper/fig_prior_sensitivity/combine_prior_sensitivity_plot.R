############################################################
# Combine TYPE I/II/III prior-sensitivity CSVs into one figure
# Usage: Rscript combine_prior_sensitivity_plot.R <parent_or_csv_dir> [OUT_DIR]
############################################################

args <- commandArgs(trailingOnly = TRUE)
in_dir <- if (length(args) >= 1L) args[[1]] else "."
OUT_DIR <- if (length(args) >= 2L) args[[2]] else in_dir
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

library(ggplot2)

csvs <- list.files(in_dir, pattern = "^prior_sensitivity_.*\\.csv$",
                   full.names = TRUE, recursive = TRUE)
csvs <- csvs[!grepl("combined|_meta", basename(csvs))]
if (!length(csvs)) stop("No prior_sensitivity_*.csv under ", in_dir)

tab <- do.call(rbind, lapply(csvs, utils::read.csv, stringsAsFactors = FALSE))
tab <- tab[is.finite(tab$kappa) & is.finite(tab$ARMSE) & is.finite(tab$ACP), ]
tab$Type <- factor(tab$Type, levels = c("TYPE I", "TYPE II", "TYPE III"))
stopifnot(nrow(tab) > 0)

df <- rbind(
  data.frame(hyper = tab$kappa, value = tab$ARMSE, Type = tab$Type, Metric = "ARMSE"),
  data.frame(hyper = tab$kappa, value = tab$ACP, Type = tab$Type, Metric = "ACP")
)
df$Metric <- factor(df$Metric, levels = c("ARMSE", "ACP"))
hyper <- sort(unique(tab$kappa))
ref95 <- data.frame(Metric = factor("ACP", levels = c("ARMSE", "ACP")), yint = 0.95)

n_tag <- if (any(grepl("_n200", csvs))) {
  "n200"
} else if (any(grepl("_n100", csvs))) {
  "n100"
} else if (median(tab$ARMSE, na.rm = TRUE) < 0.14) {
  "n200"
} else {
  "n100"
}

# N=200 ARMSE sits near 0.11; N=100 panel used ~0.12--0.30
ylim_armse <- if (identical(n_tag, "n200")) c(0.08, 0.22) else c(0.12, 0.30)
pad <- rbind(
  data.frame(
    Metric = factor("ARMSE", levels = c("ARMSE", "ACP")),
    hyper = min(hyper), value = ylim_armse,
    Type = factor("TYPE I", levels = c("TYPE I", "TYPE II", "TYPE III"))
  ),
  data.frame(
    Metric = factor("ACP", levels = c("ARMSE", "ACP")),
    hyper = min(hyper), value = c(min(tab$ACP, na.rm = TRUE), 1.0),
    Type = factor("TYPE I", levels = c("TYPE I", "TYPE II", "TYPE III"))
  )
)

pal <- c("TYPE I" = "#1B9E91", "TYPE II" = "#7B3FA0", "TYPE III" = "#E69500")
shapes <- c("TYPE I" = 16, "TYPE II" = 17, "TYPE III" = 15)
ltypes <- c("TYPE I" = "solid", "TYPE II" = "dashed", "TYPE III" = "dotdash")

p <- ggplot(df, aes(hyper, value, color = Type, shape = Type, linetype = Type)) +
  geom_hline(data = ref95, aes(yintercept = yint), color = "red", linewidth = 0.5) +
  geom_blank(data = pad) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2.2, fill = "white", stroke = 0.5) +
  facet_wrap(~Metric, ncol = 1, scales = "free_y") +
  scale_x_log10(breaks = hyper, labels = hyper) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 5),
    labels = scales::label_number(accuracy = 0.01),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_color_manual(values = pal) +
  scale_shape_manual(values = shapes) +
  scale_linetype_manual(values = ltypes) +
  labs(x = "Prior hyperparameter", y = NULL, color = NULL, shape = NULL, linetype = NULL) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.82, 0.88),
    legend.background = element_rect(fill = scales::alpha("white", 0.7), color = NA),
    legend.key.width = grid::unit(1.2, "cm"),
    strip.background = element_rect(fill = "grey85", color = "grey40"),
    strip.text = element_text(face = "bold", hjust = 0.02),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3)
  ) +
  geom_text(
    data = data.frame(
      Metric = factor("ACP", levels = c("ARMSE", "ACP")),
      hyper = min(hyper) * 1.2, value = 0.95
    ),
    aes(x = hyper, y = value, label = "0.95"),
    inherit.aes = FALSE, color = "red", vjust = -0.6, hjust = 0, size = 3.2
  )

pdf_file <- file.path(OUT_DIR, sprintf("prior_sensitivity_scen2_mixture_%s_combined.pdf", n_tag))
ggsave(pdf_file, p, width = 6.5, height = 6.5, device = grDevices::pdf)
utils::write.csv(
  tab,
  file.path(OUT_DIR, sprintf("prior_sensitivity_scen2_mixture_%s_combined.csv", n_tag)),
  row.names = FALSE
)
cat("Wrote ", pdf_file, "\n", sep = "")
