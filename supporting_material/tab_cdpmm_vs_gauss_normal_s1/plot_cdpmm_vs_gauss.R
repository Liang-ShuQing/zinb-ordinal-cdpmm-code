############################################################
# Visualize CDPMM vs Gaussian RE recovery
# Reads joint_param_summary_*.csv; no refitting.
#
# Accepts either:
#   - prog2_n200n400_.../s1_n200_normal/   (preferred new layout)
#   - prog2_mirror/prog2_s1_n100_.../      (legacy)
#
# Layout:
#   (a) CP advantage heatmap: CP(CDPMM) - CP(Gauss)
#   (b) RMSE advantage heatmap: RMSE(Gauss) - RMSE(CDPMM)
# Facets: random-effects truth (normal / mixture); columns: Scenario x N
# Focus: Sigma11/22/33 + rho12/13/23; green = CDPMM better
#
# Usage:
#   Rscript plot_cdpmm_vs_gauss.R [PROG2_DIR] [OUT_DIR]
############################################################

args <- commandArgs(trailingOnly = TRUE)

# Prefer running from 模拟研究工具/
here <- normalizePath(".", winslash = "/", mustWork = TRUE)
parent <- dirname(here)
if (basename(here) == "模拟研究工具" || basename(here) == "sim_tools") {
  proj_root <- parent
} else if (dir.exists(file.path(here, "模拟研究结果"))) {
  proj_root <- here
} else {
  proj_root <- parent
}

mirror_default <- file.path(proj_root, "模拟研究结果", "prog2_n200n400_c10k_20260817")
if (!dir.exists(mirror_default)) {
  mirror_default <- file.path(proj_root, "模拟研究结果", "prog2_mirror")
}
out_default <- file.path(proj_root, "模拟研究结果")

mirror_dir <- if (length(args) >= 1L && nzchar(args[[1]])) args[[1]] else mirror_default
OUT_DIR <- if (length(args) >= 2L && nzchar(args[[2]])) args[[2]] else out_default
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

PARAM_KEEP <- c("Sigma11", "Sigma22", "Sigma33", "rho1", "rho2", "rho3")

dirs <- list.dirs(mirror_dir, recursive = FALSE, full.names = TRUE)
dirs_new <- dirs[grepl("^s[12]_n[0-9]+_(normal|mixture)$", basename(dirs))]
dirs_old <- dirs[grepl("prog2_s[12]_n(100|200|400)_(normal|mixture)_", basename(dirs))]
dirs <- if (length(dirs_new)) dirs_new else dirs_old
if (!length(dirs)) stop("No prog2 cell dirs found in ", mirror_dir)

parse_tag <- function(nm) {
  # New: s1_n200_mixture ; legacy: prog2_s1_n100_mixture_nsim200_20260803
  m <- regexec("(?:^|prog2_)s([12])_n([0-9]+)_(normal|mixture)", nm)
  g <- regmatches(nm, m)[[1]]
  if (length(g) < 4L) return(NULL)
  list(scenario = as.integer(g[2]), n = as.integer(g[3]), re_dist = g[4])
}

read_one <- function(path, method) {
  d <- utils::read.csv(path, stringsAsFactors = FALSE)
  d <- d[d$Parameter %in% PARAM_KEEP, , drop = FALSE]
  d$Method <- method
  d
}

rows <- list()
idx <- 1L
for (d in dirs) {
  meta <- parse_tag(basename(d))
  if (is.null(meta)) next
  f_c <- list.files(d, pattern = "^joint_param_summary_cdpmm_.*\\.csv$", full.names = TRUE)
  f_g <- list.files(d, pattern = "^joint_param_summary_gauss_.*\\.csv$", full.names = TRUE)
  if (!length(f_c) || !length(f_g)) next
  tc <- read_one(f_c[[1]], "CDPMM")
  tg <- read_one(f_g[[1]], "Gaussian")
  both <- rbind(tc, tg)
  both$Scenario <- meta$scenario
  both$N <- meta$n
  both$Truth <- meta$re_dist
  rows[[idx]] <- both
  idx <- idx + 1L
}
tab <- bind_rows(rows)
if (!nrow(tab)) stop("Failed to assemble CDPMM/Gaussian summaries.")

design_levels <- unique(sprintf("S%d\nN=%d", tab$Scenario, tab$N))
design_ord <- order(
  as.integer(sub("^S([0-9]+).*", "\\1", design_levels)),
  as.integer(sub(".*N=([0-9]+)$", "\\1", gsub("\n", "", design_levels)))
)
design_levels <- design_levels[design_ord]

tab <- tab %>%
  mutate(
    Truth = factor(Truth, levels = c("normal", "mixture"),
                   labels = c("Normal truth", "Mixture truth")),
    Design = factor(sprintf("S%d\nN=%d", Scenario, N), levels = design_levels),
    Param = factor(Parameter, levels = PARAM_KEEP),
    Method = factor(Method, levels = c("CDPMM", "Gaussian"))
  )

param_labels <- c(
  Sigma11 = "sigma[1]^2",
  Sigma22 = "sigma[2]^2",
  Sigma33 = "sigma[3]^2",
  rho1 = "rho[12]",
  rho2 = "rho[13]",
  rho3 = "rho[23]"
)

col_pos <- "#1B9E77"
col_neg <- "#D95F02"

theme_heat <- theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    legend.position = "bottom",
    legend.title = element_text(size = 9),
    axis.text.x = element_text(size = 9, lineheight = 0.95),
    axis.text.y = element_text(size = 10),
    plot.tag = element_text(face = "bold", size = 12)
  )

# Wide differences: green => CDPMM better
tab_diff <- tab %>%
  select(Truth, Design, Param, Method, CP, RMSE) %>%
  pivot_wider(names_from = Method, values_from = c(CP, RMSE)) %>%
  mutate(
    CP_diff = CP_CDPMM - CP_Gaussian,
    RMSE_diff = RMSE_Gaussian - RMSE_CDPMM
  )

p_cp <- ggplot(tab_diff, aes(x = Design, y = Param, fill = CP_diff)) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(aes(label = sprintf("%+.2f", CP_diff)), size = 3.15, color = "grey10") +
  facet_wrap(~Truth, nrow = 1) +
  scale_y_discrete(labels = function(x) parse(text = param_labels[x]), limits = rev(PARAM_KEEP)) +
  scale_fill_gradient2(
    low = col_neg, mid = "grey96", high = col_pos, midpoint = 0,
    name = expression(CP[CDPMM] - CP[Gauss])
  ) +
  labs(x = NULL, y = NULL) +
  theme_heat +
  theme(legend.position = "right")

p_rmse <- ggplot(tab_diff, aes(x = Design, y = Param, fill = RMSE_diff)) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(aes(label = sprintf("%+.3f", RMSE_diff)), size = 3.0, color = "grey10") +
  facet_wrap(~Truth, nrow = 1) +
  scale_y_discrete(labels = function(x) parse(text = param_labels[x]), limits = rev(PARAM_KEEP)) +
  scale_fill_gradient2(
    low = col_neg, mid = "grey96", high = col_pos, midpoint = 0,
    name = "RMSE(Gauss) - RMSE(CDPMM)"
  ) +
  labs(x = NULL, y = NULL) +
  theme_heat +
  theme(legend.position = "right")

p <- (p_cp / p_rmse) +
  plot_annotation(
    tag_levels = "a",
    tag_prefix = "(",
    tag_suffix = ")"
  )

out_pdf <- file.path(OUT_DIR, "cdpmm_vs_gauss_cp_bias.pdf")
ggplot2::ggsave(
  filename = out_pdf,
  plot = p,
  width = 10.8,
  height = 7.6,
  device = grDevices::pdf
)

cat("Wrote:\n  ", out_pdf, "\n", sep = "")
cat(sprintf("CP advantage range: [%.3f, %.3f]\n",
            min(tab_diff$CP_diff, na.rm = TRUE),
            max(tab_diff$CP_diff, na.rm = TRUE)))
cat(sprintf("RMSE advantage range: [%.3f, %.3f]\n",
            min(tab_diff$RMSE_diff, na.rm = TRUE),
            max(tab_diff$RMSE_diff, na.rm = TRUE)))
cat(sprintf("Assembled cells: %d rows (long), %d tiles\n", nrow(tab), nrow(tab_diff)))