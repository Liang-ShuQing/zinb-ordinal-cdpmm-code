############################################################
# Fig.1a-style: distribution of HSPNIT by SHLT (HRS)
# Matches first-paper panel (a): Frequency(%) of HSPNIT
# faceted by the second outcome (here ordinal SHLT, 3 levels).
# Usage: Rscript plot_joint_assoc_hspnit_shlt_fig1a.R [DATA_CSV] [OUT_DIR]
############################################################

args <- commandArgs(trailingOnly = TRUE)

code_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile, winslash = "/")),
  error = function(e) normalizePath(".", winslash = "/", mustWork = TRUE)
)
if (!nzchar(Sys.getenv("R_SCRIPT_PATH", ""))) {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg)) {
    code_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/"))
  }
}

default_sidecar <- function(code_dir) {
  parent <- dirname(code_dir)
  base <- basename(code_dir)
  if (grepl("[\u4e00-\u9fff]", base)) {
    file.path(parent, paste0(base, "结果"))
  } else {
    file.path(parent, paste0(base, "_results"))
  }
}

data_path <- if (length(args) >= 1L && nzchar(args[[1]])) {
  args[[1]]
} else {
  cand <- file.path(code_dir, "Data_model_complete.csv")
  if (file.exists(cand)) cand else {
    stop("Data_model_complete.csv not found: pass its path as the first argument.")
  }
}
OUT_DIR <- if (length(args) >= 2L && nzchar(args[[2]])) {
  args[[2]]
} else {
  default_sidecar(code_dir)
}
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

paper_fig_dir <- file.path(OUT_DIR, "paper_figures")
dir.create(paper_fig_dir, recursive = TRUE, showWarnings = FALSE)

cat("DATA =", data_path, "\n")
cat("OUT_DIR =", OUT_DIR, "\n")

dat_raw <- utils::read.csv(data_path, stringsAsFactors = FALSE)
if ("wave" %in% names(dat_raw)) {
  dat_raw <- dat_raw[order(dat_raw$id, dat_raw$wave), ]
} else {
  dat_raw <- dat_raw[order(dat_raw$id), ]
}

dat_raw$SHLT_num <- as.numeric(dat_raw$SHLT)
dat_raw$SHLT_recoded <- NA_integer_
dat_raw$SHLT_recoded[dat_raw$SHLT_num %in% c(1, 2)] <- 1L
dat_raw$SHLT_recoded[dat_raw$SHLT_num == 3] <- 2L
dat_raw$SHLT_recoded[dat_raw$SHLT_num %in% c(4, 5)] <- 3L

binary_code_vars <- intersect(c("HIBP", "DIAB", "LUNG"), names(dat_raw))
if (length(binary_code_vars)) {
  keep_bin <- rep(TRUE, nrow(dat_raw))
  for (var in binary_code_vars) {
    keep_bin <- keep_bin & !is.na(dat_raw[[var]]) & dat_raw[[var]] <= 1
  }
  dat_raw <- dat_raw[keep_bin, ]
}

covariate_vars <- intersect(c("SMOKEV", "HIBP", "DIAB", "LUNG"), names(dat_raw))
cc_vars <- c("HSPNIT", "SHLT_recoded", covariate_vars)
complete <- stats::complete.cases(dat_raw[, cc_vars, drop = FALSE])
dat <- dat_raw[complete, ]

y1 <- as.integer(dat$HSPNIT)
y2 <- as.integer(dat$SHLT_recoded)
stopifnot(all(y2 %in% 1:3))

shlt_lab <- c("Category 1", "Category 2", "Category 3")
shlt_f <- factor(y2, levels = 1:3, labels = shlt_lab)

# Frequency(%) of each HSPNIT value within SHLT (first-paper 1a style)
df <- data.frame(HSPNIT = y1, SHLT = shlt_f)
tab <- as.data.frame(table(HSPNIT = df$HSPNIT, SHLT = df$SHLT), responseName = "n")
tab$HSPNIT <- as.integer(as.character(tab$HSPNIT))
n_by_shlt <- tapply(tab$n, tab$SHLT, sum)
tab$freq <- tab$n / as.numeric(n_by_shlt[as.character(tab$SHLT)])
tab <- tab[tab$n > 0, ]

cat("n_obs =", nrow(df), "  max HSPNIT =", max(y1), "\n")
print(tapply(y1 == 0, shlt_f, mean))

if (!requireNamespace("ggplot2", quietly = TRUE)) stop("需要 ggplot2")
if (!requireNamespace("scales", quietly = TRUE)) stop("需要 scales")
library(ggplot2)

# Low end: first-paper dense ticks (readable under asinh).
# Tail: every-40 overcrowds under asinh; keep early tail, then jump.
# Practical set: 0,1,3,5,7,9, 41,81,121, 201, max
h_max <- max(tab$HSPNIT)
x_breaks <- c(0L, 1L, 3L, 5L, 7L, 9L, 41L, 81L, 121L, 201L)
x_breaks <- sort(unique(c(x_breaks[x_breaks < h_max], h_max)))
cat("x_breaks =", paste(x_breaks, collapse = ", "), "\n")

# Bar edges in original HSPNIT units (width ~1 night); asinh then expands near 0
tab$xmin <- pmax(tab$HSPNIT - 0.45, 0)
tab$xmax <- tab$HSPNIT + 0.45

asinh_tr <- scales::trans_new(
  name = "asinh",
  transform = asinh,
  inverse = sinh,
  domain = c(0, Inf)
)

col_shlt <- c(
  "Category 1" = "#4C78A8",
  "Category 2" = "#F2C14E",
  "Category 3" = "#E45756"
)
# Soft tint for (a) histograms so columns stay aligned with (b) but look distinct
col_a_fill <- c(
  "Category 1" = "#7A9BB8",
  "Category 2" = "#C9B06A",
  "Category 3" = "#C47A78"
)

theme_joint <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      strip.background = element_rect(fill = "grey94", color = "grey75", linewidth = 0.4),
      strip.text = element_text(face = "bold", size = 9.5),
      plot.title = element_text(face = "bold", size = 11, hjust = 0),
      axis.title = element_text(size = 10),
      plot.margin = margin(4, 6, 2, 6)
    )
}

p_a <- ggplot(tab, aes(ymin = 0, ymax = freq, fill = SHLT)) +
  geom_rect(aes(xmin = xmin, xmax = xmax), color = NA, alpha = 0.92) +
  facet_wrap(~SHLT, nrow = 1) +
  scale_fill_manual(values = col_a_fill, guide = "none") +
  scale_y_continuous(
    labels = function(x) sprintf("%.0f%%", 100 * x),
    expand = expansion(mult = c(0, 0.06)),
    limits = c(0, NA)
  ) +
  scale_x_continuous(
    trans = asinh_tr,
    breaks = x_breaks,
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  labs(
    title = "(a)  P(HSPNIT | SHLT)",
    x = "Hospital nights (HSPNIT)",
    y = "Proportion (%)"
  ) +
  theme_joint() +
  theme(axis.text.x = element_text(size = 7))

out_pdf <- file.path(OUT_DIR, "joint_hspnit_shlt_fig1a.pdf")
ggsave(out_pdf, p_a, width = 10.5, height = 3.55, device = grDevices::pdf)
file.copy(out_pdf, file.path(paper_fig_dir, "joint_hspnit_shlt_fig1a.pdf"), overwrite = TRUE)

# Panel (b): same 3-column SHLT layout, reverse conditioning P(SHLT | HSPNIT)
hsp_grp <- cut(
  y1,
  breaks = c(-Inf, 0, 1, 2, 3, 4, Inf),
  labels = c("0", "1", "2", "3", "4", "5+"),
  right = TRUE
)
df_b <- data.frame(
  HSPNIT = factor(hsp_grp, levels = c("0", "1", "2", "3", "4", "5+")),
  SHLT = shlt_f
)
tab_b <- prop.table(table(df_b$HSPNIT, df_b$SHLT), margin = 1)
df_b_plot <- as.data.frame.table(tab_b, responseName = "prop")
names(df_b_plot) <- c("HSPNIT", "SHLT", "prop")
df_b_plot$SHLT <- factor(df_b_plot$SHLT, levels = shlt_lab)

df_b_sum <- aggregate(prop ~ HSPNIT, data = df_b_plot, FUN = sum)
stopifnot(all(abs(df_b_sum$prop - 1) < 1e-8))

p_b <- ggplot(df_b_plot, aes(x = HSPNIT, y = prop, fill = SHLT)) +
  geom_col(width = 0.68, color = NA, alpha = 0.88) +
  facet_wrap(~SHLT, nrow = 1) +
  scale_fill_manual(values = col_shlt, guide = "none") +
  scale_y_continuous(
    labels = function(x) sprintf("%.0f%%", 100 * x),
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    title = "(b)  P(SHLT | HSPNIT)",
    x = "Hospital nights (HSPNIT)",
    y = "Proportion (%)"
  ) +
  theme_joint() +
  theme(
    panel.grid.major.y = element_line(color = "grey88", linewidth = 0.35)
  )

if (!requireNamespace("patchwork", quietly = TRUE)) {
  stop("需要 patchwork")
}
library(patchwork)
p_ab <- p_a / p_b + plot_layout(heights = c(1.05, 1.12))
out_ab <- file.path(OUT_DIR, "joint_hspnit_shlt_fig1ab.pdf")
ggsave(out_ab, p_ab, width = 10.5, height = 6.35, device = grDevices::pdf)
file.copy(out_ab, file.path(paper_fig_dir, "joint_hspnit_shlt_fig1ab.pdf"), overwrite = TRUE)

cat("Saved:\n  ", out_pdf, "\n  ", out_ab, "\n", sep = "")
