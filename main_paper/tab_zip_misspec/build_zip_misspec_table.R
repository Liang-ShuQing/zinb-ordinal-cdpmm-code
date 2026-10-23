#!/usr/bin/env Rscript
# Build the main-text LaTeX table for the ZIP data-generating-process robustness check
# (Scenario 2 + Truth 2, N=200, S=500) in the SAME layout as Tables 1-9:
#   Parameter | True Value | Joint: Bias RMSE CP Avg.L | Separate: Bias RMSE CP Avg.L | Variance Ratio (Gain)
#
# Sources (already pulled from the cluster):
#   - slurm-46727005.out                      (job zm500_zip; joint + ordinal-only + ZINB-only blocks, with AL)
#   - joint_param_summary_*.csv               (joint Bias/RMSE/CP)
#   - variance_ratio_*.csv                    (Variance_Ratio_Mean)
#   - criteria_winrate_*.csv                  (WAIC/LOOIC win rates for the block header)
#
# Rows: alpha1-4, beta2-5, gamma1-3, rho12/13/23  (14 identified parameters).
# beta1 and r are deliberately OMITTED: they are not separately identified when the
# generated counts carry no overdispersion (see the text).
#
# Bolding follows the project rule: bold Bias if |Bias| >= 0.10, bold CP if CP < 0.90.
# Output: zip_misspec_table.tex in the sidecar results folder (LaTeX body only, no float wrapper).

args <- commandArgs(trailingOnly = TRUE)

script_dir <- local({
  f <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(f)) {
    dirname(normalizePath(sub("^--file=", "", f[[1]]), winslash = "/", mustWork = FALSE))
  } else {
    getwd()
  }
})
repo_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = FALSE)

run_dir <- if (length(args) >= 1L) args[[1]] else file.path(repo_root, "hpc_jobs", "zeromech_misspec_s500_n200_wuzhen128c1_20261016", "zip")
out_dir <- if (length(args) >= 2L) args[[2]] else file.path(script_dir, "results")
BIAS_POOR <- 0.10
CP_POOR   <- 0.90

# ---------------------------------------------------------------- parsing ----
# valid parameter names in the printed summary tables (guards against parsing other sections)
PARAM_RE <- "^(alpha|beta|gamma|delta)[0-9]+$|^r$|^rho[0-9]+$|^Sigma[0-9]+$|^sigma2_b[0-9]+$"

parse_rows <- function(lines) {
  out <- list()
  for (ln in lines) {
    if (!grepl("^\\s*[A-Za-z]", ln)) next
    f <- strsplit(trimws(ln), "[[:space:]]+")[[1]]
    if (length(f) < 7) next
    if (!grepl(PARAM_RE, f[1])) next
    num <- suppressWarnings(as.numeric(f[2:7]))
    if (anyNA(num)) next
    out[[length(out) + 1L]] <- data.frame(
      Parameter = f[1], True_Value = num[1], Bias = num[2], RMSE = num[3],
      CI_Lower = num[4], CI_Upper = num[5], CP = num[6],
      AL = num[5] - num[4], stringsAsFactors = FALSE
    )
  }
  if (!length(out)) return(NULL)
  do.call(rbind, out)
}

# Section-delimited blocks: every "===..." line starts a new section; the three
# parameter-summary sections (joint / ordinal-only / ZINB-only) are the ones that yield rows.
slurm_file <- list.files(run_dir, pattern = "^slurm-.*\\.out$", full.names = TRUE)[1]
lines <- readLines(slurm_file, warn = FALSE, encoding = "UTF-8")
sec <- grep("^=+", lines)
secs <- lapply(seq_along(sec), function(i) {
  start <- sec[i] + 1L
  end <- if (i < length(sec)) sec[i + 1L] - 1L else length(lines)
  parse_rows(lines[start:end])
})
secs <- secs[!vapply(secs, is.null, logical(1))]
stopifnot(length(secs) >= 3L)

# joint fixed effects (alpha/delta/beta/gamma/r) and joint correlations (rho) are printed
# in two separate sections; both belong to the joint block.
is_joint_fe <- vapply(secs, function(d) all(c("alpha1", "delta1", "beta1") %in% d$Parameter), logical(1))
is_joint_re <- vapply(secs, function(d) all(grepl("^rho", d$Parameter)), logical(1))
is_ord      <- vapply(secs, function(d) any(d$Parameter == "sigma2_b3") && !any(grepl("^(alpha|beta)", d$Parameter)), logical(1))
is_zinb     <- vapply(secs, function(d) all(c("alpha1", "beta1") %in% d$Parameter) && !any(grepl("^gamma", d$Parameter)), logical(1))
stopifnot(sum(is_joint_fe) == 1L, sum(is_joint_re) == 1L, sum(is_ord) == 1L, sum(is_zinb) == 1L)
blocks <- list(
  joint = rbind(secs[[which(is_joint_fe)]], secs[[which(is_joint_re)]]),
  ord   = secs[[which(is_ord)]],
  zinb  = secs[[which(is_zinb)]]
)

joint <- blocks$joint
vr <- read.csv(list.files(run_dir, pattern = "^variance_ratio_.*\\.csv$", full.names = TRUE)[1],
               stringsAsFactors = FALSE)
wr <- read.csv(list.files(run_dir, pattern = "^criteria_winrate_.*\\.csv$", full.names = TRUE)[1],
               stringsAsFactors = FALSE)
stopifnot(nrow(wr) == 1L)

# --------------------------------------------------------------- helpers ----
fmt_num <- function(x, d = 3L) {
  if (is.na(x) || !is.finite(x)) return("---")
  sprintf(paste0("%.", d, "f"), x)
}
fmt_signed <- function(x, d = 3L) {
  if (is.na(x) || !is.finite(x)) return("---")
  s <- fmt_num(abs(x), d)
  if (x < 0) paste0("$-", s, "$") else s
}
maybe_bold <- function(s, bold) if (isTRUE(bold) && s != "---") paste0("\\textbf{", s, "}") else s
bias_poor <- function(x) is.finite(x) && abs(x) >= BIAS_POOR
cp_poor   <- function(x) is.finite(x) && x < CP_POOR
fmt_true <- function(x, d) {
  if (is.na(x)) return("---")
  if (abs(x - round(x)) < 1e-8) sprintf("%g", round(x)) else fmt_signed(x, d)
}
fmt_vr <- function(r) {
  if (is.na(r) || !is.finite(r)) return("---")
  sprintf("%.3f (%.1f\\%%)", r, 100 * (r - 1))
}
tex_param <- function(nm) {
  if (grepl("^alpha", nm)) return(sprintf("$\\alpha_{%s}$", sub("alpha", "", nm)))
  if (grepl("^beta",  nm)) return(sprintf("$\\beta_{%s}$",  sub("beta",  "", nm)))
  if (grepl("^gamma", nm)) return(sprintf("$\\gamma_{%s}$", sub("gamma", "", nm)))
  if (nm == "rho12") return("$\\rho_{12}$")
  if (nm == "rho13") return("$\\rho_{13}$")
  if (nm == "rho23") return("$\\rho_{23}$")
  nm
}
pick <- function(df, nm) {
  if (is.null(df)) return(NULL)
  i <- match(nm, df$Parameter)
  if (is.na(i)) NULL else df[i, , drop = FALSE]
}

# --------------------------------------------------------------- table ------
rows <- c(paste0("alpha", 1:4), paste0("beta", 2:5), paste0("gamma", 1:3),
          "rho12", "rho13", "rho23")

body <- character(0)
for (nm in rows) {
  j <- pick(joint, nm)
  if (is.null(j)) next
  is_rho <- grepl("^rho", nm)
  # separate-model columns: alpha/beta from the ZINB-only fit, gamma from the ordinal-only fit
  s <- if (is_rho) NULL else if (grepl("^gamma", nm)) pick(blocks$ord, nm) else pick(blocks$zinb, nm)
  vr_i <- if (is_rho) NA_real_ else {
    k <- match(nm, vr$Parameter)
    if (is.na(k)) NA_real_ else vr$Variance_Ratio_Mean[k]
  }
  td <- if (is_rho) 3L else 1L
  cells <- c(
    tex_param(nm),
    fmt_true(j$True_Value, td),
    maybe_bold(fmt_signed(j$Bias), bias_poor(j$Bias)),
    fmt_num(j$RMSE),
    maybe_bold(fmt_num(j$CP), cp_poor(j$CP)),
    fmt_num(j$AL),
    if (is.null(s)) "---" else maybe_bold(fmt_signed(s$Bias), bias_poor(s$Bias)),
    if (is.null(s)) "---" else fmt_num(s$RMSE),
    if (is.null(s)) "---" else maybe_bold(fmt_num(s$CP), cp_poor(s$CP)),
    if (is.null(s)) "---" else fmt_num(s$AL),
    fmt_vr(vr_i)
  )
  body <- c(body, paste0(paste(cells, collapse = " & "), " \\\\"))
}

corpus <- data.frame(
  Parameter = rows, stringsAsFactors = FALSE)
corpus$Joint_Bias <- vapply(rows, function(nm) { j <- pick(joint, nm); if (is.null(j)) NA_real_ else j$Bias }, numeric(1))
corpus$Joint_CP   <- vapply(rows, function(nm) { j <- pick(joint, nm); if (is.null(j)) NA_real_ else j$CP   }, numeric(1))
corpus$Sep_Bias   <- vapply(rows, function(nm) {
  if (grepl("^rho", nm)) return(NA_real_)
  s <- if (grepl("^gamma", nm)) pick(blocks$ord, nm) else pick(blocks$zinb, nm)
  if (is.null(s)) NA_real_ else s$Bias }, numeric(1))
corpus$Sep_CP     <- vapply(rows, function(nm) {
  if (grepl("^rho", nm)) return(NA_real_)
  s <- if (grepl("^gamma", nm)) pick(blocks$ord, nm) else pick(blocks$zinb, nm)
  if (is.null(s)) NA_real_ else s$CP }, numeric(1))

tex <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  paste0("\\caption{The parameter estimation results under Scenario~2 with Truth~2 and a ",
         "ZIP count-part data-generating process ($N=200$, $S=500$).}"),
  "\\label{tab:zip_misspec}",
  "\\begin{adjustbox}{max width=\\textwidth,center}",
  "\\scriptsize",
  "\\setlength{\\tabcolsep}{2.8pt}",
  "\\begin{tabular}{lcccccccccc}",
  "\\toprule",
  " & & \\multicolumn{4}{c}{\\textbf{Joint Model}} & \\multicolumn{4}{c}{\\textbf{Separate Models}} & \\textbf{Variance Ratio} \\\\",
  "\\cmidrule(lr){3-6}\\cmidrule(lr){7-10}",
  "\\textbf{Parameter} & \\textbf{True Value} & Bias & RMSE & CP & \\AvgLhdr & Bias & RMSE & CP & \\AvgLhdr & \\textbf{(Gain)} \\\\",
  "\\midrule",
  sprintf("\\multicolumn{11}{@{}l}{\\textbf{$N=200$} (WAIC/LOOIC win rates: $%.1f\\%%/%.1f\\%%$)} \\\\",
          wr$WAIC_WinPct, wr$LOOIC_WinPct),
  "\\midrule",
  body,
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{adjustbox}",
  "\\end{table}"
)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
writeLines(tex, file.path(out_dir, "zip_misspec_table.tex"), useBytes = TRUE)

cat("=== parsed blocks ===\n")
cat(sprintf("  sections with rows: %s\n", paste(vapply(secs, nrow, integer(1)), collapse = ", ")))
cat(sprintf("  joint rows=%d  ord rows=%d  zinb rows=%d\n",
            nrow(blocks$joint), nrow(blocks$ord), nrow(blocks$zinb)))
cat(sprintf("  win rates: DIC %.1f%%  WAIC %.1f%%  LOOIC %.1f%%\n",
            wr$DIC_WinPct, wr$WAIC_WinPct, wr$LOOIC_WinPct))
cat("\n=== bolding triggers ===\n")
cat("  joint Bias >= 0.10 : ",
    paste(corpus$Parameter[is.finite(corpus$Joint_Bias) & abs(corpus$Joint_Bias) >= BIAS_POOR], collapse = ", "), "\n", sep = "")
cat("  joint CP < 0.90    : ",
    paste(corpus$Parameter[is.finite(corpus$Joint_CP) & corpus$Joint_CP < CP_POOR], collapse = ", "), "\n", sep = "")
cat("  separate Bias>=0.10: ",
    paste(corpus$Parameter[is.finite(corpus$Sep_Bias) & abs(corpus$Sep_Bias) >= BIAS_POOR], collapse = ", "), "\n", sep = "")
cat("  separate CP < 0.90 : ",
    paste(corpus$Parameter[is.finite(corpus$Sep_CP) & corpus$Sep_CP < CP_POOR], collapse = ", "), "\n", sep = "")
cat("\n=== LaTeX written to ===\n  ", file.path(out_dir, "zip_misspec_table.tex"), "\n", sep = "")
