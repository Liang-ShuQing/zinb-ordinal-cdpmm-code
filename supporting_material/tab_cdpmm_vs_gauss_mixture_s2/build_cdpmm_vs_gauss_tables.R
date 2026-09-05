############################################################
# Build CDPMM vs Gaussian LaTeX tables from prog2 HPC CSVs.
# Includes fixed effects, r, Sigma, and rho.
# Bolds Bias/RMSE/CP when estimation is poor:
#   |Bias| >= 0.10  or  CP < 0.90
# Reads folders like s1_n200_normal/ under PROG2_ROOT.
# Skips incomplete cells (no joint_param_summary_*.csv).
############################################################
args <- commandArgs(trailingOnly = TRUE)

here <- normalizePath(".", winslash = "/", mustWork = TRUE)
parent <- dirname(here)
proj_root <- if (basename(here) %in% c("模拟研究工具", "sim_tools")) parent else here

PROG2_ROOT <- if (length(args) >= 1L && nzchar(args[[1]])) {
  args[[1]]
} else {
  file.path(proj_root, "模拟研究结果", "prog2_n200n400_c10k_20260817")
}
OUT_TEX <- if (length(args) >= 2L && nzchar(args[[2]])) {
  args[[2]]
} else {
  file.path(proj_root, "paper", "sim_cdpmm_vs_gauss_tables.tex")
}

BIAS_POOR <- 0.10
CP_POOR <- 0.90

# Preferred display order; missing params in a cell are skipped.
PARAM_ORDER <- c(
  paste0("alpha", 1:8),
  paste0("beta", 1:8),
  paste0("gamma", 1:8),
  "r",
  "Sigma11", "Sigma22", "Sigma33",
  "rho1", "rho2", "rho3"
)

tex_param <- function(nm) {
  if (nm == "r") return("$r$")
  if (grepl("^alpha", nm)) return(sprintf("$\\alpha_{%s}$", sub("alpha", "", nm)))
  if (grepl("^beta", nm)) return(sprintf("$\\beta_{%s}$", sub("beta", "", nm)))
  if (grepl("^gamma", nm)) return(sprintf("$\\gamma_{%s}$", sub("gamma", "", nm)))
  if (nm == "Sigma11") return("$\\sigma_1^2$")
  if (nm == "Sigma22") return("$\\sigma_2^2$")
  if (nm == "Sigma33") return("$\\sigma_3^2$")
  if (nm == "rho1") return("$\\rho_{12}$")
  if (nm == "rho2") return("$\\rho_{13}$")
  if (nm == "rho3") return("$\\rho_{23}$")
  nm
}

fmt_num <- function(x, d = 3L) {
  if (is.na(x) || !is.finite(x)) return("---")
  sprintf(paste0("%.", d, "f"), x)
}
fmt_signed_raw <- function(x, d = 3L) {
  if (is.na(x) || !is.finite(x)) return("---")
  s <- fmt_num(abs(x), d)
  if (x < 0) paste0("$-", s, "$") else s
}
fmt_true <- function(x) {
  if (is.na(x) || !is.finite(x)) return("---")
  if (abs(x - round(x)) < 1e-8) return(sprintf("%.1f", x))
  if (abs(10 * x - round(10 * x)) < 1e-8) return(sprintf("%.1f", x))
  fmt_num(x, 3L)
}
maybe_bold <- function(s, bold) {
  if (!isTRUE(bold) || identical(s, "---")) return(s)
  paste0("\\textbf{", s, "}")
}

is_poor <- function(bias, cp) {
  # retained for diagnostics; bolding is applied separately to Bias and CP only
  bad_bias <- is.finite(bias) && abs(bias) >= BIAS_POOR
  bad_cp <- is.finite(cp) && cp < CP_POOR
  isTRUE(bad_bias || bad_cp)
}

parse_cell <- function(nm) {
  m <- regexec("^s([12])_n([0-9]+)_(normal|mixture)$", nm)
  g <- regmatches(nm, m)[[1]]
  if (length(g) < 4L) return(NULL)
  list(scenario = as.integer(g[2]), n = as.integer(g[3]), re_dist = g[4])
}

read_pair <- function(dir) {
  fc <- list.files(dir, pattern = "^joint_param_summary_cdpmm_.*\\.csv$", full.names = TRUE)
  fg <- list.files(dir, pattern = "^joint_param_summary_gauss_.*\\.csv$", full.names = TRUE)
  if (!length(fc) || !length(fg)) return(NULL)
  list(cdpmm = utils::read.csv(fc[[1]], stringsAsFactors = FALSE),
       gauss = utils::read.csv(fg[[1]], stringsAsFactors = FALSE))
}

dirs <- list.dirs(PROG2_ROOT, recursive = FALSE, full.names = TRUE)
cells <- list()
for (d in dirs) {
  meta <- parse_cell(basename(d))
  if (is.null(meta)) next
  pair <- read_pair(d)
  if (is.null(pair)) {
    message("Skip incomplete: ", basename(d))
    next
  }
  meta$dir <- d
  meta$pair <- pair
  cells[[length(cells) + 1L]] <- meta
}
if (!length(cells)) stop("No complete prog2 cells in ", PROG2_ROOT)

block_tex <- function(cdpmm, gauss) {
  present <- intersect(PARAM_ORDER, intersect(cdpmm$Parameter, gauss$Parameter))
  lines <- character(0)
  for (p in present) {
    rc <- cdpmm[cdpmm$Parameter == p, , drop = FALSE]
    rg <- gauss[gauss$Parameter == p, , drop = FALSE]
    if (!nrow(rc) || !nrow(rg)) next
    tv <- rc$True_Value[1]
    bc <- rc$Bias[1]; rmsec <- rc$RMSE[1]; cpc <- rc$CP[1]
    bg <- rg$Bias[1]; rmseg <- rg$RMSE[1]; cpg <- rg$CP[1]
    bold_bc <- is.finite(bc) && abs(bc) >= BIAS_POOR
    bold_cpc <- is.finite(cpc) && cpc < CP_POOR
    bold_bg <- is.finite(bg) && abs(bg) >= BIAS_POOR
    bold_cpg <- is.finite(cpg) && cpg < CP_POOR
    lines <- c(lines, sprintf(
      "%s & %s & %s & %s & %s & %s & %s & %s \\\\",
      tex_param(p),
      fmt_true(tv),
      maybe_bold(fmt_signed_raw(bc), bold_bc),
      fmt_num(rmsec),
      maybe_bold(fmt_num(cpc), bold_cpc),
      maybe_bold(fmt_signed_raw(bg), bold_bg),
      fmt_num(rmseg),
      maybe_bold(fmt_num(cpg), bold_cpg)
    ))
  }
  paste(lines, collapse = "\n")
}

build_table <- function(re_dist, scenario, label, caption) {
  want <- Filter(function(z) {
    identical(z$re_dist, re_dist) && identical(z$scenario, as.integer(scenario))
  }, cells)
  if (!length(want)) {
    message("No cells for ", re_dist, " Scenario ", scenario)
    return("")
  }
  ord <- order(vapply(want, `[[`, integer(1), "n"))
  want <- want[ord]
  parts <- character(0)
  for (i in seq_along(want)) {
    z <- want[[i]]
    hdr <- sprintf(
      "\\multicolumn{8}{@{}l}{\\textbf{$N=%d$}} \\\\",
      z$n
    )
    if (i > 1L) parts <- c(parts, "\\midrule", hdr, "\\midrule")
    else parts <- c(parts, hdr, "\\midrule")
    parts <- c(parts, block_tex(z$pair$cdpmm, z$pair$gauss))
  }
  paste0(
    "\\begin{table}[htbp]\n",
    "\\centering\n",
    "\\caption{", caption, "}\n",
    "\\label{", label, "}\n",
    "\\begin{adjustbox}{max width=\\textwidth,center}\n",
    "\\scriptsize\n",
    "\\setlength{\\tabcolsep}{3.5pt}\n",
    "\\begin{tabular}{lccccccc}\n",
    "\\toprule\n",
    " & & \\multicolumn{3}{c}{\\textbf{CDPMM Joint}} & \\multicolumn{3}{c}{\\textbf{Gaussian Joint}} \\\\\n",
    "\\cmidrule(lr){3-5}\\cmidrule(lr){6-8}\n",
    "\\textbf{Parameter} & \\textbf{True Value} & Bias & RMSE & CP & Bias & RMSE & CP \\\\\n",
    "\\midrule\n",
    paste(parts, collapse = "\n"), "\n",
    "\\bottomrule\n",
    "\\end{tabular}\n",
    "\\end{adjustbox}\n",
    "\\end{table}\n"
  )
}

note <- paste0(
  " Bold Bias entries mark $|\\mathrm{Bias}|\\ge ", sprintf("%.2f", BIAS_POOR),
  "$; bold CP entries mark $\\mathrm{CP}<", sprintf("%.2f", CP_POOR), "$ (RMSE is not bolded)."
)
cap_one <- function(re_lab, scen) {
  paste0(
    "Bias, RMSE, and coverage for fixed effects, dispersion, and random-effects ",
    "variances/correlations under a CDPMM versus Gaussian joint prior ",
    "(", re_lab, "; Scenario~", scen, ", $N\\in\\{200,400\\}$).",
    note
  )
}

tex <- paste(
  "% Auto-generated by build_cdpmm_vs_gauss_tables.R (CDPMM vs Gaussian joint).",
  "% Incomplete cells (no CSV) are skipped.",
  "% Four tables: normal/mixture x Scenario 1/2; FE + r + Sigma + rho.",
  "% Bolds Bias if |Bias|>=0.10 and CP if CP<0.90; RMSE never bolded.",
  "% Requires: \\usepackage{booktabs,adjustbox}",
  "",
  build_table("normal", 1L, "tab:cdpmm_vs_gauss_normal_s1",
              cap_one("multivariate normal random effects", 1L)),
  build_table("normal", 2L, "tab:cdpmm_vs_gauss_normal_s2",
              cap_one("multivariate normal random effects", 2L)),
  build_table("mixture", 1L, "tab:cdpmm_vs_gauss_mixture_s1",
              cap_one("multivariate mixture random effects", 1L)),
  build_table("mixture", 2L, "tab:cdpmm_vs_gauss_mixture_s2",
              cap_one("multivariate mixture random effects", 2L)),
  sep = "\n"
)

dir.create(dirname(OUT_TEX), recursive = TRUE, showWarnings = FALSE)
con <- file(OUT_TEX, open = "wt", encoding = "UTF-8")
on.exit(close(con), add = TRUE)
writeLines(tex, con, useBytes = FALSE)
message("Wrote ", OUT_TEX)
message("Complete cells: ", paste(vapply(cells, function(z) {
  sprintf("s%d_n%d_%s", z$scenario, z$n, z$re_dist)
}, character(1)), collapse = ", "))
message("Bold rule: Bias if |Bias| >= ", BIAS_POOR, "; CP if CP < ", CP_POOR, "; RMSE never bolded")
