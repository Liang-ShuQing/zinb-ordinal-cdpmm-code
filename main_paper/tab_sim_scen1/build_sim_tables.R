############################################################
# Build Scenario 1/2 LaTeX tables from HPC outputs
# Columns: Bias, RMSE, CP, Avg.L (avg 95% CI length); plus variance ratio
# Header uses \\AvgLhdr (defined in paper.tex / supporting_material.tex)
# Default layout: N=200 & N=400 under ROOT (s*_n{200,400}_* or scen*)
############################################################
args <- commandArgs(trailingOnly = TRUE)
ROOT <- if (length(args) >= 1L) args[[1]] else {
  file.path(dirname(normalizePath(".")), "模拟研究结果")
}
OUT_TEX <- if (length(args) >= 2L) args[[2]] else {
  file.path(dirname(normalizePath(".")), "paper", "sim_tables_generated.tex")
}

normalize_param_name <- function(nm) {
  if (nm == "Sigma[1,1]") return("Sigma11")
  if (nm == "Sigma[2,2]") return("Sigma22")
  if (nm == "Sigma[3,3]") return("Sigma33")
  if (nm == "rho12") return("rho1")
  if (nm == "rho13") return("rho2")
  if (nm == "rho23") return("rho3")
  nm
}

# Parse printed summary rows: name true bias rmse ci_lo ci_hi cp
parse_rows <- function(chunk) {
  out <- list()
  for (ln in chunk) {
    ln <- trimws(ln)
    if (!grepl("^(alpha|beta|gamma|r|sigma2|Sigma\\[|rho)", ln)) next
    parts <- strsplit(ln, "\\s+")[[1]]
    if (length(parts) < 7) next
    nm <- normalize_param_name(parts[1])
    vals <- suppressWarnings(as.numeric(parts[2:7]))
    if (anyNA(vals)) next
    out[[length(out) + 1L]] <- data.frame(
      Parameter = nm,
      True_Value = vals[1],
      Bias = vals[2],
      RMSE = vals[3],
      AL = vals[5] - vals[4],
      CP = vals[6],
      stringsAsFactors = FALSE
    )
  }
  if (!length(out)) return(NULL)
  tab <- do.call(rbind, out)
  tab[!duplicated(tab$Parameter, fromLast = TRUE), , drop = FALSE]
}

parse_slurm_blocks <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  # Section anchors (Chinese headers from prog1)
  i_joint <- grep("参数估计汇总（联合模型）", lines)[1]
  i_re <- grep("随机效应方差汇总（联合模型）", lines)[1]
  i_ord <- grep("参数估计汇总（有序单独模型）", lines)[1]
  i_zinb <- grep("参数估计汇总（ZINB单独模型）", lines)[1]
  i_mc <- grep("模型比较汇总", lines)[1]
  if (anyNA(c(i_joint, i_ord, i_zinb))) {
    return(list(joint = NULL, ord = NULL, zinb = NULL))
  }
  end_fe <- if (!is.na(i_re)) i_re - 1L else i_ord - 1L
  end_re <- i_ord - 1L
  end_ord <- i_zinb - 1L
  end_zinb <- if (!is.na(i_mc)) i_mc - 1L else length(lines)

  joint_fe <- parse_rows(lines[i_joint:end_fe])
  joint_re <- if (!is.na(i_re)) parse_rows(lines[i_re:end_re]) else NULL
  joint <- if (is.null(joint_re)) joint_fe else rbind(joint_fe, joint_re)

  ord <- parse_rows(lines[i_ord:end_ord])
  zinb <- parse_rows(lines[i_zinb:end_zinb])
  if (!is.null(ord)) {
    keep <- grepl("^gamma", ord$Parameter) | ord$Parameter == "sigma2_b3"
    ord <- ord[keep, , drop = FALSE]
  }
  if (!is.null(zinb)) {
    keep <- grepl("^(alpha|beta)", zinb$Parameter) |
      zinb$Parameter %in% c("r", "sigma2_b1", "sigma2_b2")
    zinb <- zinb[keep, , drop = FALSE]
  }
  list(joint = joint, ord = ord, zinb = zinb)
}

BIAS_POOR <- 0.10
CP_POOR <- 0.90

fmt_num <- function(x, d = 3) {
  if (is.na(x) || !is.finite(x)) return("---")
  sprintf(paste0("%.", d, "f"), x)
}
fmt_signed <- function(x, d = 3) {
  if (is.na(x) || !is.finite(x)) return("---")
  s <- fmt_num(abs(x), d)
  if (x < 0) paste0("$-", s, "$") else s
}
maybe_bold <- function(s, bold) {
  if (!isTRUE(bold) || identical(s, "---")) return(s)
  paste0("\\textbf{", s, "}")
}
bias_poor <- function(x) is.finite(x) && abs(x) >= BIAS_POOR
cp_poor <- function(x) is.finite(x) && x < CP_POOR
fmt_true <- function(x, d = 3) {
  if (d <= 1L) return(fmt_signed(x, 1L))
  if (abs(x - round(x)) < 1e-8) sprintf("%g", round(x)) else fmt_signed(x, d)
}
fmt_vr <- function(r) {
  if (is.na(r) || !is.finite(r)) return("---")
  sprintf("%.3f (%.1f\\%%)", r, 100 * (r - 1))
}
tex_param <- function(nm) {
  if (nm == "r") return("$r$")
  if (grepl("^alpha", nm)) return(sprintf("$\\alpha_{%s}$", sub("alpha", "", nm)))
  if (grepl("^beta", nm)) return(sprintf("$\\beta_{%s}$", sub("beta", "", nm)))
  if (grepl("^gamma", nm)) return(sprintf("$\\gamma_{%s}$", sub("gamma", "", nm)))
  if (nm %in% c("Sigma11", "sigma2_b1")) return("$\\sigma_1^2$")
  if (nm %in% c("Sigma22", "sigma2_b2")) return("$\\sigma_2^2$")
  if (nm %in% c("Sigma33", "sigma2_b3")) return("$\\sigma_3^2$")
  if (nm %in% c("rho1", "rho12")) return("$\\rho_{12}$")
  if (nm %in% c("rho2", "rho13")) return("$\\rho_{13}$")
  if (nm %in% c("rho3", "rho23")) return("$\\rho_{23}$")
  nm
}

build_block <- function(joint, vr, sep_zinb, sep_ord, n_label,
                        true_digits = 3, true_digits_sigma_rho = NULL,
                        joint_al = NULL) {
  fe_rows <- c(
    grep("^alpha", joint$Parameter, value = TRUE),
    grep("^beta", joint$Parameter, value = TRUE),
    grep("^gamma", joint$Parameter, value = TRUE)
  )
  sort_fe <- function(nms, prefix) {
    idx <- as.integer(sub(prefix, "", nms))
    nms[order(idx)]
  }
  fe_rows <- c(
    sort_fe(grep("^alpha", fe_rows, value = TRUE), "alpha"),
    sort_fe(grep("^beta", fe_rows, value = TRUE), "beta"),
    sort_fe(grep("^gamma", fe_rows, value = TRUE), "gamma")
  )
  rows <- c(
    fe_rows,
    "r",
    "Sigma11", "Sigma22", "Sigma33",
    "rho1", "rho2", "rho3"
  )
  if (is.null(true_digits_sigma_rho)) true_digits_sigma_rho <- true_digits
  true_d_for <- function(nm) {
    if (grepl("^(Sigma|rho)", nm)) true_digits_sigma_rho else true_digits
  }
  vr_lookup <- function(nm) {
    if (nm == "Sigma11") nm2 <- "sigma2_b1"
    else if (nm == "Sigma22") nm2 <- "sigma2_b2"
    else if (nm == "Sigma33") nm2 <- "sigma2_b3"
    else if (nm %in% c("rho1", "rho2", "rho3")) return(NA_real_)
    else nm2 <- nm
    i <- match(nm2, vr$Parameter)
    if (is.na(i)) NA_real_ else vr$Variance_Ratio_Mean[i]
  }
  al_joint_lookup <- function(nm) {
    if (is.null(joint_al)) return(NA_real_)
    i <- match(nm, joint_al$Parameter)
    if (is.na(i)) NA_real_ else joint_al$AL[i]
  }
  sep_lookup <- function(nm) {
    # returns Bias, RMSE, AL, CP
    if (grepl("^gamma", nm) || nm == "Sigma33") {
      tab <- sep_ord
      key <- if (nm == "Sigma33") "sigma2_b3" else nm
    } else if (grepl("^(alpha|beta)", nm) || nm %in% c("r", "Sigma11", "Sigma22")) {
      tab <- sep_zinb
      key <- if (nm == "Sigma11") "sigma2_b1" else if (nm == "Sigma22") "sigma2_b2" else nm
    } else {
      return(c(NA, NA, NA, NA))
    }
    if (is.null(tab)) return(c(NA, NA, NA, NA))
    i <- match(key, tab$Parameter)
    if (is.na(i)) return(c(NA, NA, NA, NA))
    al <- if ("AL" %in% names(tab)) tab$AL[i] else NA_real_
    c(tab$Bias[i], tab$RMSE[i], al, tab$CP[i])
  }

  lines <- character()
  lines <- c(lines, sprintf("\\multicolumn{11}{@{}l}{\\textbf{$N=%s$}} \\\\", n_label))
  lines <- c(lines, "\\midrule")
  for (nm in rows) {
    j <- joint[joint$Parameter == nm, , drop = FALSE]
    if (!nrow(j)) next
    sep <- sep_lookup(nm)
    is_rho <- grepl("^rho", nm)
    td <- true_d_for(nm)
    jb <- j$Bias[1]; jrm <- j$RMSE[1]; jcp <- j$CP[1]
    jal <- al_joint_lookup(nm)
    # prefer AL from CSV if present
    if ("AL" %in% names(j) && is.finite(j$AL[1])) jal <- j$AL[1]
    sb <- sep[1]; srm <- sep[2]; sal <- sep[3]; scp <- sep[4]
    if (is_rho) {
      lines <- c(lines, sprintf(
        "%s & %s & %s & %s & %s & %s & --- & --- & --- & --- & --- \\\\",
        tex_param(nm), fmt_true(j$True_Value[1], td),
        maybe_bold(fmt_signed(jb, 3), bias_poor(jb)),
        fmt_num(jrm, 3),
        maybe_bold(fmt_num(jcp, 3), cp_poor(jcp)),
        fmt_num(jal, 3)
      ))
    } else {
      lines <- c(lines, sprintf(
        "%s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\",
        tex_param(nm), fmt_true(j$True_Value[1], td),
        maybe_bold(fmt_signed(jb, 3), bias_poor(jb)),
        fmt_num(jrm, 3),
        maybe_bold(fmt_num(jcp, 3), cp_poor(jcp)),
        fmt_num(jal, 3),
        maybe_bold(fmt_signed(sb, 3), bias_poor(sb)),
        fmt_num(srm, 3),
        maybe_bold(fmt_num(scp, 3), cp_poor(scp)),
        fmt_num(sal, 3),
        fmt_vr(vr_lookup(nm))
      ))
    }
  }
  lines
}

load_setting <- function(dir) {
  jfiles <- list.files(dir, pattern = "^joint_param_summary_.*\\.csv$", full.names = TRUE)
  vfiles <- list.files(dir, pattern = "^variance_ratio_.*\\.csv$", full.names = TRUE)
  sfiles <- list.files(dir, pattern = "^slurm-.*\\.out$", full.names = TRUE)
  stopifnot(length(jfiles) == 1L, length(vfiles) == 1L, length(sfiles) == 1L)
  joint <- read.csv(jfiles, stringsAsFactors = FALSE)
  # attach AL from slurm joint block if CSV lacks it
  sep <- parse_slurm_blocks(sfiles)
  if (!is.null(sep$joint)) {
    if (!("AL" %in% names(joint))) {
      joint$AL <- NA_real_
    }
    for (i in seq_len(nrow(sep$joint))) {
      nm <- sep$joint$Parameter[i]
      j <- match(nm, joint$Parameter)
      if (!is.na(j)) joint$AL[j] <- sep$joint$AL[i]
    }
  }
  vr <- read.csv(vfiles, stringsAsFactors = FALSE)
  list(joint = joint, vr = vr, sep = sep, joint_al = sep$joint, path = dir)
}

make_table <- function(settings, ns, scen, label, re_dist = c("mixture", "normal", "mvt"),
                       true_digits = NULL) {
  re_dist <- match.arg(re_dist)
  stopifnot(length(settings) == length(ns), length(settings) >= 1L)
  bold_note <- paste0(
    " Bold Bias entries mark $|\\mathrm{Bias}|\\ge ", sprintf("%.2f", BIAS_POOR),
    "$; bold CP entries mark $\\mathrm{CP}<", sprintf("%.2f", CP_POOR),
    "$ (RMSE and Avg.L are not bolded). Avg.L is the average length of the $95\\%$ credible interval."
  )
  scen_lab <- function(scen) {
    if (scen == 1L) return("Scenario~1 ($\\rho=0.5$")
    if (scen == 2L) return("Scenario~2")
    if (scen == 3L) {
      return("Scenario~3 (shared visit-level covariates with time and Gamma margins")
    }
    stop("scen must be 1, 2, or 3")
  }
  re_lab <- if (identical(re_dist, "normal")) {
    "multivariate random effects"
  } else if (identical(re_dist, "mvt")) {
    "multivariate-$t$ random effects, $df=8$"
  } else {
    "multivariate mixture random effects"
  }
  # Close paren opened in scen_lab for Sc1/Sc3; Sc2 has no open paren.
  if (scen %in% c(1L, 3L)) {
    cap <- paste0("Parameter estimation results of ", scen_lab(scen), "; ", re_lab, ").", bold_note)
  } else {
    cap <- paste0("Parameter estimation results of ", scen_lab(scen), " (", re_lab, ").", bold_note)
  }
  if (is.null(true_digits)) true_digits <- if (scen == 3L) 2L else 1L
  true_sr <- if (identical(re_dist, "mixture")) 3L else 1L
  bodies <- list()
  for (i in seq_along(settings)) {
    s <- settings[[i]]
    bodies[[i]] <- build_block(
      s$joint, s$vr, s$sep$zinb, s$sep$ord, as.character(ns[i]),
      true_digits, true_sr, joint_al = s$joint_al
    )
  }
  body_lines <- bodies[[1]]
  if (length(bodies) >= 2L) {
    for (i in 2:length(bodies)) {
      body_lines <- c(body_lines, "\\midrule", bodies[[i]])
    }
  }
  c(
    sprintf("\\begin{table}[htbp]"),
    "\\centering",
    sprintf("\\caption{%s}", cap),
    sprintf("\\label{%s}", label),
    "\\begin{adjustbox}{max width=\\textwidth,center}",
    "\\scriptsize",
    "\\setlength{\\tabcolsep}{2.8pt}",
    "\\begin{tabular}{lcccccccccc}",
    "\\toprule",
    " & & \\multicolumn{4}{c}{\\textbf{Joint Model}} & \\multicolumn{4}{c}{\\textbf{Separate Models}} & \\textbf{Variance Ratio} \\\\",
    "\\cmidrule(lr){3-6}\\cmidrule(lr){7-10}",
    "\\textbf{Parameter} & \\textbf{True Value} & Bias & RMSE & CP & \\AvgLhdr & Bias & RMSE & CP & \\AvgLhdr & \\textbf{(Gain)} \\\\",
    "\\midrule",
    body_lines,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{adjustbox}",
    "\\end{table}",
    ""
  )
}

has_dir <- function(nm) dir.exists(file.path(ROOT, nm))
use_tab14 <- has_dir("s1_n200_normal") && has_dir("s1_n400_normal") &&
  has_dir("s1_n200_mixture") && has_dir("s1_n400_mixture")
use_n200n400 <- has_dir("scen1_n200") && has_dir("scen1_n400") &&
  has_dir("scen1_n200_normal") && has_dir("scen1_n400_normal")
use_n100n200 <- has_dir("scen1_n100") && has_dir("scen1_n200") &&
  has_dir("scen1_n100_normal") && has_dir("scen1_n200_normal")
has300 <- has_dir("scen1_n300") && has_dir("scen1_n300_normal")

if (use_tab14) {
  s1_a_m <- load_setting(file.path(ROOT, "s1_n200_mixture"))
  s1_b_m <- load_setting(file.path(ROOT, "s1_n400_mixture"))
  s2_a_m <- load_setting(file.path(ROOT, "s2_n200_mixture"))
  s2_b_m <- load_setting(file.path(ROOT, "s2_n400_mixture"))
  s1_a_n <- load_setting(file.path(ROOT, "s1_n200_normal"))
  s1_b_n <- load_setting(file.path(ROOT, "s1_n400_normal"))
  s2_a_n <- load_setting(file.path(ROOT, "s2_n200_normal"))
  s2_b_n <- load_setting(file.path(ROOT, "s2_n400_normal"))
  ns_vec <- c(200L, 400L)
  n_labs <- c("200", "400")
} else if (use_n200n400) {
  s1_a_m <- load_setting(file.path(ROOT, "scen1_n200"))
  s1_b_m <- load_setting(file.path(ROOT, "scen1_n400"))
  s2_a_m <- load_setting(file.path(ROOT, "scen2_n200"))
  s2_b_m <- load_setting(file.path(ROOT, "scen2_n400"))
  s1_a_n <- load_setting(file.path(ROOT, "scen1_n200_normal"))
  s1_b_n <- load_setting(file.path(ROOT, "scen1_n400_normal"))
  s2_a_n <- load_setting(file.path(ROOT, "scen2_n200_normal"))
  s2_b_n <- load_setting(file.path(ROOT, "scen2_n400_normal"))
  ns_vec <- c(200L, 400L)
  n_labs <- c("200", "400")
} else if (use_n100n200) {
  s1_a_m <- load_setting(file.path(ROOT, "scen1_n100"))
  s1_b_m <- load_setting(file.path(ROOT, "scen1_n200"))
  s2_a_m <- load_setting(file.path(ROOT, "scen2_n100"))
  s2_b_m <- load_setting(file.path(ROOT, "scen2_n200"))
  s1_a_n <- load_setting(file.path(ROOT, "scen1_n100_normal"))
  s1_b_n <- load_setting(file.path(ROOT, "scen1_n200_normal"))
  s2_a_n <- load_setting(file.path(ROOT, "scen2_n100_normal"))
  s2_b_n <- load_setting(file.path(ROOT, "scen2_n200_normal"))
  ns_vec <- c(100L, 200L)
  n_labs <- c("100", "200")
  if (has300) {
    s1_c_m <- load_setting(file.path(ROOT, "scen1_n300"))
    s1_c_n <- load_setting(file.path(ROOT, "scen1_n300_normal"))
    s2_c_m <- load_setting(file.path(ROOT, "scen2_n300"))
    s2_c_n <- load_setting(file.path(ROOT, "scen2_n300_normal"))
    ns_vec <- c(100L, 200L, 300L)
    n_labs <- c("100", "200", "300")
  }
} else {
  stop("ROOT must contain tab14 (s*_n{200,400}_*) or scen*_n200+n400 or scen*_n100+n200 folders: ", ROOT)
}

cat("Mode N=", paste(ns_vec, collapse = ","), "\n")
cat("Parsed sep zinb rows mix s2 Na:", NROW(s2_a_m$sep$zinb),
    " ord:", NROW(s2_a_m$sep$ord),
    " joint AL rows:", NROW(s2_a_m$joint_al), "\n")
cat("Sample joint AL alpha1:", {
  j <- s1_a_n$joint
  if ("AL" %in% names(j)) j$AL[j$Parameter == "alpha1"][1] else NA
}, "\n")

pack2 <- function(sa, sb, sc = NULL) {
  if (is.null(sc)) list(sa, sb) else list(sa, sb, sc)
}
s1_c_n <- if (exists("s1_c_n")) s1_c_n else NULL
s2_c_n <- if (exists("s2_c_n")) s2_c_n else NULL
s1_c_m <- if (exists("s1_c_m")) s1_c_m else NULL
s2_c_m <- if (exists("s2_c_m")) s2_c_m else NULL

tex <- c(
  "% Auto-generated by build_sim_tables.R — do not edit by hand",
  "% Bolds Bias if |Bias|>=0.10 and CP if CP<0.90; RMSE/Avg.L never bolded.",
  "% Columns: Bias, RMSE, CP, Avg.L (\\AvgLhdr). Avg.L = average 95% CI length.",
  "% --- Normal random-effects truth ---",
  make_table(pack2(s1_a_n, s1_b_n, s1_c_n), ns_vec, 1L, "tab:sim_scen1_normal", "normal"),
  make_table(pack2(s2_a_n, s2_b_n, s2_c_n), ns_vec, 2L, "tab:sim_scen2_normal", "normal"),
  "% --- Mixture random-effects truth ---",
  make_table(pack2(s1_a_m, s1_b_m, s1_c_m), ns_vec, 1L, "tab:sim_scen1", "mixture"),
  make_table(pack2(s2_a_m, s2_b_m, s2_c_m), ns_vec, 2L, "tab:sim_main", "mixture")
)
dir.create(dirname(OUT_TEX), recursive = TRUE, showWarnings = FALSE)
writeLines(tex, OUT_TEX, useBytes = TRUE)
cat("Wrote", OUT_TEX, "\n")

avg_cp <- function(j) {
  keep <- !grepl("^delta", j$Parameter)
  mean(j$CP[keep], na.rm = TRUE)
}
avg_rmse <- function(j) {
  keep <- !grepl("^delta", j$Parameter)
  mean(j$RMSE[keep], na.rm = TRUE)
}
avg_vr <- function(vr) mean(vr$Variance_Ratio_TrimMean, na.rm = TRUE)
report_one <- function(tag, nlab, s) {
  cat(sprintf("%s N%s: CP=%.3f RMSE=%.3f VRtrim=%.3f\n",
              tag, nlab, avg_cp(s$joint), avg_rmse(s$joint), avg_vr(s$vr)))
}
report_one("Normal S1", n_labs[1], s1_a_n)
report_one("Normal S1", n_labs[2], s1_b_n)
if (!is.null(s1_c_n)) report_one("Normal S1", n_labs[3], s1_c_n)
report_one("Normal S2", n_labs[1], s2_a_n)
report_one("Normal S2", n_labs[2], s2_b_n)
if (!is.null(s2_c_n)) report_one("Normal S2", n_labs[3], s2_c_n)
report_one("Mixture S1", n_labs[1], s1_a_m)
report_one("Mixture S1", n_labs[2], s1_b_m)
if (!is.null(s1_c_m)) report_one("Mixture S1", n_labs[3], s1_c_m)
report_one("Mixture S2", n_labs[1], s2_a_m)
report_one("Mixture S2", n_labs[2], s2_b_m)
if (!is.null(s2_c_m)) report_one("Mixture S2", n_labs[3], s2_c_m)
