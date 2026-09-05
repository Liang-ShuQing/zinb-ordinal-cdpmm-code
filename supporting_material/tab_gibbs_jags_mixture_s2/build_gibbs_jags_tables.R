############################################################
# Build Gibbs (main sim) vs JAGS vs Stan (scheme 1) LaTeX tables.
# Input: 实际数据分析结果/gibbs_jags_sch1_20260822/<cell>/
#   gibbs_joint_param_summary.csv
#   jags_param_summary.csv
#   stan_param_summary.csv  (optional; missing -> ---)
# Four tables: normal/mixture x Scenario 1/2; $N in {200,400}.
# Bolds Bias if |Bias| >= 0.10; CP if CP < 0.90; RMSE never bolded.
############################################################
args <- commandArgs(trailingOnly = TRUE)

script_dir <- tryCatch({
  normalizePath(dirname(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))])), winslash = "/", mustWork = TRUE)
}, error = function(e) normalizePath(".", winslash = "/", mustWork = TRUE))
parent <- dirname(script_dir)
proj_root <- if (basename(script_dir) %in% c("模拟研究工具", "sim_tools")) parent else {
  cand <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
  if (dir.exists(file.path(cand, "paper"))) cand else normalizePath(".", winslash = "/", mustWork = TRUE)
}

DATA_ROOT <- if (length(args) >= 1L && nzchar(args[[1]])) {
  args[[1]]
} else {
  file.path(proj_root, "实际数据分析结果", "gibbs_jags_sch1_20260822")
}
OUT_TEX <- if (length(args) >= 2L && nzchar(args[[2]])) {
  args[[2]]
} else {
  file.path(proj_root, "paper", "supporting_material", "gibbs_jags_tables.tex")
}

BIAS_POOR <- 0.10
CP_POOR <- 0.90

SKIP_PARAMS <- c(
  "delta1", "delta2", "delta3", "delta4",
  "Sigma21", "Sigma31", "Sigma32"
)

tex_param <- function(nm) {
  if (nm == "r") return("$r$")
  if (grepl("^alpha", nm)) return(sprintf("$\\alpha_{%s}$", sub("alpha", "", nm)))
  if (grepl("^beta", nm)) return(sprintf("$\\beta_{%s}$", sub("beta", "", nm)))
  if (grepl("^gamma", nm)) return(sprintf("$\\gamma_{%s}$", sub("gamma", "", nm)))
  if (nm == "Sigma11") return("$\\sigma_1^2$")
  if (nm == "Sigma22") return("$\\sigma_2^2$")
  if (nm == "Sigma33") return("$\\sigma_3^2$")
  if (nm %in% c("rho1", "rho12")) return("$\\rho_{12}$")
  if (nm %in% c("rho2", "rho13")) return("$\\rho_{13}$")
  if (nm %in% c("rho3", "rho23")) return("$\\rho_{23}$")
  nm
}

norm_param <- function(nm) {
  nm <- as.character(nm)
  out <- nm
  out[nm %in% c("rho1", "rho12")] <- "rho12"
  out[nm %in% c("rho2", "rho13")] <- "rho13"
  out[nm %in% c("rho3", "rho23")] <- "rho23"
  out
}

param_sort_key <- function(nm) {
  nm <- norm_param(nm)
  if (grepl("^alpha", nm)) return(100 + as.integer(sub("alpha", "", nm)))
  if (grepl("^beta", nm)) return(200 + as.integer(sub("beta", "", nm)))
  if (grepl("^gamma", nm)) return(300 + as.integer(sub("gamma", "", nm)))
  if (nm == "r") return(400)
  if (nm == "Sigma11") return(500)
  if (nm == "Sigma22") return(501)
  if (nm == "Sigma33") return(502)
  if (nm == "rho12") return(600)
  if (nm == "rho13") return(601)
  if (nm == "rho23") return(602)
  9999
}

fmt_num <- function(x, d = 3L) {
  if (is.na(x) || !is.finite(x)) return("---")
  sprintf(paste0("%.", d, "f"), x)
}
maybe_bold <- function(s, bold) {
  if (!isTRUE(bold) || identical(s, "---")) return(s)
  paste0("\\textbf{", s, "}")
}
bias_poor <- function(x) is.finite(x) && abs(x) >= BIAS_POOR
cp_poor <- function(x) is.finite(x) && x < CP_POOR
fmt_signed <- function(x, d = 3L) {
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

parse_cell <- function(nm) {
  m <- regexec("^s([12])_n([0-9]+)_(normal|mixture)$", nm)
  g <- regmatches(nm, m)[[1]]
  if (length(g) < 4L) return(NULL)
  list(scenario = as.integer(g[2]), n = as.integer(g[3]), re_dist = g[4])
}

read_gibbs <- function(path) {
  if (!file.exists(path)) return(NULL)
  d <- utils::read.csv(path, stringsAsFactors = FALSE)
  if ("Parameter" %in% names(d)) {
    data.frame(
      param = norm_param(d$Parameter),
      True = d$True_Value,
      Bias = d$Bias,
      RMSE = d$RMSE,
      CP = d$CP,
      stringsAsFactors = FALSE
    )
  } else {
    NULL
  }
}

read_method <- function(path) {
  if (!file.exists(path)) return(NULL)
  d <- utils::read.csv(path, stringsAsFactors = FALSE)
  if ("param" %in% names(d)) {
    data.frame(
      param = norm_param(d$param),
      True = d$True,
      Bias = d$Bias,
      RMSE = d$RMSE,
      CP = d$CP,
      stringsAsFactors = FALSE
    )
  } else {
    NULL
  }
}

empty_metrics <- function() {
  data.frame(param = character(0), True = numeric(0), Bias = numeric(0),
             RMSE = numeric(0), CP = numeric(0), stringsAsFactors = FALSE)
}

dirs <- list.dirs(DATA_ROOT, recursive = FALSE, full.names = TRUE)
cells <- list()
for (d in dirs) {
  meta <- parse_cell(basename(d))
  if (is.null(meta)) next
  g <- read_gibbs(file.path(d, "gibbs_joint_param_summary.csv"))
  j <- read_method(file.path(d, "jags_param_summary.csv"))
  s <- read_method(file.path(d, "stan_param_summary.csv"))
  if (is.null(g) || is.null(j)) {
    message("Skip incomplete (need Gibbs+JAGS): ", basename(d))
    next
  }
  meta$dir <- d
  meta$gibbs <- g
  meta$jags <- j
  meta$stan <- if (is.null(s)) empty_metrics() else s
  meta$has_stan <- !is.null(s) && nrow(s) > 0L
  cells[[length(cells) + 1L]] <- meta
}
if (!length(cells)) stop("No complete cells in ", DATA_ROOT)
n_stan <- sum(vapply(cells, function(z) isTRUE(z$has_stan), logical(1)))
message("Cells with Stan: ", n_stan, " / ", length(cells))

method_triplet <- function(df, p) {
  if (is.null(df) || !nrow(df)) {
    return(list(Bias = NA_real_, RMSE = NA_real_, CP = NA_real_))
  }
  r <- df[df$param == p, , drop = FALSE]
  if (!nrow(r)) return(list(Bias = NA_real_, RMSE = NA_real_, CP = NA_real_))
  list(Bias = r$Bias[1], RMSE = r$RMSE[1], CP = r$CP[1])
}

fmt_method_cols <- function(m) {
  sprintf(
    "%s & %s & %s",
    maybe_bold(fmt_signed(m$Bias), bias_poor(m$Bias)),
    fmt_num(m$RMSE),
    maybe_bold(fmt_num(m$CP), cp_poor(m$CP))
  )
}

read_timing <- function(cell_dir, method) {
  fname <- switch(method,
    gibbs = "gibbs_timing_summary.csv",
    jags = "jags_timing_summary.csv",
    stan = "stan_timing_summary.csv",
    NULL
  )
  if (is.null(fname)) return(NULL)
  path <- file.path(cell_dir, fname)
  if (!file.exists(path)) return(NULL)
  t <- utils::read.csv(path, stringsAsFactors = FALSE)
  if (!nrow(t)) return(NULL)
  t[1, , drop = FALSE]
}

timing_triplet <- function(cell_dir) {
  out <- list(
    Gibbs = list(sec_iter = NA_real_, ess_sec = NA_real_),
    JAGS = list(sec_iter = NA_real_, ess_sec = NA_real_),
    Stan = list(sec_iter = NA_real_, ess_sec = NA_real_)
  )
  gt <- read_timing(cell_dir, "gibbs")
  if (!is.null(gt)) {
    chain_len <- as.numeric(gt$Chain[1])
    if (is.finite(chain_len) && chain_len > 0) {
      out$Gibbs$sec_iter <- as.numeric(gt$Mean_Sec[1]) / chain_len
    }
    out$Gibbs$ess_sec <- as.numeric(gt$ESS_per_sec[1])
  }
  jt <- read_timing(cell_dir, "jags")
  if (!is.null(jt)) {
    iter_denom <- as.numeric(jt$Iter[1]) * as.numeric(jt$Chains[1])
    if (is.finite(iter_denom) && iter_denom > 0) {
      out$JAGS$sec_iter <- as.numeric(jt$Mean_Sec[1]) / iter_denom
    }
    out$JAGS$ess_sec <- as.numeric(jt$ESS_per_sec[1])
  }
  st <- read_timing(cell_dir, "stan")
  if (!is.null(st)) {
    iter_denom <- as.numeric(st$Sampling[1]) * as.numeric(st$Chains[1])
    if (is.finite(iter_denom) && iter_denom > 0) {
      out$Stan$sec_iter <- as.numeric(st$Mean_Sec[1]) / iter_denom
    }
    out$Stan$ess_sec <- as.numeric(st$ESS_per_sec[1])
  }
  out
}

fmt_timing_val <- function(x, d = 3L) {
  if (is.na(x) || !is.finite(x)) return("---")
  sprintf(paste0("%.", d, "f"), x)
}

timing_block_tex <- function(cell_dir) {
  tt <- timing_triplet(cell_dir)
  methods <- c("Gibbs", "JAGS", "Stan")
  sec_vals <- vapply(methods, function(m) fmt_timing_val(tt[[m]]$sec_iter), character(1))
  ess_vals <- vapply(methods, function(m) fmt_timing_val(tt[[m]]$ess_sec), character(1))
  paste(
    sprintf(
      "\\multicolumn{2}{@{}l}{\\textit{Sec/iter}} & \\multicolumn{3}{c}{%s} & \\multicolumn{3}{c}{%s} & \\multicolumn{3}{c}{%s} \\\\",
      sec_vals[["Gibbs"]], sec_vals[["JAGS"]], sec_vals[["Stan"]]
    ),
    sprintf(
      "\\multicolumn{2}{@{}l}{\\textit{ESS/sec}} & \\multicolumn{3}{c}{%s} & \\multicolumn{3}{c}{%s} & \\multicolumn{3}{c}{%s} \\\\",
      ess_vals[["Gibbs"]], ess_vals[["JAGS"]], ess_vals[["Stan"]]
    ),
    sep = "\n"
  )
}

block_tex <- function(gibbs, jags, stan) {
  gparams <- setdiff(unique(gibbs$param), SKIP_PARAMS)
  jparams <- setdiff(unique(jags$param), SKIP_PARAMS)
  present <- intersect(gparams, jparams)
  if (!is.null(stan) && nrow(stan)) {
    sparams <- setdiff(unique(stan$param), SKIP_PARAMS)
    present <- intersect(present, sparams)
  }
  present <- present[order(vapply(present, param_sort_key, numeric(1)))]
  lines <- character(0)
  for (p in present) {
    rg <- gibbs[gibbs$param == p, , drop = FALSE]
    if (!nrow(rg)) next
    tv <- rg$True[1]
    mg <- method_triplet(gibbs, p)
    mj <- method_triplet(jags, p)
    ms <- method_triplet(stan, p)
    lines <- c(lines, sprintf(
      "%s & %s & %s & %s & %s \\\\",
      tex_param(p),
      fmt_true(tv),
      fmt_method_cols(mg),
      fmt_method_cols(mj),
      fmt_method_cols(ms)
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
  ncol <- 11L
  for (i in seq_along(want)) {
    z <- want[[i]]
    hdr <- sprintf("\\multicolumn{%d}{@{}l}{\\textbf{$N=%d$}} \\\\", ncol, z$n)
    if (i > 1L) parts <- c(parts, "\\midrule", hdr, "\\midrule")
    else parts <- c(parts, hdr, "\\midrule")
    parts <- c(parts, block_tex(z$gibbs, z$jags, z$stan))
    parts <- c(parts, timing_block_tex(z$dir))
  }
  paste0(
    "\\begin{table}[htbp]\n",
    "\\centering\n",
    "\\caption{", caption, "}\n",
    "\\label{", label, "}\n",
    "\\begin{adjustbox}{max width=\\textwidth,center}\n",
    "\\scriptsize\n",
    "\\setlength{\\tabcolsep}{2.2pt}\n",
    "\\begin{tabular}{lcccccccccc}\n",
    "\\toprule\n",
    " & & \\multicolumn{3}{c}{\\textbf{Gibbs}} & \\multicolumn{3}{c}{\\textbf{JAGS}} & \\multicolumn{3}{c}{\\textbf{Stan}} \\\\\n",
    "\\cmidrule(lr){3-5}\\cmidrule(lr){6-8}\\cmidrule(lr){9-11}\n",
    "\\textbf{Parameter} & \\textbf{True Value} & Bias & RMSE & CP & Bias & RMSE & CP & Bias & RMSE & CP \\\\\n",
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
  "$; bold CP entries mark $\\mathrm{CP}<", sprintf("%.2f", CP_POOR),
  "$ (RMSE is not bolded).",
  " Sec/iter divides mean wall-clock seconds per replicate by the number of Gibbs iterations or JAGS/Stan sampling iterations ($\\times$ chains); ESS/sec is mean effective sample size divided by wall-clock seconds."
)
cap_one <- function(re_lab, scen) {
  paste0(
    "Bias, RMSE, and coverage under Gibbs (main Monte Carlo), JAGS, and Stan ",
    "(", re_lab, "; Scenario~", scen, ", $N\\in\\{200,400\\}$, $S=500$). ",
    "Gibbs: chain $10{,}000$, burn-in $5{,}000$, thin $5$. ",
    "JAGS: two chains, adapt/burn/sampling $500$, thin $1$ ($\\approx 1{,}000$ retained draws). ",
    "Stan: two chains, warmup/sampling $500$ each ($\\approx 1{,}000$ retained draws).",
    note
  )
}

tex <- paste(
  "% Auto-generated by build_gibbs_jags_tables.R",
  "% Gibbs from tab14 main sim; JAGS/Stan from scheme-1 HPC (20260822).",
  "",
  build_table("normal", 1L, "tab:gibbs_jags_normal_s1",
              cap_one("multivariate normal random effects", 1L)),
  build_table("normal", 2L, "tab:gibbs_jags_normal_s2",
              cap_one("multivariate normal random effects", 2L)),
  build_table("mixture", 1L, "tab:gibbs_jags_mixture_s1",
              cap_one("multivariate mixture random effects", 1L)),
  build_table("mixture", 2L, "tab:gibbs_jags_mixture_s2",
              cap_one("multivariate mixture random effects", 2L)),
  sep = "\n"
)

dir.create(dirname(OUT_TEX), recursive = TRUE, showWarnings = FALSE)
con <- file(OUT_TEX, open = "wt", encoding = "UTF-8")
on.exit(close(con), add = TRUE)
writeLines(tex, con, useBytes = FALSE)
message("Wrote ", OUT_TEX)
message("Cells: ", paste(vapply(cells, function(z) {
  sprintf("s%d_n%d_%s%s", z$scenario, z$n, z$re_dist,
          if (isTRUE(z$has_stan)) "" else " [no Stan]")
}, character(1)), collapse = ", "))
