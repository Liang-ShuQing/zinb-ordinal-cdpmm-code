#!/usr/bin/env Rscript
# Aggregate one Figure-9 stress cell per row.
#
# Usage:
#   Rscript aggregate_stress_cells.R <job_root>
#
# <job_root> is the directory holding one sub-folder per cell, named
#   n<N>_hd<HD>
# each containing
#   joint_param_summary_scen2_mixture_joint_vs_sep_n<N>_chain10000[_hd<HD>].csv
#
# Output: <job_root>/stress_surface_long_S500_14cells.csv with columns
#   N, HD, ARMSE, ACP, nullMaxBias, note
#
# The two P = 0 cells at N = 200/400 are skipped on purpose: they come from the
# main-text design (same spec, S = 500) and are appended when the figure is drawn.
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) args[[1]] else Sys.getenv("JOB_ROOT", unset = ".")
if (!dir.exists(root)) stop("missing job root: ", root)

keep19 <- c("alpha1", "alpha2", "alpha3", "alpha4",
            "beta1", "beta2", "beta3", "beta4", "beta5",
            "gamma1", "gamma2", "gamma3", "r",
            "Sigma11", "Sigma22", "Sigma33", "rho1", "rho2", "rho3")

# Null (true coefficient 0) covariates added by HD_EXTRA:
# alpha5+, beta6+, gamma4+ (the base design has 4 alphas, 5 betas, 3 gammas).
is_null <- function(p) {
  n <- suppressWarnings(as.integer(gsub("[^0-9]", "", p)))
  (grepl("^alpha", p) & n >= 5) | (grepl("^beta", p) & n >= 6) | (grepl("^gamma", p) & n >= 4)
}

rows <- list()
unfinished <- character(0)
for (N in c(50, 100, 200, 400)) {
  for (HD in c(0, 5, 20, 30)) {
    if (N %in% c(200, 400) && HD == 0) next          # taken from the main design
    tag <- sprintf("n%d_hd%d", N, HD)
    sfx <- if (HD > 0) sprintf("_hd%d", HD) else ""
    f <- file.path(root, tag,
                   sprintf("joint_param_summary_scen2_mixture_joint_vs_sep_n%d_chain10000%s.csv",
                           N, sfx))
    if (!file.exists(f)) { unfinished <- c(unfinished, tag); next }
    d <- read.csv(f, stringsAsFactors = FALSE)
    b <- d[d$Parameter %in% keep19, ]
    nl <- d[is_null(d$Parameter), ]
    if (nrow(b) != 19L) { unfinished <- c(unfinished, tag); next }
    rows[[length(rows) + 1L]] <- data.frame(
      N = N, HD = HD,
      ARMSE = round(mean(b$RMSE), 4),
      ACP = round(mean(b$CP), 4),
      nullMaxBias = if (nrow(nl)) round(max(abs(nl$Bias)), 4) else 0,
      note = "")
  }
}
out <- do.call(rbind, rows)
out <- out[order(out$N, out$HD), ]
rownames(out) <- NULL
f_out <- file.path(root, "stress_surface_long_S500_14cells.csv")
write.csv(out, f_out, row.names = FALSE)
cat("wrote:", f_out, "\n")
print(out, row.names = FALSE)
if (length(unfinished)) cat("\nUNFINISHED:", paste(unfinished, collapse = ", "), "\n")
