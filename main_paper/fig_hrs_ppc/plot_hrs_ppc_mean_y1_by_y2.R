# HRS PPC: E(HSPNIT | SHLT) from saved MCMC RDS (post-process tool).
args <- commandArgs(trailingOnly = TRUE)
OUT_DIR <- if (length(args) >= 1) args[[1]] else {
  code_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)
  parent <- dirname(code_dir)
  base <- basename(code_dir)
  if (grepl("[\u4e00-\u9fff]", base)) file.path(parent, paste0(base, "\u7ed3\u679c")) else file.path(parent, paste0(base, "_results"))
}
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

proj <- if (length(args) >= 2) args[[2]] else dirname(normalizePath(".", winslash = "/", mustWork = TRUE))
ppc_dir <- file.path(proj, "hpc_jobs", "hrs_ppc_g8k23_xh5_20260915")
rds_path <- if (length(args) >= 3) args[[3]] else {
  file.path(ppc_dir, "combined_samples_slim.rds")
}
data_path <- if (length(args) >= 4) args[[4]] else {
  file.path(ppc_dir, "Data_model_complete.csv")
}
idx_path <- if (length(args) >= 5) args[[5]] else {
  file.path(ppc_dir, "hrs_ppc_reps_summary.rds")
}

if (!requireNamespace("mvtnorm", quietly = TRUE)) stop("Need package mvtnorm")
inv_logit <- function(x) stats::plogis(x)

cat("Loading RDS...\n")
fit <- readRDS(rds_path)
cs <- fit$combined_samples
stopifnot(!is.null(cs$Alpha), !is.null(cs$Beta), !is.null(cs$Gamma),
          !is.null(cs$Delta), !is.null(cs$Sigma), !is.null(cs$R))

covariate_vars <- c("SMOKEV", "HIBP", "DIAB", "LUNG")
dat_raw <- read.csv(data_path)
if ("wave" %in% colnames(dat_raw)) {
  dat_raw <- dat_raw[order(dat_raw$id, dat_raw$wave), ]
} else {
  dat_raw <- dat_raw[order(dat_raw$id), ]
}
dat_raw$SHLT_num <- as.numeric(dat_raw$SHLT)
dat_raw$SHLT_recoded <- NA_integer_
dat_raw$SHLT_recoded[dat_raw$SHLT_num %in% c(1, 2)] <- 1L
dat_raw$SHLT_recoded[dat_raw$SHLT_num == 3] <- 2L
dat_raw$SHLT_recoded[dat_raw$SHLT_num %in% c(4, 5)] <- 3L
keep_bin <- !is.na(dat_raw$HIBP) & dat_raw$HIBP <= 1 &
  !is.na(dat_raw$DIAB) & dat_raw$DIAB <= 1 &
  !is.na(dat_raw$LUNG) & dat_raw$LUNG <= 1
dat_raw <- dat_raw[keep_bin, ]
cc <- stats::complete.cases(dat_raw[, c("HSPNIT", "SHLT_recoded", covariate_vars)])
dat <- dat_raw[cc, ]

y1 <- dat$HSPNIT
y2 <- as.integer(dat$SHLT_recoded)
ID <- dat$id
C <- length(unique(y2))
X <- cbind(Intercept = 1, as.matrix(dat[, covariate_vars, drop = FALSE]))
id <- as.numeric(factor(ID))
ord_rows <- order(id)
id <- id[ord_rows]; y1 <- y1[ord_rows]; y2 <- y2[ord_rows]
X <- X[ord_rows, , drop = FALSE]
n <- length(unique(id)); N <- length(id)
nis <- as.numeric(table(id))
cat(sprintf("Data: n=%d N=%d C=%d\n", n, N, C))
stopifnot(n == fit$settings$n || TRUE) # warn only
if (!identical(as.integer(fit$settings$n), as.integer(n))) {
  cat(sprintf("NOTE: RDS n=%s vs data n=%s (using data nis for replicates)\n",
              as.character(fit$settings$n), as.character(n)))
}

S_all <- nrow(cs$Alpha)
idx_saved <- readRDS(idx_path)$idx_ppc
if (is.null(idx_saved) || any(idx_saved < 1L) || any(idx_saved > S_all)) {
  stop("idx_ppc in ", idx_path, " is missing or out of range for S_all=", S_all)
}
idx_ppc <- as.integer(idx_saved)
S_ppc <- length(idx_ppc)
cat(sprintf("PPC draws: %d / %d (idx from Fig. 11 summary)\n", S_ppc, S_all))

levels_y2 <- sort(unique(y2))
obs_mean <- vapply(levels_y2, function(c) mean(y1[y2 == c]), numeric(1))
rep_mean <- matrix(NA_real_, S_ppc, length(levels_y2))

for (ii in seq_len(S_ppc)) {
  s <- idx_ppc[ii]
  alpha_s <- cs$Alpha[s, ]
  beta_s  <- cs$Beta[s, ]
  gamma_s <- cs$Gamma[s, ]
  delta_s <- cs$Delta[s, ]
  r_s     <- cs$R[s]
  Sigma_s <- matrix(cs$Sigma[s, ], 3L, 3L)
  b_s <- tryCatch(
    mvtnorm::rmvnorm(n, sigma = Sigma_s),
    error = function(e) mvtnorm::rmvnorm(n, sigma = Sigma_s + diag(1e-4, 3L))
  )
  b1s <- rep(b_s[, 1], times = nis)
  b2s <- rep(b_s[, 2], times = nis)
  b3s <- rep(b_s[, 3], times = nis)
  pi_i <- inv_logit(as.numeric(X %*% alpha_s + b1s))
  u_i  <- stats::rbinom(N, 1L, pmin(pmax(pi_i, 1e-10), 1 - 1e-10))
  phi_i <- inv_logit(as.numeric(X %*% beta_s + b2s))
  mu_i  <- pmax(r_s * phi_i / (1 - phi_i), 1e-10)
  y1rep <- as.integer(u_i * stats::rnbinom(N, size = r_s, mu = mu_i))
  eta_ord <- as.numeric(X %*% gamma_s + b3s)
  l_rep <- eta_ord + stats::rlogis(N)
  y2rep <- as.integer(cut(l_rep, breaks = c(-Inf, delta_s, Inf),
                          labels = FALSE, include.lowest = TRUE))
  bad <- which(is.na(y2rep))
  if (length(bad)) y2rep[bad] <- ifelse(l_rep[bad] <= delta_s[1], 1L, C)
  for (j in seq_along(levels_y2)) {
    idx <- which(y2rep == levels_y2[j])
    rep_mean[ii, j] <- if (length(idx)) mean(y1rep[idx]) else NA_real_
  }
  if (ii %% 100L == 0L || ii == S_ppc) cat(sprintf("  replicate %d / %d\n", ii, S_ppc))
}

pred_mean <- colMeans(rep_mean, na.rm = TRUE)
lo_mean <- apply(rep_mean, 2, stats::quantile, probs = 0.025, na.rm = TRUE)
hi_mean <- apply(rep_mean, 2, stats::quantile, probs = 0.975, na.rm = TRUE)

lab <- c("Favorable", "Moderate", "Unfavorable")
if (length(levels_y2) != 3L) lab <- paste0("SHLT=", levels_y2)

summ <- data.frame(
  SHLT = lab,
  level = levels_y2,
  Obs = as.numeric(obs_mean),
  Pred = as.numeric(pred_mean),
  Lo = as.numeric(lo_mean),
  Hi = as.numeric(hi_mean)
)
write.csv(summ, file.path(OUT_DIR, "hrs_ppc_mean_y1_by_y2.csv"), row.names = FALSE)
saveRDS(list(summary = summ, rep_mean = rep_mean, idx_ppc = idx_ppc,
             rds_path = rds_path, data_path = data_path,
             settings = list(S_ppc = S_ppc, S_all = S_all, n = n, N = N, C = C)),
        file.path(OUT_DIR, "hrs_ppc_mean_y1_by_y2.rds"))

col_obs <- "#0072B2"
col_pred <- "#E69F00"
ylim_top <- max(c(obs_mean, hi_mean), na.rm = TRUE) * 1.12
out_pdf <- file.path(OUT_DIR, "hrs_ppc_mean_y1_by_y2.pdf")
grDevices::pdf(out_pdf, width = 6.2, height = 4.6)
op <- graphics::par(mar = c(4.2, 4.2, 1.2, 1.2))
bp <- graphics::barplot(
  obs_mean, names.arg = lab, col = col_obs, border = NA,
  ylim = c(0, ylim_top),
  xlab = "Self-rated health (SHLT)",
  ylab = expression(Mean ~ of ~ HSPNIT)
)
graphics::arrows(x0 = bp, y0 = lo_mean, y1 = hi_mean,
                 angle = 90, code = 3, length = 0.06, col = col_pred, lwd = 1.5)
graphics::points(bp, pred_mean, pch = 16, col = col_pred, cex = 1.2)
graphics::par(op)
grDevices::dev.off()
cat("Wrote", out_pdf, "\n")
print(summ)