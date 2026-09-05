############################################################
# Scenario-2 illustrative 3-chain diagnostics (CDPMM joint)
# Plots ONLY alpha, beta, gamma, r:
#   - three_chain5000.pdf
#   - GelmanRubinConvergencePlot5000.pdf
# Uses helpers from prog1_cdpmm_joint_vs_separate.R (SKIP_MAIN_SIM).
############################################################

args <- commandArgs(trailingOnly = TRUE)
OUT_DIR <- if (length(args) >= 1) args[[1]] else {
  code_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)
  parent <- dirname(code_dir)
  base <- basename(code_dir)
  if (grepl("[\u4e00-\u9fff]", base)) {
    file.path(parent, paste0(base, "结果"))
  } else {
    file.path(parent, paste0(base, "_results"))
  }
}
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SKIP_MAIN_SIM <- TRUE
sim_candidates <- c("prog1_cdpmm_joint_vs_separate.R", "code.R", "simulation_code_hpc.R")
sim_file <- sim_candidates[file.exists(sim_candidates)][1]
if (is.na(sim_file)) stop("Cannot find prog1_cdpmm_joint_vs_separate.R")
source(sim_file, encoding = "UTF-8", local = FALSE)

library(coda)
library(ggplot2)
library(parallel)
library(doParallel)
library(foreach)

# Diagnostics settings (3 chains; burn=0 for EPSR trajectory from iter 1)
# Default CHAIN=10000 to align with main Monte Carlo length; override via env.
n_chains <- 3L
chain_diag <- as.integer(Sys.getenv("CHAIN", unset = "10000"))
burn_diag <- 0L
thin_diag <- 1L
seed_data <- as.integer(Sys.getenv("SEED_DATA", unset = "2025"))
scenario_diag <- 2L
re_dist_diag <- Sys.getenv("RE_DIST", unset = "mixture")
if (!nzchar(re_dist_diag)) re_dist_diag <- "mixture"
epsr_step <- as.integer(Sys.getenv("EPSR_STEP", unset = "25"))
epsr_min <- as.integer(Sys.getenv("EPSR_MIN", unset = "50"))
if (!is.finite(chain_diag) || chain_diag < 100L) chain_diag <- 10000L

# Optional: also copy figures into paper/figures/
paper_fig_dir <- Sys.getenv("PAPER_FIG_DIR", unset = "")
if (!nzchar(paper_fig_dir)) {
  cand <- file.path(dirname(normalizePath(".", winslash = "/", mustWork = TRUE)),
                    "paper", "figures")
  if (dir.exists(cand)) paper_fig_dir <- cand
}

cat("=== 3-chain EPSR / trace diagnostics ===\n")
cat(sprintf("OUT_DIR=%s\n", OUT_DIR))
cat(sprintf("scenario=%d, re_dist=%s, n=%d, chain=%d, burn=%d, thin=%d\n",
            scenario_diag, re_dist_diag, n, chain_diag, burn_diag, thin_diag))

############################################################
# 1. One representative Scenario-2 dataset
############################################################
dat <- generate_data(
  n = n,
  nis = nis_fixed,
  seed = seed_data,
  random_nis = random_nis,
  nis_range = nis_range,
  scenario = scenario_diag,
  C = C,
  re_dist = re_dist_diag,
  mvt_df = mvt_df
)

p_ord <- ncol(dat$X_ordinal)
p_zero <- ncol(dat$X_zero)
p_count <- ncol(dat$X_count)
n_subj <- dat$n
C_use <- dat$C

# Notation: alpha=zero (p_zero), beta=count (p_count), gamma=ordinal (p_ord)
cat(sprintf("Data: N=%d, n=%d | p_alpha=%d, p_beta=%d, p_gamma=%d\n",
            dat$N, n_subj, p_zero, p_count, p_ord))

############################################################
# 2. Dispersed fixed initial values (3 chains)
############################################################
make_init <- function(scale_fe, r0, sigma_diag, rho0, b_level) {
  Sig <- matrix(rho0, 3, 3)
  diag(Sig) <- sigma_diag
  list(
    alpha = scale_fe * c(0.1, -0.1, 0.05, 0.08, -0.05)[seq_len(p_zero)],
    beta  = scale_fe * c(0.1, 0.05, -0.1, 0.05, -0.08, 0.06)[seq_len(p_count)],
    gamma = scale_fe * c(0.1, -0.1, 0.1, -0.05)[seq_len(p_ord)],
    delta = seq(0, 3, length.out = C_use - 1L),
    r = r0,
    Sigma = Sig,
    b = matrix(b_level, n_subj, 3L)
  )
}

fixed_inits <- list(
  make_init(scale_fe = 1.0, r0 = 1.2, sigma_diag = 0.8, rho0 = 0.2, b_level = 0.1),
  make_init(scale_fe = 3.0, r0 = 2.0, sigma_diag = 1.2, rho0 = 0.4, b_level = 0.3),
  make_init(scale_fe = 6.0, r0 = 3.0, sigma_diag = 1.5, rho0 = 0.55, b_level = 0.5)
)

############################################################
# 3. Run 3 chains (parallel when possible)
############################################################
run_one_chain <- function(chain_id) {
  set.seed(2025L + as.integer(chain_id))
  cat(sprintf("\n--- Chain %d start ---\n", chain_id))
  t0 <- proc.time()[3]
  fit <- fit_joint_model(
    dat = dat,
    chain = chain_diag,
    burn = burn_diag,
    thin = thin_diag,
    delta_min = delta_min,
    delta_max = delta_max,
    G = G_mix,
    init = fixed_inits[[chain_id]]
  )
  elapsed <- proc.time()[3] - t0
  cat(sprintf("--- Chain %d done (%.1f min) ---\n", chain_id, elapsed / 60))
  list(
    chain_id = chain_id,
    Alpha = fit$alpha_samples,
    Beta  = fit$beta_samples,
    Gamma = fit$gamma_samples,
    R     = as.numeric(fit$r_samples),
    elapsed_min = elapsed / 60
  )
}

n_cores_env <- suppressWarnings(as.integer(Sys.getenv(
  c("N_CORES", "SLURM_NTASKS_PER_NODE", "SLURM_CPUS_PER_TASK")
)))
n_cores_env <- n_cores_env[is.finite(n_cores_env) & n_cores_env >= 1L]
n_cores <- if (length(n_cores_env)) n_cores_env[[1]] else {
  max(1L, parallel::detectCores(logical = TRUE) - 1L)
}
n_cores <- min(n_chains, as.integer(n_cores))
cat(sprintf("Running %d chains with %d fork workers...\n", n_chains, n_cores))

if (n_cores <= 1L) {
  chain_list <- lapply(seq_len(n_chains), run_one_chain)
} else {
  # Fork (Linux/HPC): no PSOCK makeCluster
  doParallel::registerDoParallel(cores = n_cores)
  on.exit(try(doParallel::stopImplicitCluster(), silent = TRUE), add = TRUE)
  chain_list <- foreach::foreach(
    chain_id = seq_len(n_chains),
    .packages = c("BayesLogit", "mvtnorm", "MCMCpack", "truncnorm", "coda", "loo")
  ) %dopar% {
    run_one_chain(chain_id)
  }
  doParallel::stopImplicitCluster()
  on.exit(NULL)
}

saveRDS(
  list(dat = dat, chains = chain_list,
       settings = list(scenario = scenario_diag, re_dist = re_dist_diag,
                       chain = chain_diag, burn = burn_diag, thin = thin_diag,
                       seed_data = seed_data)),
  file.path(OUT_DIR, "three_chain_diag_raw.rds")
)
cat("Saved three_chain_diag_raw.rds\n")

############################################################
# 4. Parameter specs (alpha / beta / gamma / r only)
############################################################
param_specs <- list()
for (i in seq_len(p_zero)) {
  param_specs[[length(param_specs) + 1L]] <- list(
    key = paste0("alpha", i), type = "Alpha", idx = i,
    label_chr = sprintf("alpha[%d]", i), family = "alpha"
  )
}
for (i in seq_len(p_count)) {
  param_specs[[length(param_specs) + 1L]] <- list(
    key = paste0("beta", i), type = "Beta", idx = i,
    label_chr = sprintf("beta[%d]", i), family = "beta"
  )
}
for (i in seq_len(p_ord)) {
  param_specs[[length(param_specs) + 1L]] <- list(
    key = paste0("gamma", i), type = "Gamma", idx = i,
    label_chr = sprintf("gamma[%d]", i), family = "gamma"
  )
}
param_specs[[length(param_specs) + 1L]] <- list(
  key = "r", type = "R", idx = 1L,
  label_chr = "r", family = "r"
)

extract_vec <- function(ch, type, idx) {
  if (type == "Alpha") return(ch$Alpha[, idx])
  if (type == "Beta")  return(ch$Beta[, idx])
  if (type == "Gamma") return(ch$Gamma[, idx])
  if (type == "R")     return(as.numeric(ch$R))
  stop("unknown type")
}

############################################################
# 5. Three-chain trace plots
############################################################
cat("\nDrawing three-chain trace plots...\n")

trace_df <- do.call(rbind, lapply(param_specs, function(spec) {
  do.call(rbind, lapply(seq_along(chain_list), function(k) {
    vals <- extract_vec(chain_list[[k]], spec$type, spec$idx)
    data.frame(
      iteration = seq_along(vals),
      value = vals,
      chain = factor(paste0("Chain ", k), levels = paste0("Chain ", seq_len(n_chains))),
      param_key = spec$key,
      stringsAsFactors = FALSE
    )
  }))
}))

labeller_chr <- setNames(
  vapply(param_specs, `[[`, "", "label_chr"),
  vapply(param_specs, `[[`, "", "key")
)
trace_df$param_key <- factor(trace_df$param_key, levels = names(labeller_chr))

n_params <- length(param_specs)
ncol_tr <- 4L
nrow_tr <- ceiling(n_params / ncol_tr)

p_trace <- ggplot(trace_df, aes(x = iteration, y = value, color = chain)) +
  geom_line(linewidth = 0.35, alpha = 0.85) +
  facet_wrap(~ param_key, scales = "free_y", ncol = ncol_tr,
             labeller = as_labeller(labeller_chr, label_parsed)) +
  scale_color_manual(values = c("Chain 1" = "blue", "Chain 2" = "green3", "Chain 3" = "red")) +
  labs(x = "Iteration", y = "Value", color = NULL) +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 10),
    panel.grid.minor = element_blank(),
    plot.margin = margin(8, 8, 8, 8)
  )

trace_file <- file.path(OUT_DIR, sprintf("three_chain%d.pdf", chain_diag))
ggsave(trace_file, p_trace,
       width = 11, height = max(6.5, 2.0 * nrow_tr),
       device = grDevices::pdf)
cat("Wrote:", trace_file, "\n")

############################################################
# 6. EPSR (R-hat) trajectories
############################################################
cat("\nComputing EPSR trajectories...\n")

create_mcmc_list_one <- function(type, idx) {
  coda::as.mcmc.list(lapply(chain_list, function(ch) {
    coda::as.mcmc(matrix(extract_vec(ch, type, idx), ncol = 1L))
  }))
}

calculate_rhat_evolution <- function(type, idx, min_length = epsr_min, step_size = epsr_step) {
  mlist <- create_mcmc_list_one(type, idx)
  n_samples <- nrow(as.matrix(mlist[[1]]))
  seq_lengths <- seq(min_length, n_samples, by = step_size)
  if (tail(seq_lengths, 1) < n_samples) seq_lengths <- c(seq_lengths, n_samples)

  out <- data.frame(Chain_Length = integer(), EPSR = numeric())
  for (L in seq_lengths) {
    trunc <- coda::as.mcmc.list(lapply(mlist, function(ch) {
      coda::as.mcmc(as.matrix(ch)[seq_len(L), , drop = FALSE])
    }))
    gd <- tryCatch(
      coda::gelman.diag(trunc, autoburnin = FALSE, multivariate = FALSE),
      error = function(e) NULL
    )
    if (is.null(gd)) next
    psrf <- gd$psrf
    if (is.null(dim(psrf))) psrf <- matrix(psrf, nrow = 1)
    rh <- as.numeric(psrf[1, 1])
    if (is.finite(rh)) {
      out <- rbind(out, data.frame(Chain_Length = L, EPSR = rh))
    }
  }
  out
}

epsr_df <- do.call(rbind, lapply(param_specs, function(spec) {
  tmp <- calculate_rhat_evolution(spec$type, spec$idx)
  if (nrow(tmp) == 0) return(NULL)
  tmp$Parameter <- spec$key
  tmp
}))
epsr_df$Parameter <- factor(epsr_df$Parameter, levels = vapply(param_specs, `[[`, "", "key"))

param_labels <- setNames(
  lapply(vapply(param_specs, `[[`, "", "label_chr"), function(s) parse(text = s)[[1]]),
  vapply(param_specs, `[[`, "", "key")
)

n_sel <- length(param_specs)
cols <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
  "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62",
  "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#B3B3B3"
)
ltys <- rep(c("solid", "dashed", "dotted", "dotdash"), length.out = n_sel)
color_scale <- setNames(cols[seq_len(n_sel)], levels(epsr_df$Parameter))
linetype_scale <- setNames(ltys, levels(epsr_df$Parameter))

y_max <- max(2.0, ceiling(max(epsr_df$EPSR, na.rm = TRUE) * 10) / 10)
y_max <- min(y_max, 5)

p_epsr <- ggplot(epsr_df, aes(x = Chain_Length, y = EPSR,
                              color = Parameter, linetype = Parameter)) +
  geom_line(linewidth = 0.65, alpha = 0.9) +
  geom_hline(yintercept = 1.0, linetype = "solid", color = "darkgreen", linewidth = 0.65) +
  scale_color_manual(values = color_scale, labels = param_labels) +
  scale_linetype_manual(values = linetype_scale, labels = param_labels) +
  scale_x_continuous(breaks = seq(0, chain_diag, by = 1000),
                     limits = c(0, chain_diag * 1.02),
                     expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(limits = c(0.95, y_max),
                     breaks = seq(1.0, y_max, by = 0.2)) +
  labs(x = "Iteration", y = "EPSR (R-hat)", color = "Parameter", linetype = "Parameter") +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 9),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

epsr_file <- file.path(
  OUT_DIR, sprintf("GelmanRubinConvergencePlot%d.pdf", chain_diag)
)
ggsave(epsr_file, p_epsr, width = 12, height = 7, device = grDevices::pdf)
cat("Wrote:", epsr_file, "\n")

# Final EPSR table
final_epsr <- do.call(rbind, lapply(param_specs, function(spec) {
  mlist <- create_mcmc_list_one(spec$type, spec$idx)
  gd <- tryCatch(
    coda::gelman.diag(mlist, autoburnin = FALSE, multivariate = FALSE),
    error = function(e) NULL
  )
  rh <- if (is.null(gd)) NA_real_ else as.numeric(gd$psrf[1, 1])
  data.frame(Parameter = spec$key, EPSR = round(rh, 4),
             Converged_1.2 = is.finite(rh) && rh < 1.2)
}))
write.csv(final_epsr, file.path(OUT_DIR, "EPSR_final_alpha_beta_gamma_r.csv"),
          row.names = FALSE)
print(final_epsr)
cat(sprintf("EPSR < 1.2: %d / %d\n",
            sum(final_epsr$Converged_1.2, na.rm = TRUE), nrow(final_epsr)))

############################################################
# 7. Copy into paper/figures if available
############################################################
if (nzchar(paper_fig_dir) && dir.exists(paper_fig_dir)) {
  file.copy(trace_file, file.path(paper_fig_dir, basename(trace_file)), overwrite = TRUE)
  file.copy(epsr_file, file.path(paper_fig_dir, basename(epsr_file)), overwrite = TRUE)
  cat("Copied PDFs to:", paper_fig_dir, "\n")
}

cat("\n=== Done ===\n")
cat("Trace:", trace_file, "\n")
cat("EPSR: ", epsr_file, "\n")

