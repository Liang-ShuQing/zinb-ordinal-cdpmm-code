###############################################################################
# Posterior Predictive Checks (PPC) for the joint ZINB--ordinal CDPMM model
# 布局对齐 参考代码/ppc_plots.R（上一篇 binary--ZINB）
#
#   本篇记号: y1 = ZINB 计数, y2 = 有序多分类 (C 类)
#
#   Fig.4  计数分布: 真实 y1 各取值比例(蓝柱) vs 后验预测比例(红点+2.5/97.5误差棒)
#   Fig.5  cor(y1, y2) PPC: 观测相关(蓝竖线) vs 后验预测相关分布(红直方图)
#   Fig.6  cor(I(y1>0), y2) PPC: 同上
#   Fig.7  按 y2 分组的 y1 均值 PPC: 观测均值(蓝) vs 后验预测 2.5/97.5 区间(红)
#   合并图  Fig.4/5/6 纵向拼成一页 (a)(b)(c); 另输出 Fig.4+6 (a)(b) 论文用简版
#
# 说明:
#   - source(prog1) 复用 generate_data / fit_joint_model（SKIP_MAIN_SIM）
#   - 不同零膨胀率 = 改变零膨胀子模型【截距】zi_intercept（越小 -> 在险概率越低 -> 零率越高）
#   - 后验预测复制: 从后验样本抽 Sigma 生成 RE, 再联合生成 (y1_rep, y2_rep)
#   - 默认: Scenario 2 + mixture RE；可用 Env 覆盖
#   - Env: CHAIN, BURN, THIN, N, SCENARIOS, ZI_INTERCEPTS, K_MAX, RE_DIST, OUT_DIR
#   - CLI: 第一个参数为 OUT_DIR（默认并列 模拟研究结果/）
###############################################################################

args <- commandArgs(trailingOnly = TRUE)
OUT_DIR <- if (length(args) >= 1L && nzchar(args[[1]])) {
  args[[1]]
} else {
  env_out <- Sys.getenv("OUT_DIR", unset = "")
  if (nzchar(env_out)) {
    env_out
  } else {
    code_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)
    parent <- dirname(code_dir)
    base <- basename(code_dir)
    if (grepl("[\u4e00-\u9fff]", base)) {
      file.path(parent, paste0(base, "结果"))
    } else {
      file.path(parent, paste0(base, "_results"))
    }
  }
}
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SKIP_MAIN_SIM <- TRUE
sim_candidates <- c("prog1_cdpmm_joint_vs_separate.R", "code.R")
sim_file <- sim_candidates[file.exists(sim_candidates)][1]
if (is.na(sim_file)) stop("Cannot find prog1_cdpmm_joint_vs_separate.R")
source(sim_file, encoding = "UTF-8", local = FALSE)

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# =============================================================================
# 0. 可修改参数（可用环境变量覆盖；本地冒烟可调小 CHAIN）
# =============================================================================
# 候选组: REP_ID=1..K -> SEED_BASE = 100 + 1000*(REP_ID-1)
# 也可直接设 SEED_BASE / INITIAL_SEED
rep_id <- suppressWarnings(as.integer(Sys.getenv("REP_ID", unset = "")))
seed_base_env <- suppressWarnings(as.integer(Sys.getenv("SEED_BASE", unset = "")))
init_seed_env <- suppressWarnings(as.integer(Sys.getenv("INITIAL_SEED", unset = "")))
if (length(seed_base_env) == 1L && is.finite(seed_base_env)) {
  initial_seed <- as.integer(seed_base_env)
} else if (length(init_seed_env) == 1L && is.finite(init_seed_env)) {
  initial_seed <- as.integer(init_seed_env)
} else if (length(rep_id) == 1L && is.finite(rep_id) && rep_id >= 1L) {
  initial_seed <- as.integer(100L + 1000L * (rep_id - 1L))
} else {
  initial_seed <- 100L
  rep_id <- NA_integer_
}
n_ppc <- as.integer(Sys.getenv("N", unset = "100"))
if (!is.finite(n_ppc) || n_ppc < 1L) n_ppc <- 100L

random_nis <- TRUE
nis_range  <- 1:20
nis_fixed  <- 10L
C_ord      <- 5L
re_dist_ppc <- Sys.getenv("RE_DIST", unset = "mixture")

scen_env <- Sys.getenv("SCENARIOS", unset = "2")
scenarios <- as.integer(strsplit(scen_env, "[,; ]+")[[1]])
scenarios <- scenarios[is.finite(scenarios)]
if (length(scenarios) < 1L) scenarios <- 2L

zi_env <- Sys.getenv("ZI_INTERCEPTS", unset = "2.0,0.5,-1.0,-2.0")
zi_intercepts <- as.numeric(strsplit(zi_env, "[,; ]+")[[1]])
zi_intercepts <- zi_intercepts[is.finite(zi_intercepts)]
if (length(zi_intercepts) < 1L) zi_intercepts <- c(2.0, 0.5, -1.0, -2.0)

chain_ppc <- as.integer(Sys.getenv("CHAIN", unset = "5000"))
burn_ppc  <- as.integer(Sys.getenv("BURN", unset = "2000"))
thin_ppc  <- as.integer(Sys.getenv("THIN", unset = "5"))
if (!is.finite(chain_ppc) || chain_ppc < 10L) chain_ppc <- 5000L
if (!is.finite(burn_ppc) || burn_ppc < 0L) burn_ppc <- 2000L
if (!is.finite(thin_ppc) || thin_ppc < 1L) thin_ppc <- 5L

K_max <- as.integer(Sys.getenv("K_MAX", unset = "15"))
if (!is.finite(K_max) || K_max < 1L) K_max <- 15L

safe_mean <- function(x) {
  if (length(x) == 0) NA_real_ else mean(x, na.rm = TRUE)
}

# =============================================================================
# 1. 后验预测复制: 联合生成 y1_rep (ZINB) 与 y2_rep (ordinal)
# =============================================================================
generate_joint_ppc_reps <- function(fit, dat) {
  S <- nrow(fit$alpha_samples)
  N <- dat$N
  n <- dat$n
  nis <- dat$nis
  C <- dat$C
  Xz <- dat$X_zero
  Xc <- dat$X_count
  Xo <- dat$X_ordinal

  Y1rep <- matrix(NA_integer_, S, N)
  Y2rep <- matrix(NA_integer_, S, N)

  for (s in seq_len(S)) {
    alpha_s <- fit$alpha_samples[s, ]
    beta_s  <- fit$beta_samples[s, ]
    gamma_s <- fit$gamma_samples[s, ]
    delta_s <- fit$delta_samples[s, ]
    r_s     <- fit$r_samples[s]
    Sigma_s <- matrix(fit$Sigma_samples[s, ], 3, 3)
    # 用后验隐含边际协方差抽样 RE（与参考代码对 IW/Sigma 的用法一致）
    b_s <- tryCatch(
      rmvnorm(n, sigma = Sigma_s),
      error = function(e) rmvnorm(n, sigma = Sigma_s + diag(1e-4, 3))
    )

    b1 <- rep(b_s[, 1], times = nis)
    b2 <- rep(b_s[, 2], times = nis)
    b3 <- rep(b_s[, 3], times = nis)

    pi_i <- inv_logit(as.numeric(Xz %*% alpha_s + b1))
    u_i  <- rbinom(N, 1, pmin(pmax(pi_i, 1e-10), 1 - 1e-10))
    phi_i <- inv_logit(as.numeric(Xc %*% beta_s + b2))
    mu_i  <- pmax(r_s * phi_i / (1 - phi_i), 1e-10)
    Y1rep[s, ] <- as.integer(u_i * rnbinom(N, size = r_s, mu = mu_i))

    eta_ord <- as.numeric(Xo %*% gamma_s + b3)
    l_rep <- eta_ord + stats::rlogis(N)
    Y2rep[s, ] <- as.integer(cut(
      l_rep,
      breaks = c(-Inf, delta_s, Inf),
      labels = FALSE,
      include.lowest = TRUE
    ))
    # 数值保护: cut 偶发 NA 时落到两端类别
    bad <- which(is.na(Y2rep[s, ]))
    if (length(bad) > 0L) {
      Y2rep[s, bad] <- ifelse(l_rep[bad] <= delta_s[1], 1L, C)
    }
  }
  list(Y1rep = Y1rep, Y2rep = Y2rep)
}

compute_association_ppc <- function(y1, y2, Y1rep, Y2rep, C) {
  S <- nrow(Y1rep)
  cats <- seq_len(C)

  rep_cor_y1_y2 <- vapply(seq_len(S), function(s) {
    safe_cor(Y1rep[s, ], Y2rep[s, ])
  }, numeric(1))
  rep_cor_y1pos_y2 <- vapply(seq_len(S), function(s) {
    safe_cor(as.numeric(Y1rep[s, ] > 0), Y2rep[s, ])
  }, numeric(1))

  mean_mat <- matrix(NA_real_, S, C)
  for (s in seq_len(S)) {
    for (c in cats) {
      mean_mat[s, c] <- safe_mean(Y1rep[s, Y2rep[s, ] == c])
    }
  }

  obs_mean <- vapply(cats, function(c) safe_mean(y1[y2 == c]), numeric(1))

  list(
    obs = list(
      cor_y1_y2     = safe_cor(y1, y2),
      cor_y1pos_y2  = safe_cor(as.numeric(y1 > 0), y2),
      mean_y1_by_y2 = obs_mean
    ),
    rep = list(
      cor_y1_y2     = rep_cor_y1_y2,
      cor_y1pos_y2  = rep_cor_y1pos_y2,
      mean_y1_by_y2 = mean_mat
    )
  )
}

summarize_rep <- function(x) {
  c(
    pred = mean(x, na.rm = TRUE),
    lo   = as.numeric(stats::quantile(x, 0.025, na.rm = TRUE)),
    hi   = as.numeric(stats::quantile(x, 0.975, na.rm = TRUE))
  )
}

# =============================================================================
# 2. 对每个情景 × 零膨胀率: 生成数据 -> 拟合 -> 后验预测 -> 统计量
# =============================================================================
dist_df <- list()
cor_rep_df <- list()
mean_rep_df <- list()
assoc_summary_df <- list()
idx <- 1L
ks <- 0:K_max

cat(sprintf(
  "PPC settings: n=%d, scenarios=%s, zi_intercepts=%s, chain=%d, burn=%d, thin=%d, re_dist=%s\n",
  n_ppc, paste(scenarios, collapse = ","),
  paste(zi_intercepts, collapse = ","),
  chain_ppc, burn_ppc, thin_ppc, re_dist_ppc
))
cat(sprintf("REP_ID=%s, initial_seed=%d\n",
            if (is.finite(rep_id)) as.character(rep_id) else "NA",
            initial_seed))
cat("OUT_DIR =", OUT_DIR, "\n")

for (sc in scenarios) {
  for (li in seq_along(zi_intercepts)) {
    cat(sprintf("\n==== Scenario %d | 零膨胀截距 = %.2f (%d/%d) ====\n",
                sc, zi_intercepts[li], li, length(zi_intercepts)))

    dat <- generate_data(
      n = n_ppc,
      seed = initial_seed + 10L * as.integer(sc) + li,
      scenario = sc,
      C = C_ord,
      re_dist = re_dist_ppc,
      random_nis = random_nis,
      nis_range = nis_range,
      nis = nis_fixed,
      zi_intercept = zi_intercepts[li]
    )
    obs_zero <- mean(dat$y1 == 0)
    cat(sprintf("观测零比例 = %.1f%%, N = %d\n", obs_zero * 100, dat$N))

    fit <- fit_joint_model(dat, chain = chain_ppc, burn = burn_ppc, thin = thin_ppc)
    reps <- generate_joint_ppc_reps(fit, dat)
    Y1rep <- reps$Y1rep
    Y2rep <- reps$Y2rep
    assoc <- compute_association_ppc(dat$y1, dat$y2, Y1rep, Y2rep, dat$C)

    panel_lab <- sprintf("Scenario%d \u00B7 Zero rate %.1f%%", sc, obs_zero * 100)

    obs_prop <- sapply(ks, function(k) mean(dat$y1 == k))
    rep_prop <- sapply(ks, function(k) rowMeans(Y1rep == k))

    dist_df[[idx]] <- data.frame(
      Scenario  = paste0("Scenario ", sc),
      ZiInt     = zi_intercepts[li],
      ZeroRate  = obs_zero,
      k         = ks,
      obs       = obs_prop,
      pred      = colMeans(rep_prop),
      lo        = apply(rep_prop, 2, quantile, 0.025),
      hi        = apply(rep_prop, 2, quantile, 0.975)
    )

    cor_rep_df[[idx]] <- bind_rows(
      data.frame(
        Panel = panel_lab, Statistic = "cor(y1, y2)",
        value = assoc$rep$cor_y1_y2, obs = assoc$obs$cor_y1_y2
      ),
      data.frame(
        Panel = panel_lab, Statistic = "cor(I(y1>0), y2)",
        value = assoc$rep$cor_y1pos_y2, obs = assoc$obs$cor_y1pos_y2
      )
    )

    mean_summ <- t(apply(assoc$rep$mean_y1_by_y2, 2, summarize_rep))
    mean_rep_df[[idx]] <- data.frame(
      Panel    = panel_lab,
      y2_level = seq_len(dat$C),
      obs      = assoc$obs$mean_y1_by_y2,
      pred     = mean_summ[, "pred"],
      lo       = mean_summ[, "lo"],
      hi       = mean_summ[, "hi"]
    )

    s_cor12  <- summarize_rep(assoc$rep$cor_y1_y2)
    s_corpos <- summarize_rep(assoc$rep$cor_y1pos_y2)
    assoc_summary_df[[idx]] <- data.frame(
      Panel = panel_lab,
      Statistic = c(
        "cor(y1,y2)", "cor(I(y1>0),y2)",
        paste0("mean(y1|y2=", seq_len(dat$C), ")")
      ),
      Obs = c(assoc$obs$cor_y1_y2, assoc$obs$cor_y1pos_y2, assoc$obs$mean_y1_by_y2),
      Pred = c(s_cor12["pred"], s_corpos["pred"], mean_summ[, "pred"]),
      Lo = c(s_cor12["lo"], s_corpos["lo"], mean_summ[, "lo"]),
      Hi = c(s_cor12["hi"], s_corpos["hi"], mean_summ[, "hi"])
    )

    idx <- idx + 1L
  }
}

dist_df <- bind_rows(dist_df)

# 面板标签: 每个面板顶部仅写该面板的实际零膨胀率
dist_df$Scenario <- factor(dist_df$Scenario, levels = paste0("Scenario ", scenarios))
dist_df$ScenShort <- factor(paste0("Scenario", as.integer(sub("Scenario ", "", dist_df$Scenario))))
dist_df$ZeroLab  <- sprintf("Zero rate %.1f%%", dist_df$ZeroRate * 100)
dist_df$Panel    <- paste0(dist_df$ScenShort, " \u00B7 ", dist_df$ZeroLab)

panel_levels <- dist_df %>%
  distinct(Scenario, ZiInt, Panel) %>%
  arrange(Scenario, desc(ZiInt)) %>%
  pull(Panel)
dist_df$Panel <- factor(dist_df$Panel, levels = panel_levels)

# =============================================================================
# 3. 作图: Fig.4/5/6 合并单页 (a)(b)(c), 并可选输出单张图
# =============================================================================
col_obs  <- "#1B9E91"
col_pred <- "#D7301F"
n_zi <- length(zi_intercepts)

p_dist <- ggplot(dist_df, aes(x = k)) +
  geom_col(aes(y = obs), fill = col_obs, width = 0.7, alpha = 0.85) +
  geom_errorbar(aes(ymin = lo, ymax = hi), color = col_pred, width = 0.3, linewidth = 0.4) +
  geom_point(aes(y = pred), color = col_pred, size = 1.3) +
  facet_wrap(~ Panel, ncol = n_zi) +
  labs(x = "Count value", y = "Proportion") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )

cor_rep_df <- bind_rows(cor_rep_df)
mean_rep_df <- bind_rows(mean_rep_df)
assoc_summary_df <- bind_rows(assoc_summary_df)

write.csv(
  assoc_summary_df,
  file.path(OUT_DIR, "ppc_association_summary.csv"),
  row.names = FALSE
)

build_cor_ppc_plot <- function(stat_label, show_caption = FALSE) {
  sub_df <- cor_rep_df[cor_rep_df$Statistic == stat_label, ]
  sub_df$Panel <- factor(sub_df$Panel, levels = panel_levels)

  ggplot(sub_df, aes(x = value)) +
    geom_histogram(bins = 30, fill = col_pred, color = "white", alpha = 0.75) +
    geom_vline(aes(xintercept = obs), color = col_obs, linewidth = 0.9) +
    facet_wrap(~ Panel, ncol = n_zi, scales = "free") +
    labs(
      x = stat_label,
      y = "Density",
      caption = if (show_caption) {
        "Blue line: observed; Red histogram: posterior predictive"
      } else {
        NULL
      }
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold")
    )
}

p_cor_y1_y2 <- build_cor_ppc_plot("cor(y1, y2)")
p_cor_y1pos_y2 <- build_cor_ppc_plot("cor(I(y1>0), y2)")

p_combined <- (p_dist / p_cor_y1_y2 / p_cor_y1pos_y2) +
  plot_layout(guides = "collect", heights = c(1, 1, 1)) +
  plot_annotation(
    tag_levels = "a",
    tag_prefix = "(",
    tag_suffix = ")",
    theme = theme(
      plot.tag = element_text(face = "bold", size = 14),
      plot.caption = element_text(size = 9, hjust = 0)
    )
  )

combined_width <- 14
combined_height <- 17
ggplot2::ggsave(
  filename = file.path(OUT_DIR, "ppc_fig4_5_6_combined.pdf"),
  plot = p_combined,
  width = combined_width,
  height = combined_height,
  device = grDevices::pdf,
  limitsize = FALSE
)

# 论文用简版: 仅分布 PPC + cor(I(y1>0), y2), 标注 (a)(b)
p_combined_ab <- (p_dist / p_cor_y1pos_y2) +
  plot_layout(guides = "collect", heights = c(1, 1)) +
  plot_annotation(
    tag_levels = "a",
    tag_prefix = "(",
    tag_suffix = ")",
    caption = "Teal: observed; Red: posterior predictive.",
    theme = theme(
      plot.tag = element_text(face = "bold", size = 14),
      plot.caption = element_text(size = 9, hjust = 0)
    )
  )

combined_ab_height <- 12
ggplot2::ggsave(
  filename = file.path(OUT_DIR, "ppc_fig4_6_combined.pdf"),
  plot = p_combined_ab,
  width = combined_width,
  height = combined_ab_height,
  device = grDevices::pdf,
  limitsize = FALSE
)

# 可选: 仍输出单张图
ggplot2::ggsave(
  file.path(OUT_DIR, "ppc_fig4_distribution.pdf"),
  p_dist, width = 14, height = 6, device = grDevices::pdf
)
ggplot2::ggsave(
  file.path(OUT_DIR, "ppc_fig5_cor_y1_y2.pdf"),
  p_cor_y1_y2, width = 14, height = 6, device = grDevices::pdf
)
ggplot2::ggsave(
  file.path(OUT_DIR, "ppc_fig6_cor_y1pos_y2.pdf"),
  p_cor_y1pos_y2, width = 14, height = 6, device = grDevices::pdf
)

# =============================================================================
# 4. 按 y2 分组的 y1 均值 PPC
# =============================================================================
mean_rep_df$y2_label <- factor(
  mean_rep_df$y2_level,
  levels = sort(unique(mean_rep_df$y2_level)),
  labels = paste0("y2 = ", sort(unique(mean_rep_df$y2_level)))
)
mean_rep_df$Panel <- factor(mean_rep_df$Panel, levels = panel_levels)

p_mean_by_y2 <- ggplot(mean_rep_df, aes(x = y2_label)) +
  geom_col(aes(y = obs), fill = col_obs, width = 0.55, alpha = 0.85) +
  geom_errorbar(aes(ymin = lo, ymax = hi), color = col_pred, width = 0.15, linewidth = 0.5) +
  geom_point(aes(y = pred), color = col_pred, size = 2.2) +
  facet_wrap(~ Panel, ncol = n_zi, scales = "free_y") +
  labs(
    x = NULL,
    y = expression("Mean of " ~ y[1]),
    caption = "Bars: observed; Red point/interval: posterior predictive mean and 95% interval"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

ggplot2::ggsave(
  file.path(OUT_DIR, "ppc_fig7_mean_y1_by_y2.pdf"),
  p_mean_by_y2, width = 14, height = 6, device = grDevices::pdf
)

# 辅助筛选指标（仅供排序参考，最终以目视为准）
dist_cover <- mean(dist_df$obs >= dist_df$lo & dist_df$obs <= dist_df$hi, na.rm = TRUE)
mean_cover <- mean(mean_rep_df$obs >= mean_rep_df$lo & mean_rep_df$obs <= mean_rep_df$hi, na.rm = TRUE)
cor_in_interval <- assoc_summary_df %>%
  dplyr::filter(Statistic %in% c("cor(y1,y2)", "cor(I(y1>0),y2)")) %>%
  dplyr::summarise(rate = mean(Obs >= Lo & Obs <= Hi, na.rm = TRUE)) %>%
  dplyr::pull(rate)
meta <- data.frame(
  rep_id = if (is.finite(rep_id)) rep_id else NA_integer_,
  initial_seed = initial_seed,
  scenario = paste(scenarios, collapse = ","),
  re_dist = re_dist_ppc,
  n = n_ppc,
  chain = chain_ppc,
  burn = burn_ppc,
  thin = thin_ppc,
  dist_obs_in_95 = dist_cover,
  mean_obs_in_95 = mean_cover,
  cor_obs_in_95 = cor_in_interval,
  score = 0.5 * dist_cover + 0.25 * mean_cover + 0.25 * cor_in_interval,
  finished_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
)
write.csv(meta, file.path(OUT_DIR, "ppc_candidate_meta.csv"), row.names = FALSE)

cat("\n完成: 已输出以下文件到", OUT_DIR, "\n")
cat("  - ppc_fig4_6_combined.pdf  (Fig.4+6 合并, 标注 (a)(b), 论文简版)\n")
cat("  - ppc_fig4_5_6_combined.pdf (Fig.4/5/6 合并, 标注 (a)(b)(c))\n")
cat("  - ppc_fig4_distribution / ppc_fig5_cor_y1_y2 / ppc_fig6_cor_y1pos_y2\n")
cat("  - ppc_fig7_mean_y1_by_y2.pdf (按 y2 的 y1 均值 PPC)\n")
cat("  - ppc_association_summary.csv / ppc_candidate_meta.csv\n")
cat(sprintf("  score=%.3f (dist_cover=%.3f, mean_cover=%.3f, cor_cover=%.3f)\n",
            meta$score, meta$dist_obs_in_95, meta$mean_obs_in_95, meta$cor_obs_in_95))
