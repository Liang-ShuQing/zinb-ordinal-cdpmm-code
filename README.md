# Code by paper figure / table

**Repository:** https://github.com/Liang-ShuQing/zinb-ordinal-cdpmm-code

Companion code for the manuscript **Bayesian Semiparametric Joint Modeling of Longitudinal Zero-Inflated Count and Ordinal Outcomes** (ZINB–cumulative logit joint model with CDPMM random effects).

This package is organized **one folder per figure or table**.

## R package: zinbcdpmm

Reusable R package for the joint ZINB–ordinal model with CDPMM random effects (Gibbs sampler). Public API: `simulate_zinb_ordinal()`, `fit_zinb_ordinal()`, plus `print` / `summary` methods.

**Install**

```r
# install.packages("remotes")
remotes::install_github("Liang-ShuQing/zinb-ordinal-cdpmm-code", subdir = "zinbcdpmm")
```

**Minimal usage** (short chains for a quick demo)

```r
library(zinbcdpmm)
set.seed(1)
dat <- simulate_zinb_ordinal(n = 30, nis = 4, scenario = 2, re_dist = "normal")
fit <- fit_zinb_ordinal(dat, chain = 200, burn = 100, thin = 5, G = 4)
print(fit)
summary(fit)
```

The figure/table folders below remain for paper replication; `zinbcdpmm/` is the reusable fitting interface. Separate ZINB/ordinal fitters exist internally but are not exported yet. HRS microdata are not included.

## What is included

- Runnable R / Stan / JAGS / LaTeX **source scripts** only
- Shared model cores under `00_shared/`

## What is **not** included

- Monte Carlo result RDS/CSV archives (large)
- HRS microdata (`Data_model_complete.csv` etc.) — obtain from [HRS](https://hrs.isr.umich.edu/) under their terms
- Compiled PDFs of figures (regenerate with the scripts)

Full replication of \(S=500\) tables typically requires an HPC cluster.

## Environment (typical)

- R ≥ 4.2
- Packages: `coda`, `ggplot2`, `foreach`, `doParallel`, `mvtnorm`, `MASS`, `truncnorm`, `Matrix`, …
- Supporting Material Stan column: `cmdstanr` + CmdStan
- Supporting Material JAGS column: `rjags`

On Linux HPC prefer **fork** `registerDoParallel(cores = …)` (avoid PSOCK at high core counts).

## Main paper index

| Paper label | Folder |
|-------------|--------|
| `fig:graphical_model` | [main_paper/fig_graphical_model](main_paper/fig_graphical_model/) |
| `fig:trace_sim` | [main_paper/fig_trace_sim](main_paper/fig_trace_sim/) |
| `fig:epsr_sim` | [main_paper/fig_epsr_sim](main_paper/fig_epsr_sim/) |
| `fig:ppc_combined` | [main_paper/fig_ppc_combined](main_paper/fig_ppc_combined/) |
| `fig:re_density_panel` | [main_paper/fig_re_density_panel](main_paper/fig_re_density_panel/) |
| `fig:winrate_by_rho` | [main_paper/fig_winrate_by_rho](main_paper/fig_winrate_by_rho/) |
| `fig:G_sensitivity` | [main_paper/fig_G_sensitivity](main_paper/fig_G_sensitivity/) |
| `fig:prior_sensitivity` | [main_paper/fig_prior_sensitivity](main_paper/fig_prior_sensitivity/) |
| `tab:sim_scen1_normal` | [main_paper/tab_sim_scen1_normal](main_paper/tab_sim_scen1_normal/) |
| `tab:sim_scen2_normal` | [main_paper/tab_sim_scen2_normal](main_paper/tab_sim_scen2_normal/) |
| `tab:sim_scen1` | [main_paper/tab_sim_scen1](main_paper/tab_sim_scen1/) |
| `tab:sim_main` | [main_paper/tab_sim_main](main_paper/tab_sim_main/) |
| `tab:baseline` | [main_paper/tab_baseline](main_paper/tab_baseline/) |
| `fig:hrs_trace` | [main_paper/fig_hrs_trace](main_paper/fig_hrs_trace/) |
| `tab:hrs_posterior` | [main_paper/tab_hrs_posterior](main_paper/tab_hrs_posterior/) |

## Supporting Material index

| Label | Folder |
|-------|--------|
| `tab:cdpmm_vs_gauss_normal_s1` | [supporting_material/tab_cdpmm_vs_gauss_normal_s1](supporting_material/tab_cdpmm_vs_gauss_normal_s1/) |
| `tab:cdpmm_vs_gauss_normal_s2` | [supporting_material/tab_cdpmm_vs_gauss_normal_s2](supporting_material/tab_cdpmm_vs_gauss_normal_s2/) |
| `tab:cdpmm_vs_gauss_mixture_s1` | [supporting_material/tab_cdpmm_vs_gauss_mixture_s1](supporting_material/tab_cdpmm_vs_gauss_mixture_s1/) |
| `tab:cdpmm_vs_gauss_mixture_s2` | [supporting_material/tab_cdpmm_vs_gauss_mixture_s2](supporting_material/tab_cdpmm_vs_gauss_mixture_s2/) |
| `tab:gibbs_jags_normal_s1` | [supporting_material/tab_gibbs_jags_normal_s1](supporting_material/tab_gibbs_jags_normal_s1/) |
| `tab:gibbs_jags_normal_s2` | [supporting_material/tab_gibbs_jags_normal_s2](supporting_material/tab_gibbs_jags_normal_s2/) |
| `tab:gibbs_jags_mixture_s1` | [supporting_material/tab_gibbs_jags_mixture_s1](supporting_material/tab_gibbs_jags_mixture_s1/) |
| `tab:gibbs_jags_mixture_s2` | [supporting_material/tab_gibbs_jags_mixture_s2](supporting_material/tab_gibbs_jags_mixture_s2/) |

Shared cores for CDPMM vs Gauss / Gibbs–JAGS–Stan sit in the `*_normal_s1` folders and in `00_shared/`; sibling table folders point there via their README.

All eight Gibbs | JAGS | Stan Monte Carlo cells (\(N\in\{200,400\}\), Scenarios 1–2, normal and mixture) are complete; rebuild LaTeX tables with `build_gibbs_jags_tables.R` after regenerating summaries.
