# Shared model cores

| File | Role |
|------|------|
| `prog1_cdpmm_joint_vs_separate.R` | Main CDPMM joint Gibbs sampler + data generation (simulation) |
| `prog2_cdpmm_vs_gaussian_joint.R` | CDPMM vs homogeneous Gaussian joint comparator |
| `joint_cdpmm.stan` | Stan model for Supporting Material computational benchmark |
| `joint_cdpmm.jags` | JAGS model for Supporting Material computational benchmark |
| `build_sim_tables.R` | Main-text Monte Carlo table builder (scrapes `joint_param_summary_*.csv`, `variance_ratio_*.csv`, separate-model summaries; prints LaTeX) |

Many figure/table folders already embed a local copy of `prog1_…` so they run from that directory. Prefer the local copy when present; otherwise copy from here.

`build_sim_tables.R` is the authoritative copy: the Scenario 1/2 tables ([`../main_paper/tab_sim_scen1_normal`](../main_paper/tab_sim_scen1_normal/), `tab_sim_scen2_normal`, `tab_sim_scen1`, `tab_sim_main`) and the Scenario 3 / MVT builders (`build_sim_tables_scen3.R`, `build_sim_tables_mvt8.R`) all reuse its helpers.
