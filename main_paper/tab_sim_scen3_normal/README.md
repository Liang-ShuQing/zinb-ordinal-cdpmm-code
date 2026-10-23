# `tab:sim_scen3_normal`

Monte Carlo table: Scenario 3, **multivariate normal** random-effects truth, \(N\in\{200,400\}\), \(S=500\).

## Code

- `build_sim_tables_scen3.R` — one run writes all three Scenario 3 tables (`tab:sim_scen3_normal`, `tab:sim_scen3`, `tab:sim_scen3_mvt`); this folder documents the normal-truth row
- `build_sim_tables.R` — shared scraping / formatting helpers (also in [`../../00_shared`](../../00_shared/))
- Driver: [`../../00_shared/prog1_cdpmm_joint_vs_separate.R`](../../00_shared/prog1_cdpmm_joint_vs_separate.R) with `SCENARIO=3`, `RE_DIST=normal`, `G=8`

## Settings

| Item | Value |
|------|--------|
| Scenario | 3 |
| RE_DIST | normal |
| \(N\) | 200 / 400 |
| \(S\) | 500 |
| Chains | 10000 iter, 5000 burn-in, thin 5, \(G=8\) |

## Run

```bash
Rscript build_sim_tables_scen3.R /path/to/hpc_results ./sim_tables_scen3.tex
```

The results folder must contain `s3_n200_normal/`, `s3_n400_normal/`, `s3_n200_mixture/`,
`s3_n400_mixture/`, `s3_n200_mvt8/`, `s3_n400_mvt8/`, each with `joint_param_summary_*.csv`,
`variance_ratio_*.csv` and the separate-model summaries. Those Monte Carlo archives are **not** shipped.

See also: [`../tab_sim_scen3`](../tab_sim_scen3/), [`../tab_sim_scen3_mvt`](../tab_sim_scen3_mvt/).
