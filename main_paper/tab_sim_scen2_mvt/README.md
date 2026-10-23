# `tab:sim_scen2_mvt`

Monte Carlo table: Scenario 2, **multivariate-\(t\) (\(df=8\))** random-effects truth, \(N\in\{200,400\}\), \(S=500\).

## Code

- `build_sim_tables_mvt8.R` — one run writes `tab:sim_scen1_mvt` and `tab:sim_scen2_mvt`
- `build_sim_tables.R` — shared scraping / formatting helpers
- Driver: [`../../00_shared/prog1_cdpmm_joint_vs_separate.R`](../../00_shared/prog1_cdpmm_joint_vs_separate.R) with `SCENARIO=2`, `RE_DIST=mvt`, `MVT_DF=8`, `G=8`

## Settings

| Item | Value |
|------|--------|
| Scenario | 2 |
| RE_DIST | mvt (\(df=8\)) |
| \(N\) | 200 / 400 |
| \(S\) | 500 |
| Chains | 10000 iter, 5000 burn-in, thin 5, \(G=8\) |

## Run

```bash
Rscript build_sim_tables_mvt8.R /path/to/hpc_results ./sim_tables_mvt8.tex
```

The results folder must contain `s1_n200_mvt8/`, `s1_n400_mvt8/`, `s2_n200_mvt8/`, `s2_n400_mvt8/`.
Result archives are **not** shipped.

See also: [`../tab_sim_scen1_mvt`](../tab_sim_scen1_mvt/), [`../tab_sim_scen2_normal`](../tab_sim_scen2_normal/).
