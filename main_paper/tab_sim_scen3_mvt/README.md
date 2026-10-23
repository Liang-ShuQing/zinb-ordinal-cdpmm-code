# `tab:sim_scen3_mvt`

Monte Carlo table: Scenario 3, **multivariate-\(t\) (\(df=8\))** random-effects truth, \(N\in\{200,400\}\), \(S=500\).

## Code

- `build_sim_tables_scen3.R` — one run writes all three Scenario 3 tables (`tab:sim_scen3_normal`, `tab:sim_scen3`, `tab:sim_scen3_mvt`); this folder documents the heavy-tail row
- `build_sim_tables.R` — shared scraping / formatting helpers
- Driver: [`../../00_shared/prog1_cdpmm_joint_vs_separate.R`](../../00_shared/prog1_cdpmm_joint_vs_separate.R) with `SCENARIO=3`, `RE_DIST=mvt`, `MVT_DF=8`, `G=8`

## Settings

| Item | Value |
|------|--------|
| Scenario | 3 |
| RE_DIST | mvt (\(df=8\)) |
| \(N\) | 200 / 400 |
| \(S\) | 500 |
| Chains | 10000 iter, 5000 burn-in, thin 5, \(G=8\) |

## Run

```bash
Rscript build_sim_tables_scen3.R /path/to/hpc_results ./sim_tables_scen3.tex
```

The `mvt8` cells are read from `s3_n200_mvt8/` and `s3_n400_mvt8/`; the other four cell folders are
needed only because the same script writes the normal and mixture tables. Result archives are **not** shipped.

See also: [`../tab_sim_scen3_normal`](../tab_sim_scen3_normal/), [`../tab_sim_scen3`](../tab_sim_scen3/).
