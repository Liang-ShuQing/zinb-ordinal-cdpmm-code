# `tab:sm_stress` (Tables S24–S25)

Supporting Material design-stress tables:

| Table | Content |
|-------|---------|
| S24 | small-\(N\) cells |
| S25 | high-dimensional covariate cells (extra noise covariates) |

## Code

- `build_sm_stress.R` — reads the downloaded HPC summary CSVs and writes the LaTeX fragments

## Outputs

Committed next to this folder (bodies are inlined into the Supporting Material, no `\input`):
[`../tab_sm_stress_smallN_hd_results`](../tab_sm_stress_smallN_hd_results/) → `table_S24.tex`, `table_S25.tex`.

## Run

```bash
Rscript build_sm_stress.R /path/to/stress_results ./results
```

The result folder must contain the per-cell summary CSVs of the stress runs; those Monte Carlo
archives are **not** shipped. Output goes to a sibling `results/` folder (override with the second argument).
