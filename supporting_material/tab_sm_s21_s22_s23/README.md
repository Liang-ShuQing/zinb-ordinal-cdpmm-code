# `tab:weak_rho`, `tab:zi_severity`, `tab:rslope`

Supporting Material Tables S21–S23:

| Table | Content |
|-------|---------|
| S21 | random-slope (`tab:rslope`) extension under the multivariate-\(t\) truth |
| S22 | zero-inflation severity (`tab:zi_severity`) |
| S23 | weak random-effect correlation (`tab:weak_rho`) |

## Code

- `build_s21_s22_s23.R` — reads the downloaded HPC summary CSVs and writes the three LaTeX fragments

## Outputs

The generated fragments are committed next to this folder, because the project inlines table bodies
(no `\input`): [`../tab_sm_s21_s22_s23_results`](../tab_sm_s21_s22_s23_results/) →
`table_S21.tex`, `table_S22.tex`, `table_S23.tex`.

## Run

```bash
Rscript build_s21_s22_s23.R /path/to/hpc_jobs ./results
```

`/path/to/hpc_jobs` must hold the corresponding `s*_n*/joint_param_summary_*.csv` and
`variance_ratio_*.csv` files; \(S=500\) Monte Carlo archives are **not** shipped.
Output goes to a sibling `results/` folder (override with the second argument).
