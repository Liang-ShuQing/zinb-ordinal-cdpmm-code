# `tab:zip_misspec`

Robustness check with a ZIP count-part data-generating process (Scenario 2, Truth 2, \(N=200\), \(S=500\)).
Only the identified parameters are reported: \(\beta_1\) and \(r\) are dropped because the generated counts carry no overdispersion.

## Code

- `build_zip_misspec_table.R` — parses the joint / ordinal-only / ZINB-only summary blocks and writes the LaTeX table
- Driver: [`../../00_shared/prog1_cdpmm_joint_vs_separate.R`](../../00_shared/prog1_cdpmm_joint_vs_separate.R)

## Inputs (not shipped)

One run directory containing:

- `slurm-*.out` (joint + separate summary blocks, with average interval length)
- `joint_param_summary_*.csv`, `variance_ratio_*.csv`, `criteria_winrate_*.csv`

## Run

```bash
Rscript build_zip_misspec_table.R /path/to/zip_run ./results
```

Output: `results/zip_misspec_table.tex` (tabular body only; the float wrapper lives in the manuscript).
