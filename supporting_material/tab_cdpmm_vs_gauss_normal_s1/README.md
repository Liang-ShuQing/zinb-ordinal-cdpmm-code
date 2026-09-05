# `tab:cdpmm_vs_gauss_normal_s1`

Supporting Material: CDPMM vs Gaussian random-effects prior (normal truth; Scenario 1; \(N\in\{200,400\}\)).

## Files (canonical copy for all four CDPMM-vs-Gauss tables)

- `prog2_cdpmm_vs_gaussian_joint.R` — Monte Carlo driver
- `build_cdpmm_vs_gauss_tables.R` — LaTeX tables
- `plot_cdpmm_vs_gauss.R` — optional figures

Sibling folders only ship the builder and point here for the driver.

## Run

```bash
# Set SCENARIO=1 RE_DIST=normal N=200|400 N_SIM=500 …
Rscript prog2_cdpmm_vs_gaussian_joint.R ./output
Rscript build_cdpmm_vs_gauss_tables.R
```
