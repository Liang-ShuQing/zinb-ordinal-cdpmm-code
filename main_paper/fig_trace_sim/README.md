# `fig:trace_sim`

Three-chain MCMC trace plots for \(\boldsymbol{\alpha}\), \(\boldsymbol{\beta}\), \(\boldsymbol{\gamma}\), and \(r\) (Scenario 2, mixture truth, \(N=200\); chain length \(10{,}000\)).

## Files

- `plot_epsr_trace.R` — entry (produces both trace and EPSR PDFs)
- `prog1_cdpmm_joint_vs_separate.R` — sourced helpers (`SKIP_MAIN_SIM`)

Same entry script as `fig_epsr_sim`.

## Run

```bash
# From this folder; optional OUT_DIR as first argument
Rscript plot_epsr_trace.R ./output
```

Environment overrides (optional): `CHAIN=10000`, `SEED_DATA=2025`.

Expected outputs include `three_chain10000.pdf` (filename may vary with `CHAIN`).
