# `fig:epsr_sim`

EPSR (R-hat) trajectories for \(\boldsymbol{\alpha}\), \(\boldsymbol{\beta}\), \(\boldsymbol{\gamma}\), and \(r\) (Scenario 2, mixture truth, \(N=200\); chain length \(10{,}000\)).

## Files

- `plot_epsr_trace.R` — entry (produces both trace and EPSR PDFs)
- `prog1_cdpmm_joint_vs_separate.R` — sourced helpers

Same entry script as `fig_trace_sim`.

## Run

```bash
Rscript plot_epsr_trace.R ./output
```

Expected outputs include `GelmanRubinConvergencePlot10000.pdf` (name depends on `CHAIN`).
