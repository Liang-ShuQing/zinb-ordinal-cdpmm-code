# `fig:hrs_assoc`

Joint association between hospital nights (HSPNIT) and self-rated health (SHLT) in the HRS application:
(a) \(P(\text{HSPNIT}\mid\text{SHLT})\), (b) \(P(\text{SHLT}\mid\text{HSPNIT})\).

## Files

- `plot_joint_assoc_hspnit_shlt_fig1a.R` — builds the one-panel (`joint_hspnit_shlt_fig1a.pdf`) and the two-panel (`joint_hspnit_shlt_fig1ab.pdf`) versions used in the paper

## Data (not shipped)

The script needs the analysis-ready HRS extract (`Data_model_complete.csv`), which is **not** distributed here. Obtain the HRS files under their terms (<https://hrs.isr.umich.edu/>) and pass the CSV path as the first argument.

## Run

```bash
Rscript plot_joint_assoc_hspnit_shlt_fig1a.R /path/to/Data_model_complete.csv ./results
```

Outputs: `results/joint_hspnit_shlt_fig1ab.pdf` (paper figure), `results/joint_hspnit_shlt_fig1a.pdf`, plus copies under `results/paper_figures/`.

Paper figure: `joint_hspnit_shlt_fig1ab.pdf`.
