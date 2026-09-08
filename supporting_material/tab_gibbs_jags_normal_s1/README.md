# `tab:gibbs_jags_normal_s1`

Supporting Material: Bias / RMSE / coverage under **Gibbs | JAGS | Stan** (normal RE; Scenario 1; \(N\in\{200,400\}\)).

## Files (canonical scheme-1 templates for all four Gibbs–JAGS–Stan tables)

| File | Role |
|------|------|
| `prog1_cdpmm_joint_vs_separate.R` | Gibbs blocked sampler (main Monte Carlo column) |
| `code_stan.R` | Stan Monte Carlo driver |
| `code_jags.R` | JAGS Monte Carlo driver |
| `joint_cdpmm.stan` / `joint_cdpmm.jags` | Model files |
| `build_gibbs_jags_tables.R` | Build SM LaTeX tables from CSV summaries |
| `R_stan_example.slurm` | Example Slurm (adapt partition / cores) |

Sibling table folders only include the builder and point here.

## Scheme 1 (retained draws ≈ 1000)

- Gibbs: chain \(10{,}000\), burn \(5{,}000\), thin \(5\)
- JAGS: 2 chains × adapt/burn/sampling \(500\), thin \(1\)
- Stan: 2 chains × warmup/sampling \(500\); prefer `STAN_PARALLEL_CHAINS=1` when forking many replicates

## Run (outline)

```bash
# Per cell: set SCENARIO=1 RE_DIST=normal N=200|400 NSIM=500
Rscript code_stan.R .
Rscript code_jags.R .
# Gibbs column usually taken from the main simulation (prog1) outputs
Rscript build_gibbs_jags_tables.R
```
