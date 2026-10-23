# `tab:cdpmm_vs_gauss_mvt_s3`

Supporting Material table: CDPMM joint prior versus a **homogeneous-Gaussian** joint prior under
**Scenario 3**, multivariate-\(t\) (**\(df=8\)**) random-effects truth, \(N\in\{200,400\}\), \(S=500\).

## Code

- `build_sm_cdpmm_vs_gauss_rest.py` — one run writes five tables: `tab:cdpmm_vs_gauss_mvt_s1`,
  `tab:cdpmm_vs_gauss_mvt_s2`, `tab:cdpmm_vs_gauss_normal_s3`, `tab:cdpmm_vs_gauss_mixture_s3`,
  `tab:cdpmm_vs_gauss_mvt_s3`
- CDPMM columns are read from the main-paper joint summaries (see [`../../main_paper/tab_sim_scen3_mvt`](../../main_paper/tab_sim_scen3_mvt/));
  the Gaussian columns come from the Gaussian-only \(S=500\) runs
- Column layout and bolding match Tables S1–S4 in [`../tab_cdpmm_vs_gauss_normal_s1`](../tab_cdpmm_vs_gauss_normal_s1/)

## Run

```bash
python build_sm_cdpmm_vs_gauss_rest.py ./sim_cdpmm_vs_gauss_tables_rest.tex
```

Expected project layout (result archives are **not** shipped):

```text
hpc_jobs/scen3v2_mvt8_n200n400_s500_wuzhen_20260915/s3_n{200,400}_{normal,mixture,mvt8}/
hpc_jobs/mvt8_sc12_n200n400_s500_taiyuan_20260915/s{1,2}_n{200,400}_mvt8/
hpc_jobs/sm_gauss_sc12_mvt8_s500_taiyuan_20260918/     # Gaussian-only, Sc1/2
hpc_jobs/sm_gauss_sc3_s500_wuzhen_20260918/            # Gaussian-only, Sc3
```

Bolding: Bias is bold when \(|\text{Bias}|\ge 0.10\), CP is bold when \(\text{CP}<0.90\) (RMSE is never bolded).
