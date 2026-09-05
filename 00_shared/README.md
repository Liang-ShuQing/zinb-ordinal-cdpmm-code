# Shared model cores

| File | Role |
|------|------|
| `prog1_cdpmm_joint_vs_separate.R` | Main CDPMM joint Gibbs sampler + data generation (simulation) |
| `prog2_cdpmm_vs_gaussian_joint.R` | CDPMM vs homogeneous Gaussian joint comparator |
| `joint_cdpmm.stan` | Stan model for Supporting Material computational benchmark |
| `joint_cdpmm.jags` | JAGS model for Supporting Material computational benchmark |

Many figure/table folders already embed a local copy of `prog1_…` so they run from that directory. Prefer the local copy when present; otherwise copy from here.
