# zinbcdpmm

Bayesian joint model for longitudinal **zero-inflated negative binomial (ZINB)** and **ordinal** outcomes with **CDPMM** random effects (Liang & Du).

贝叶斯联合模型：纵向零膨胀负二项 + 有序多分类，随机效应为中心化 Dirichlet 过程混合模型（CDPMM）。

## Install / 安装

```r
# install.packages("remotes")
remotes::install_github("Liang-ShuQing/zinb-ordinal-cdpmm-code", subdir = "zinbcdpmm")
```

## Tiny example / 小例子

```r
library(zinbcdpmm)
set.seed(1)
dat <- simulate_zinb_ordinal(n = 30, nis = 4, scenario = 2, re_dist = "normal")
fit <- fit_zinb_ordinal(dat, chain = 200, burn = 100, thin = 5, G = 4)
print(fit)
summary(fit)
```

## Paper code / 论文代码仓库

Simulation and analysis scripts (this repository):
https://github.com/Liang-ShuQing/zinb-ordinal-cdpmm-code

## License

MIT © 2026 Shuqing Liang and Jiang Du
