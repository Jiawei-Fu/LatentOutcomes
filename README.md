# LatentOutcomes
R package for estimating causal effects with latent outcomes.

Author: Jiawei Fu; Donald P. Green

Contributor: Jiawei Fu; Donald P. Green

The package estimates causal effects with latent outcomes.

To install and use the latest version of the package, try the following codes:
```r
install.packages("devtools")
devtools::install_github("Jiawei-Fu/LatentOutcomes")
library(LatentOutcomes)
```
Some of users will receive warnings about time zone, it will not affect the installation.

Main function:

- `estlatent(mod = NULL, Z, Y, data = NULL, method = "sem", IV_Y = NULL, opt = TRUE)`

Example:
```r
data("test_dat", package = "LatentOutcomes")  # input data

# SEM
fit_sem <- estlatent(
  mod = "~Z",
  Z = "Z",
  Y = c("Y1", "Y2", "Y3"),
  data = test_dat,
  method = "sem"
)
summary(fit_sem)

# GMM (efficient two-step by default)
fit_gmm <- estlatent(
  mod = "~Z",
  Z = "Z",
  Y = c("Y1", "Y2", "Y3"),
  data = test_dat,
  method = "gmm",
  IV_Y = c("Z", "Y1", "Y2", "Y3"),
  opt = TRUE
)
summary(fit_gmm)

# SEM with covariates in mod
fit_sem_cov <- estlatent(
  mod = "~Z+x1+x2",
  Z = "Z",
  Y = c("Y1", "Y2", "Y3"),
  data = test_dat,
  method = "sem"
)
summary(fit_sem_cov)

# GMM with covariates
fit_gmm_cov <- estlatent(
  mod = "~Z+x1+x2",
  Z = "Z",
  Y = c("Y1", "Y2", "Y3"),
  data = test_dat,
  method = "gmm",
  IV_Y = c("Z", "Y1", "Y2", "Y3"),
  opt = TRUE
)
summary(fit_gmm_cov)
```

Notes:

- `mod` controls the latent regression specification `eta ~ ...`.
- Recommended `mod` format is RHS-only, e.g. `~Z+x1+x2` or `~Z+x1+x2:x3`.
- `mod = NULL` means treatment-only model (`~Z`).
- `:` denotes interaction terms (for example, `x2:x3`).
- `IV_Y` is for loading moments in `method = "gmm"` only.
- For `method = "sem"`, `IV_Y` is not used; SEM is estimated from the full model specification directly.
- Regression moments in GMM automatically use covariates from `mod` as IVs.
- Use `summary(fit)` to print coefficient tables.

If you have any comments, suggestions, or findings, appreciate sending me the email: jiawei.fu@duke.edu
