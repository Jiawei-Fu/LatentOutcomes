# LatentOutcomes
R package for estimating causal effects with latent outcomes.

Authors: Jiawei Fu; Donald P. Green

Contributors: Jiawei Fu; Donald P. Green

The package estimates causal effects with latent outcomes.

To install and use the latest version of the package, run:
```r
install.packages("devtools")
devtools::install_github("Jiawei-Fu/LatentOutcomes")
library(LatentOutcomes)
```
Some users may receive time-zone warnings; this does not affect installation.
If you have installed an older version, run:
```r
remove.packages("LatentOutcomes")
devtools::install_github("Jiawei-Fu/LatentOutcomes")
.rs.restartR()
```

Main functions:

- `estlatent(mod, Y, data = NULL, method = "sem", IV_Y = NULL, opt = TRUE)`
- `estlatent_ave(mod, Y, data = NULL, method = "sem", IV_Y = NULL, opt = TRUE)`
- `estlatent_robust(mod, Y, data = NULL, method = "sem", IV_Y = NULL, opt = TRUE)`

Estimation variants:

- `estlatent`: default per-outcome GMM implementation; recommended for primary analysis.
- `estlatent_ave`: average-proxy GMM implementation. It first combines multiple measurements into an average proxy; results are often similar to per-outcome GMM.
- `estlatent_robust`: two-step robust version (first-stage loading estimation, then regression; Godambe sandwich variance). It is typically less efficient but more robust to misspecification in joint GMM.

Example:
```r
data("test_dat", package = "LatentOutcomes")  # input data

# SEM
fit_sem <- estlatent(
  mod = "~Z",
  Y = c("Y1", "Y2", "Y3"),
  data = test_dat,
  method = "sem"
)
summary(fit_sem)

# GMM (efficient two-step by default)
fit_gmm <- estlatent(
  mod = "~Z",
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
  Y = c("Y1", "Y2", "Y3"),
  data = test_dat,
  method = "sem"
)
summary(fit_sem_cov)

# GMM with covariates
fit_gmm_cov <- estlatent(
  mod = "~Z+x1+x2",
  Y = c("Y1", "Y2", "Y3"),
  data = test_dat,
  method = "gmm",
  IV_Y = c("Z", "Y1", "Y2", "Y3"),
  opt = TRUE
)
summary(fit_gmm_cov)

# average-proxy GMM
fit_gmm_ave <- estlatent_ave(
  mod = "~Z+x1+x2",
  Y = c("Y1", "Y2", "Y3"),
  data = test_dat,
  method = "gmm",
  IV_Y = c("Z"),
  opt = TRUE
)
summary(fit_gmm_ave)

# Two-step robust GMM
fit_gmm_robust <- estlatent_robust(
  mod = "~Z+x1+x2",
  Y = c("Y1", "Y2", "Y3"),
  data = test_dat,
  method = "gmm",
  IV_Y = c("Z"),
  opt = TRUE
)
summary(fit_gmm_robust)
```

Notes:

- `mod` controls the latent regression specification `eta ~ ...`.
- Recommended `mod` format is RHS-only, e.g. `~Z+x1+x2`.
- `mod` is required.
- `IV_Y` is for loading moments in `method = "gmm"` only.
- For `method = "sem"`, `IV_Y` is not used; SEM is estimated from the full model specification directly.
- Regression moments in GMM automatically use regressors from `mod` as instruments.
- `IV_Y` may include variables from `data` even if they are not in `mod`.
- Use `summary(fit)` to print coefficient tables.

Model misspecification:

In the following simulation, the data-generating process includes an additional latent component `eta_alt` that enters the second and third measurements:

```r
n <- 5000
X <- 3 * rnorm(n)
ATE <- 2

eta_0 <- X + rnorm(n)
eta_1 <- ATE + X + rnorm(n)

eta_alt_0 <- X + rnorm(n)  # additional latent outcome
eta_alt_1 <- X + rnorm(n)

y1_0 <- 1.0 * eta_0 + 1 * rnorm(n)
y1_1 <- 1.0 * eta_1 + 1 * rnorm(n)

y2_0 <- 0.5 * eta_0 + 2 * rnorm(n) + 1 * eta_alt_0  # additional latent outcome affects Y2
y2_1 <- 0.5 * eta_1 + 2 * rnorm(n) + 1 * eta_alt_1

y3_0 <- 2 * eta_0 + 2 * rnorm(n) + 2 * eta_alt_0  # additional latent outcome affects Y3
y3_1 <- 2 * eta_1 + 2 * rnorm(n) + 2 * eta_alt_1

# treatment assignment
trt <- complete_ra(n)

Y1 <- y1_1 * trt + (1 - trt) * y1_0
Y2 <- y2_1 * trt + (1 - trt) * y2_0
Y3 <- y3_1 * trt + (1 - trt) * y3_0

dat <- data.frame(trt, Y1, Y2, Y3, X)
```

In this setting, adding `X` to `mod` (`~trt+X`) can make single-step GMM (jointly estimating loading and regression moments) unstable for loadings, because `X`-related regression moments can pull loading estimates.

`estlatent_robust` is designed for this case:

- stage 1 identifies loadings using loading-IV moments (`IV_Y`, for example `trt`);
- stage 2 estimates regression coefficients conditional on stage-1 loadings;
- variance is adjusted with stacked estimating equations (Godambe sandwich).

So when you move from `mod = "~trt"` to `mod = "~trt+X"` under this type of DGP, `estlatent_robust(..., method = "gmm")` is recommended. In practice, researchers can use `estlatent` first. If results from `estlatent_robust` are very different from `estlatent`, and all IVs are valid, this suggests that covariate-related regression moments may be influencing loading estimation in joint GMM. In that case, `estlatent_robust` is recommended.

If you have comments, suggestions, or findings, please email: jiawei.fu@duke.edu
