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

Main function:

- `estlatent(mod = NULL, Z, Y, data = NULL, method = "sem", IV_Y = NULL, opt = TRUE)`
- `estlatent_ave(mod = NULL, Z, Y, data = NULL, method = "sem", IV_Y = NULL, opt = TRUE)`
- `estlatent_robust(mod = NULL, Z, Y, data = NULL, method = "sem", IV_Y = NULL, opt = TRUE)`

Estimation variants:

- `estlatent`: default per-outcome GMM implementation; this is the recommended option for optimal estimation.
- `estlatent_ave`: average-proxy GMM implementation. In GMM, this variant first combines multiple measurements into an average proxy. Results are often similar to the per-outcome GMM implementation.
- `estlatent_robust`: two-step robust version (first-stage loading estimation, then regression; Godambe sandwich for variance). This is less efficient but more robust to model misspecification. See the example below.

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

# average-proxy GMM
fit_gmm_ave <- estlatent_ave(
  mod = "~Z+x1+x2",
  Z = "Z",
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
  Z = "Z",
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
- `mod = NULL` means treatment-only model (`~Z`).
- `IV_Y` is for loading moments in `method = "gmm"` only.
- For `method = "sem"`, `IV_Y` is not used; SEM is estimated from the full model specification directly.
- Regression moments in GMM automatically use covariates from `mod` as IVs.
- `IV_Y` may include variables from `data` even if they are not in `mod`.
- Use `summary(fit)` to print coefficient tables.


Model Misspecification:

In the following simulation, the data-generating process includes an additional latent component `eta_alt` that also depends on `X` and enters the second and third measurements:

```r
n <- 5000
X <- 3*rnorm(n)
ATE <- 2

eta_0 <- X + rnorm(n)
eta_1 <- ATE + X + rnorm(n)

eta_alt_0 <- X + rnorm(n)  # additional latent outcome
eta_alt_1 <- X + rnorm(n)

y1_0 <- 1.0*eta_0 + 1*rnorm(n)
y1_1 <- 1.0*eta_1 + 1*rnorm(n)

y2_0 <- 0.5*eta_0 + 2*rnorm(n) + 1*eta_alt_0 # additional latent outcome afects Y_2
y2_1 <- 0.5*eta_1 + 2*rnorm(n) + 1*eta_alt_1

y3_0 <- 2*eta_0 + 2*rnorm(n) + 2*eta_alt_0 # additional latent outcome afects Y_3
y3_1 <- 2*eta_1 + 2*rnorm(n) + 2*eta_alt_1

### treatment assignment

Z <- complete_ra(n)

Y1 <- y1_1*Z + (1-Z)*y1_0  # observed indicators
Y2 <- y2_1*Z + (1-Z)*y2_0
Y3 <- y3_1*Z + (1-Z)*y3_0

dat <- data.frame(Z,Y1,Y2,Y3,X)

```

In this setting, adding `X` to `mod` (`~Z+X`) can make single-step GMM (jointly estimating loading and regression moments) unstable for loadings, because `X`-related regression moments can pull loading estimates.

`estlatent_robust` is designed for this case:

- stage 1 identifies loadings using loading-IV moments (`IV_Y`, e.g. `Z`);
- stage 2 estimates regression coefficients conditional on stage-1 loadings;
- variance is adjusted with stacked estimating equations (Godambe sandwich).

So when you move from `mod="~Z"` to `mod="~Z+X"` under this type of DGP, `estlatent_robust(..., method="gmm")` is recommended. In practice, researchers can use `estlatent` first. If the result from `estlatent_robust` is very different from `estlatent`, and all IVs are valid, this suggests additional concern that covariates `X` may affect loading estimation in the joint GMM model. In that case, `estlatent_robust` is recommended.

If you have any comments, suggestions, or findings, please email me at: jiawei.fu@duke.edu
