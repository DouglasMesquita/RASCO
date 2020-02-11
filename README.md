
# RASCO: An R package to Alleviate the Spatial Confounding

<!-- badges: start -->

<!-- [![CRAN status](https://www.r-pkg.org/badges/version/FLAMES)](https://cran.r-project.org/package=FLAMES) -->

[![Travis build
status](https://travis-ci.org/DouglasMesquita/RASCO.svg?branch=master)](https://travis-ci.org/DouglasMesquita/RASCO)
<!-- [![Codecov test coverage](https://codecov.io/gh/DouglasMesquita/FLAMES/branch/master/graph/badge.svg)](https://codecov.io/gh/DouglasMesquita/FLAMES?branch=master) -->
<!-- badges: end -->

## Overview

RASCO is an R package that allow the practioners to fit restricted
spatial models for three class of models:

  - Generalized linear mixed models

  - Shared component models

  - Spatial frailty models

## Installation

``` r
# Install from CRAN (when available)
install.packages("RASCO")
# Or the development version from GitHub
# install.packages("devtools")
devtools::install_github("DouglasMesquita/RASCO")
```

``` r
library(RASCO)
library(spdep)

set.seed(11022020)
```

## Restricted generalized linear mixed models

``` r
##-- Spatial structure
data("neigh_RJ")

beta <- c(-0.1, 0.7)
tau <- 1

##-- Data ----
family <- "poisson"
data <- rglmm(beta = beta, tau = tau, family = family,
              confounding = "linear", neigh = neigh_RJ,
              scale = TRUE)

##-- Models ----
sglm_mod <- rsglmm(data = data, formula = Y ~ X1 + X2,
                   family = family,
                   proj = "none", nsamp = 1000)

sglmm_mod <- rsglmm(data = data, formula = Y ~ X1 + X2,
                    area = "reg", model = "besag", neigh = neigh_RJ,
                    family = family,
                    proj = "none", nsamp = 1000)

rglmm_mod <- rsglmm(data = data, formula = Y ~ X1 + X2,
                    area = "reg", model = "restricted_besag", neigh = neigh_RJ,
                    family = family,
                    proj = "rhz", nsamp = 1000)

sglm_mod$unrestricted$summary_fixed
#>                  mean    median       sd     lower    upper
#> (Intercept)  0.288384  0.288066 0.095445  0.099737 0.471440
#> X1          -0.110637 -0.111355 0.085908 -0.286088 0.043759
#> X2           0.428009  0.427437 0.094528  0.249244 0.619268
sglmm_mod$unrestricted$summary_fixed
#>                  mean    median       sd     lower    upper
#> (Intercept)  0.194440  0.194217 0.107994 -0.001995 0.423315
#> X1          -0.103173 -0.103763 0.094092 -0.280014 0.091976
#> X2           0.553270  0.525261 0.283259  0.005422 1.143170
rglmm_mod$unrestricted$summary_fixed
#>                  mean    median       sd     lower    upper
#> (Intercept)  0.196574  0.195831 0.111715 -0.023874 0.413492
#> X1          -0.102724 -0.103048 0.091587 -0.269686 0.083344
#> X2           0.549309  0.522097 0.268174 -0.008230 1.055000
rglmm_mod$restricted$summary_fixed
#>                  mean    median       sd     lower    upper
#> (Intercept)  0.196574  0.195831 0.111715 -0.023874 0.413492
#> X1          -0.118781 -0.118497 0.085893 -0.283233 0.043029
#> X2           0.435790  0.434563 0.097035  0.259315 0.628266

sglmm_mod$unrestricted$summary_hyperpar
#>                       mean   median       sd   lower    upper
#> Precision for reg 93.31269 2.227233 464.0461 0.40322 415.3316
```

The *SVIF* function provides the coefficients variance ratio between two
models. We can notice that the restricted model alleviates the variance
inflation

``` r
SVIF(rglmm_mod$restricted, rglmm_mod$unrestricted)
#>     parameter      VIF
#> 1 (Intercept) 1.000000
#> 2          X1 1.136978
#> 3          X2 7.637945
```

## Restricted shared component models

``` r
##-- Spatial structure
data("neigh_RJ")

##-- Parameters
alpha_1 <- 0.5
alpha_2 <- 0.1
beta_1 <- c(-0.5, -0.2)
beta_2 <- c(-0.8, -0.4)
tau_s <- 1
tau_1 <- tau_2 <- 10
delta <- 1.5

##-- Data
data <- rshared(alpha_1 = alpha_1, alpha_2 = alpha_2, beta_1 = beta_1, beta_2 = beta_2, delta = delta,
                tau_1 = tau_1, tau_2 = tau_2, tau_s = tau_s,
                confounding = "linear",
                neigh = neigh_RJ)

##-- Models
scm_mod <- rscm(data = data,
                formula1 = Y1 ~ X11 + X12,
                formula2 = Y2 ~ X21 + X12,
                E1 = E1, E2 = E2,
                family = c("poisson", "poisson"),
                area = "reg", neigh = neigh_RJ,
                prior_prec = c(0.5, 0.05), prior_gamma = c(0, 0.5),
                proj = "none", nsamp = 1000,
                random_effects = list(shared = TRUE, specific_1 = TRUE, specific_2 = TRUE))

rscm_mod <- rscm(data = data,
                 formula1 = Y1 ~ X11 + X12,
                 formula2 = Y2 ~ X21 + X12,
                 E1 = E1, E2 = E2,
                 family = c("poisson", "poisson"),
                 area = "reg", neigh = neigh_RJ,
                 prior_prec = c(0.5, 0.05), prior_gamma = c(0, 0.5),
                 proj = "spock", nsamp = 1000,
                 random_effects = list(shared = TRUE, specific_1 = TRUE, specific_2 = TRUE))

##-- Summary
scm_mod$summary_fixed
#>             mean    median       sd     lower     upper
#> alpha1  0.490313  0.490043 0.052476  0.394132  0.597933
#> alpha2  0.107959  0.107471 0.053962  0.009425  0.224420
#> X11_1  -0.474291 -0.473146 0.093431 -0.685415 -0.318965
#> X12_1  -0.585328 -0.573763 0.521342 -1.668230  0.399582
#> X21_2  -0.901234 -0.901356 0.055838 -1.005799 -0.791220
#> X12_2  -0.351207 -0.344064 0.223743 -0.800528  0.096918
rscm_mod$summary_fixed
#>             mean    median       sd     lower     upper
#> alpha1  0.527162  0.529550 0.051921  0.427655  0.625021
#> alpha2  0.115785  0.115874 0.050932  0.019730  0.215868
#> X11_1  -0.556725 -0.557658 0.092407 -0.731693 -0.374324
#> X12_1  -0.246953 -0.251505 0.104884 -0.456743 -0.046619
#> X21_2  -0.932813 -0.932291 0.048282 -1.018921 -0.834202
#> X12_2  -0.134994 -0.133983 0.054346 -0.244502 -0.036357

scm_mod$summary_hyperpar
#>                           mean    median       sd    lower     upper
#> Precision for psi     3.284081  3.108296 1.087815 1.447357  5.462632
#> Precision for phi1   12.022372  9.289885 9.805738 0.813852 30.863818
#> Precision for phi2   13.822226 11.077708 9.782992 2.211119 32.889614
#> Beta for psi_gamma    2.762142  2.735141 0.393967 2.018536  3.547350
#> Delta                 1.651258  1.643129 0.117720 1.424644  1.874252
#> Precision for psi SC  1.167597  1.142887 0.256402 0.721834  1.724771
rscm_mod$summary_hyperpar
#>                           mean    median        sd    lower     upper
#> Precision for psi     2.011603  1.895586  0.698430 0.848612  3.414483
#> Precision for phi1   13.182119  8.975717 13.737089 0.291429 39.161188
#> Precision for phi2   23.396619 18.056743 19.550683 0.960910 61.198257
#> Beta for psi_gamma    2.785578  2.749610  0.431860 1.982126  3.651760
#> Delta                 1.665025  1.662183  0.126486 1.416899  1.901444
#> Precision for psi SC  0.708885  0.694174  0.164693 0.444253  1.074248
```

``` r
SVIF(rscm_mod, scm_mod)
#>   parameter       VIF
#> 1    alpha1  1.021493
#> 2    alpha2  1.122521
#> 3     X11_1  1.022286
#> 4     X12_1 24.707394
#> 5     X21_2  1.337486
#> 6     X12_2 16.949767
```

## Restricted spatial frailty models

``` r
#-- Spatial structure
data("neigh_RJ")

##-- Individuals and regions
n_reg <- length(neigh_RJ)
n_id <- sample(x = 3:5, size = n_reg, replace = T)
beta <- c(0.3, -0.3)
tau <- 0.75 # Scale of spatial effect

##-- Data
data <- rsurv(n_id = n_id,
              coefs = beta, cens = 0.5, scale = FALSE,
              cens_type = "right", hazard = "weibull",
              hazard_params = hazard_params <- list(weibull = list(alpha = 1.2, variant = 0)),
              spatial = "ICAR",
              neigh = neigh_RJ, tau = tau, confounding = "linear", proj = "none")

##-- Models
weibull_mod <- rsfm(data = data,
                    formula = surv(time = L, event = status) ~ X1 + X2,
                    model = "none", family = "weibull",
                    proj = "rhz", nsamp = 1000, approach = "inla")

rsfm_mod <- rsfm(data = data, area = "reg",
                 formula = surv(time = L, event = status) ~ X1 + X2,
                 model = "restricted_besag", neigh = neigh_RJ, family = "weibull",
                 proj = "rhz", nsamp = 1000, approach = "inla")

weibull_mod$unrestricted$summary_fixed
#>                  mean    median       sd     lower     upper
#> (Intercept) -0.462781 -0.461182 0.074816 -0.603656 -0.321386
#> X1           0.370831  0.369983 0.075084  0.226143  0.516263
#> X2          -0.649794 -0.649874 0.087579 -0.821054 -0.476955
rsfm_mod$unrestricted$summary_fixed
#>                  mean    median       sd     lower     upper
#> (Intercept) -0.459245 -0.463886 0.074734 -0.606537 -0.312035
#> X1           0.368424  0.366112 0.075397  0.231683  0.518384
#> X2          -0.647384 -0.641536 0.089163 -0.825129 -0.494712
rsfm_mod$restricted$summary_fixed
#>                  mean    median       sd     lower     upper
#> (Intercept) -0.459227 -0.463886 0.074718 -0.606550 -0.311142
#> X1           0.368502  0.366044 0.075337  0.231973  0.516257
#> X2          -0.647098 -0.641864 0.087087 -0.816795 -0.490380
```

``` r
SVIF(weibull_mod$unrestricted, rsfm_mod$unrestricted)
#>     parameter       VIF
#> 1 (Intercept) 0.9978092
#> 2          X1 1.0083547
#> 3          X2 1.0365002
SVIF(weibull_mod$unrestricted, rsfm_mod$restricted)
#>     parameter      VIF
#> 1 (Intercept) 0.997382
#> 2          X1 1.006750
#> 3          X2 0.988796
```
