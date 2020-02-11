
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
#> (Intercept)  0.289335  0.289671 0.090315  0.110202 0.457302
#> X1          -0.109400 -0.108917 0.083468 -0.254921 0.074609
#> X2           0.429735  0.427567 0.094172  0.247935 0.611579
sglmm_mod$unrestricted$summary_fixed
#>                  mean    median       sd     lower    upper
#> (Intercept)  0.195294  0.194610 0.113238 -0.030374 0.403536
#> X1          -0.099749 -0.097084 0.093392 -0.280571 0.072265
#> X2           0.556317  0.527599 0.270054  0.105365 1.175044
rglmm_mod$unrestricted$summary_fixed
#>                  mean    median       sd     lower    upper
#> (Intercept)  0.194258  0.196741 0.112036 -0.001281 0.427242
#> X1          -0.104612 -0.103433 0.093364 -0.296806 0.066523
#> X2           0.546118  0.524324 0.276618  0.058632 1.191173
rglmm_mod$restricted$summary_fixed
#>                  mean    median       sd     lower    upper
#> (Intercept)  0.194258  0.196741 0.112036 -0.001281 0.427242
#> X1          -0.119187 -0.118705 0.084993 -0.288867 0.048484
#> X2           0.433062  0.434083 0.100123  0.239147 0.620994

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
#> 2          X1 1.206681
#> 3          X2 7.632963
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
#> alpha1  0.489328  0.488564 0.052356  0.389638  0.581460
#> alpha2  0.113828  0.113739 0.051601  0.016677  0.212021
#> X11_1  -0.479250 -0.475099 0.094742 -0.660327 -0.303922
#> X12_1  -0.566737 -0.573537 0.516605 -1.543959  0.448547
#> X21_2  -0.897951 -0.898970 0.053831 -0.998165 -0.792814
#> X12_2  -0.346984 -0.345222 0.222795 -0.774976  0.111338
rscm_mod$summary_fixed
#>             mean    median       sd     lower     upper
#> alpha1  0.527536  0.528523 0.050101  0.427762  0.619014
#> alpha2  0.114563  0.113566 0.052067  0.009696  0.209329
#> X11_1  -0.555447 -0.555081 0.095758 -0.754289 -0.381246
#> X12_1  -0.257609 -0.256567 0.110835 -0.447862 -0.009351
#> X21_2  -0.933821 -0.933750 0.049359 -1.033968 -0.839445
#> X12_2  -0.139817 -0.141468 0.053293 -0.239519 -0.034900

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
#>   parameter        VIF
#> 1    alpha1  1.0920440
#> 2    alpha2  0.9821801
#> 3     X11_1  0.9788924
#> 4     X12_1 21.7251788
#> 5     X21_2  1.1894117
#> 6     X12_2 17.4771472
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
#> (Intercept) -0.460212 -0.456574 0.078705 -0.616012 -0.314946
#> X1           0.367287  0.366657 0.075775  0.232325  0.523148
#> X2          -0.650069 -0.650063 0.086333 -0.818502 -0.492382
rsfm_mod$unrestricted$summary_fixed
#>                  mean    median       sd     lower     upper
#> (Intercept) -0.459556 -0.458014 0.076660 -0.604354 -0.302516
#> X1           0.366251  0.367296 0.075687  0.219687  0.507452
#> X2          -0.649631 -0.647690 0.089157 -0.835462 -0.489533
rsfm_mod$restricted$summary_fixed
#>                  mean    median       sd     lower     upper
#> (Intercept) -0.459575 -0.457969 0.076650 -0.605107 -0.303576
#> X1           0.366359  0.367269 0.075578  0.219337  0.506583
#> X2          -0.650046 -0.648947 0.084869 -0.821640 -0.490096
```

``` r
SVIF(weibull_mod$unrestricted, rsfm_mod$unrestricted)
#>     parameter       VIF
#> 1 (Intercept) 0.9487089
#> 2          X1 0.9976787
#> 3          X2 1.0664911
SVIF(weibull_mod$unrestricted, rsfm_mod$restricted)
#>     parameter       VIF
#> 1 (Intercept) 0.9484614
#> 2          X1 0.9948072
#> 3          X2 0.9663724
```
