
# RASCO: An R package to Alleviate the Spatial Confounding

<!-- badges: start -->

<!-- [![CRAN status](https://www.r-pkg.org/badges/version/FLAMES)](https://cran.r-project.org/package=FLAMES) -->

[![Travis build
status](https://travis-ci.org/DouglasMesquita/RASCO.svg?branch=master)](https://travis-ci.org/DouglasMesquita/RASCO)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/DouglasMesquita/RASCO?branch=master&svg=true)](https://ci.appveyor.com/project/DouglasMesquita/RASCO)
[![Codecov test
coverage](https://codecov.io/gh/DouglasMesquita/RASCO/branch/master/graph/badge.svg)](https://codecov.io/gh/DouglasMesquita/RASCO?branch=master)
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
#> (Intercept)  0.290110  0.292490 0.093569  0.100859 0.455292
#> X1          -0.111292 -0.112162 0.084320 -0.278974 0.048427
#> X2           0.428306  0.426124 0.093354  0.259872 0.622916
sglmm_mod$unrestricted$summary_fixed
#>                  mean    median       sd     lower    upper
#> (Intercept)  0.187974  0.192661 0.111995 -0.020120 0.426281
#> X1          -0.099619 -0.100398 0.093043 -0.289350 0.075949
#> X2           0.550247  0.525210 0.259846  0.084583 1.094470
rglmm_mod$unrestricted$summary_fixed
#>                  mean    median       sd     lower    upper
#> (Intercept)  0.195614  0.198452 0.109912 -0.027891 0.398912
#> X1          -0.105183 -0.105355 0.092058 -0.275318 0.088454
#> X2           0.576132  0.547001 0.267208  0.063450 1.136711
rglmm_mod$restricted$summary_fixed
#>                  mean    median       sd     lower    upper
#> (Intercept)  0.195614  0.198452 0.109912 -0.027891 0.398912
#> X1          -0.119478 -0.122081 0.085638 -0.282163 0.049612
#> X2           0.433391  0.436773 0.100427  0.227558 0.613935

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
#> 2          X1 1.155553
#> 3          X2 7.079424
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
                priors = list(prior_prec = c(0.5, 0.05), 
                              prior_gamma = c(0, 0.5)),
                proj = "none", nsamp = 1000,
                random_effects = list(shared = TRUE, specific_1 = TRUE, specific_2 = TRUE))
#> priors$prior_prec should be a list of tau_s, tau_1 and tau_2 entries.
#>                      Setting all priors as Gamma(0.5, 0.05).

rscm_mod <- rscm(data = data,
                 formula1 = Y1 ~ X11 + X12,
                 formula2 = Y2 ~ X21 + X12,
                 E1 = E1, E2 = E2,
                 family = c("poisson", "poisson"),
                 area = "reg", neigh = neigh_RJ,
                 priors = list(prior_prec = list(tau_s = c(0.5, 0.05), 
                                                 tau_1 = c(0.5, 0.05), 
                                                 tau_2 = c(0.5, 0.05)), 
                               prior_gamma = c(0, 0.5)),
                 proj = "spock", nsamp = 1000,
                 random_effects = list(shared = TRUE, specific_1 = TRUE, specific_2 = TRUE))

##-- Summary
scm_mod$summary_fixed
#>             mean    median       sd     lower     upper
#> alpha1  0.488101  0.489665 0.053849  0.388463  0.599593
#> alpha2  0.108551  0.109165 0.050546  0.009578  0.205503
#> X11_1  -0.476797 -0.476222 0.095401 -0.653551 -0.283833
#> X12_1  -0.582682 -0.596267 0.507903 -1.591821  0.375221
#> X21_2  -0.900741 -0.900096 0.055479 -1.012364 -0.797356
#> X12_2  -0.350218 -0.346017 0.221059 -0.840420  0.035665
rscm_mod$summary_fixed
#>             mean    median       sd     lower     upper
#> alpha1  0.526091  0.525487 0.050711  0.433985  0.622907
#> alpha2  0.116017  0.116075 0.053132  0.011229  0.216486
#> X11_1  -0.552423 -0.554767 0.095959 -0.729508 -0.370794
#> X12_1  -0.254576 -0.255629 0.107115 -0.478685 -0.064111
#> X21_2  -0.934883 -0.935380 0.049972 -1.046773 -0.847449
#> X12_2  -0.139080 -0.137913 0.055438 -0.242766 -0.032737

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
#> 1    alpha1  1.1275893
#> 2    alpha2  0.9050264
#> 3     X11_1  0.9884038
#> 4     X12_1 22.4833474
#> 5     X21_2  1.2325478
#> 6     X12_2 15.9001527
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
#> (Intercept) -0.458985 -0.456183 0.074568 -0.616267 -0.325996
#> X1           0.368564  0.369362 0.073571  0.227852  0.504217
#> X2          -0.652584 -0.653087 0.085142 -0.821409 -0.495560
rsfm_mod$unrestricted$summary_fixed
#>                  mean    median       sd     lower     upper
#> (Intercept) -0.430938 -0.429452 0.077328 -0.594945 -0.291085
#> X1           0.368916  0.370432 0.077681  0.201439  0.505782
#> X2          -0.653416 -0.662574 0.236759 -1.107735 -0.149682
rsfm_mod$restricted$summary_fixed
#>                  mean    median       sd     lower     upper
#> (Intercept) -0.434229 -0.434032 0.077096 -0.599067 -0.297298
#> X1           0.374846  0.377902 0.077744  0.225373  0.527650
#> X2          -0.687850 -0.691893 0.084236 -0.861404 -0.539015
```

``` r
SVIF(weibull_mod$unrestricted, rsfm_mod$unrestricted)
#>     parameter      VIF
#> 1 (Intercept) 1.075396
#> 2          X1 1.114850
#> 3          X2 7.732596
SVIF(weibull_mod$unrestricted, rsfm_mod$restricted)
#>     parameter       VIF
#> 1 (Intercept) 1.0689532
#> 2          X1 1.1166587
#> 3          X2 0.9788311
```
