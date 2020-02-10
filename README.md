
# RASCO: An R package to Alleviate the Spatial Confounding

<!-- badges: start -->

<!-- [![CRAN status](https://www.r-pkg.org/badges/version/FLAMES)](https://cran.r-project.org/package=FLAMES) -->

<!-- [![Travis build status](https://travis-ci.org/DouglasMesquita/FLAMES.svg?branch=master)](https://travis-ci.org/DouglasMesquita/FLAMES) -->

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
```

## Restricted shared component models

``` r
library(spdep)

set.seed(1)

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
scm_inla <- rscm(data = data,
                 formula1 = Y1 ~ X11 + X12,
                 formula2 = Y2 ~ X21 + X12,
                 E1 = E1, E2 = E2,
                 family = c("poisson", "poisson"),
                 area = "reg", neigh = neigh_RJ,
                 prior_prec = c(0.5, 0.05), prior_gamma = c(0, 0.5),
                 proj = "none", nsamp = 1000,
                 random_effects = list(shared = TRUE, specific_1 = TRUE, specific_2 = TRUE))

rscm_inla <- rscm(data = data,
                  formula1 = Y1 ~ X11 + X12,
                  formula2 = Y2 ~ X21 + X12,
                  E1 = E1, E2 = E2,
                  family = c("poisson", "poisson"),
                  area = "reg", neigh = neigh_RJ,
                  prior_prec = c(0.5, 0.05), prior_gamma = c(0, 0.5),
                  proj = "spock", nsamp = 1000,
                  random_effects = list(shared = TRUE, specific_1 = TRUE, specific_2 = TRUE))

##-- Summary
scm_inla$summary_fixed
#>             mean    median       sd     lower     upper
#> alpha1  0.488155  0.488093 0.047341  0.396239  0.584773
#> alpha2  0.139756  0.138907 0.048688  0.042943  0.232202
#> X11_1  -0.444870 -0.445402 0.071095 -0.575126 -0.301725
#> X12_1  -0.556883 -0.563534 0.459428 -1.409031  0.341550
#> X21_2  -0.798265 -0.799240 0.054362 -0.895003 -0.690848
#> X12_2  -0.315622 -0.326434 0.228091 -0.739571  0.172519
rscm_inla$summary_fixed
#>             mean    median       sd     lower     upper
#> alpha1  0.472469  0.475138 0.051145  0.375071  0.568643
#> alpha2  0.140067  0.139617 0.050419  0.041303  0.238954
#> X11_1  -0.441870 -0.438864 0.077742 -0.589851 -0.292732
#> X12_1  -0.371793 -0.374852 0.108077 -0.585048 -0.175539
#> X21_2  -0.741650 -0.741494 0.054993 -0.844736 -0.636424
#> X12_2  -0.369193 -0.369934 0.059706 -0.489755 -0.259909

scm_inla$summary_hyperpar
#>                           mean    median        sd    lower     upper
#> Precision for psi     2.782194  2.649246  0.883663 1.270608  4.553980
#> Precision for phi1   13.753305 10.339682 11.959767 0.729111 36.651889
#> Precision for phi2   16.857693 13.528922 11.886939 2.712245 40.026533
#> Beta for psi_gamma    2.329956  2.306823  0.344941 1.678362  3.016666
#> Delta                 1.528520  1.523282  0.114793 1.326475  1.767456
#> Precision for psi SC  1.176645  1.149637  0.252747 0.735257  1.709944
rscm_inla$summary_hyperpar
#>                           mean    median        sd    lower     upper
#> Precision for psi     1.084122  1.040748  0.291365 0.581587  1.672468
#> Precision for phi1   16.095214 11.557751 15.369736 0.503655 45.443979
#> Precision for phi2   28.534089 23.682577 19.891185 2.627656 67.384196
#> Beta for psi_gamma    2.137382  2.117494  0.262393 1.645071  2.663296
#> Delta                 1.458065  1.452159  0.086622 1.301890  1.625117
#> Precision for psi SC  0.503221  0.490751  0.099809 0.336539  0.703725
```

The *SVIF* function provides the coefficients variance ratio between two
models. We can notice that the restricted model alleviates the variance
inflation

``` r
SVIF(rscm_inla, scm_inla)
#>   parameter        VIF
#> 1    alpha1  0.8567784
#> 2    alpha2  0.9325141
#> 3     X11_1  0.8363089
#> 4     X12_1 18.0704246
#> 5     X21_2  0.9771833
#> 6     X12_2 14.5942017
```

## Restricted spatial frailty models

``` r
set.seed(1)

##-- Spatial structure
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
weibull_inla <- rsfm(data = data,
                     formula = surv(time = L, event = status) ~ X1 + X2,
                     model = "none", family = "weibull",
                     proj = "rhz", nsamp = 1000, approach = "inla")

rsfm_inla <- rsfm(data = data, area = "reg",
                  formula = surv(time = L, event = status) ~ X1 + X2,
                  model = "restricted_besag", neigh = neigh_RJ, family = "weibull",
                  proj = "rhz", nsamp = 1000, approach = "inla")

weibull_inla$unrestricted$summary_fixed
#>                  mean    median       sd     lower     upper
#> (Intercept) -0.632758 -0.634013 0.080221 -0.786814 -0.480236
#> X1           0.221983  0.224242 0.071727  0.075577  0.351124
#> X2          -0.957781 -0.957047 0.099364 -1.130800 -0.742370
rsfm_inla$unrestricted$summary_fixed
#>                  mean    median       sd     lower     upper
#> (Intercept) -0.546252 -0.549211 0.078801 -0.721317 -0.414075
#> X1           0.319382  0.316741 0.082555  0.158015  0.470698
#> X2          -1.119702 -1.102745 0.394210 -1.842930 -0.343968
rsfm_inla$restricted$summary_fixed
#>                  mean    median       sd     lower     upper
#> (Intercept) -0.511496 -0.515569 0.077444 -0.678699 -0.371490
#> X1           0.334976  0.334862 0.078474  0.191850  0.497439
#> X2          -1.359038 -1.358118 0.115761 -1.601813 -1.154402
```

``` r
SVIF(weibull_inla$unrestricted, rsfm_inla$unrestricted)
#>     parameter        VIF
#> 1 (Intercept)  0.9649111
#> 2          X1  1.3247119
#> 3          X2 15.7397250
SVIF(weibull_inla$unrestricted, rsfm_inla$restricted)
#>     parameter       VIF
#> 1 (Intercept) 0.9319646
#> 2          X1 1.1969782
#> 3          X2 1.3572705
```
