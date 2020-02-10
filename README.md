
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
                 Y1 = "Y1", Y2 = "Y2",
                 X1 = c("X11", "X12"), X2 = c("X21", "X12"),
                 E1 = "E1", E2 = "E2",
                 area = "reg", neigh = neigh_RJ,
                 prior_prec = c(0.5, 0.05), prior_gamma = c(0, 0.5),
                 proj = "none")

rscm_inla <- rscm(data = data,
                  Y1 = "Y1", Y2 = "Y2",
                  X1 = c("X11", "X12"), X2 = c("X21", "X12"),
                  E1 = "E1", E2 = "E2",
                  area = "reg", neigh = neigh_RJ,
                  prior_prec = c(0.5, 0.05), prior_gamma = c(0, 0.5),
                  proj = "spock")

##-- Summary
scm_inla$summary_fixed
#>             mean    median       sd     lower     upper
#> alpha1  0.488412  0.489477 0.048235  0.394513  0.579293
#> alpha2  0.143232  0.142143 0.049261  0.050156  0.237546
#> X11_1  -0.443722 -0.442796 0.073942 -0.588624 -0.295900
#> X12_1  -0.562739 -0.563071 0.451748 -1.363277  0.399427
#> X21_2  -0.794629 -0.794316 0.055399 -0.893619 -0.679871
#> X12_2  -0.315873 -0.321550 0.218762 -0.756360  0.116906
rscm_inla$summary_fixed
#>             mean    median       sd     lower     upper
#> alpha1  0.478246  0.479306 0.051048  0.379235  0.575073
#> alpha2  0.142061  0.140984 0.048946  0.048625  0.236406
#> X11_1  -0.441951 -0.443142 0.078288 -0.595883 -0.291657
#> X12_1  -0.373009 -0.371430 0.103388 -0.563260 -0.170196
#> X21_2  -0.742390 -0.742596 0.055199 -0.856003 -0.643758
#> X12_2  -0.367838 -0.367043 0.060279 -0.476309 -0.243539

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
#> 1    alpha1  0.8928266
#> 2    alpha2  1.0129127
#> 3     X11_1  0.8920557
#> 4     X12_1 19.0920344
#> 5     X21_2  1.0072596
#> 6     X12_2 13.1707860
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
#> (Intercept) -0.631158 -0.631809 0.077973 -0.770568 -0.472022
#> X1           0.218495  0.215794 0.074147  0.067264  0.351719
#> X2          -0.953558 -0.957676 0.100043 -1.140607 -0.735529
rsfm_inla$unrestricted$summary_fixed
#>                  mean    median       sd     lower     upper
#> (Intercept) -0.545976 -0.544101 0.080899 -0.704822 -0.385425
#> X1           0.325475  0.329555 0.085392  0.156862  0.492125
#> X2          -1.114874 -1.106152 0.391321 -1.954553 -0.401492
rsfm_inla$restricted$summary_fixed
#>                  mean    median       sd     lower     upper
#> (Intercept) -0.511915 -0.510954 0.079607 -0.676323 -0.364322
#> X1           0.340865  0.340927 0.082559  0.178246  0.504190
#> X2          -1.362622 -1.362327 0.117976 -1.582961 -1.117825
```

``` r
SVIF(weibull_inla$unrestricted, rsfm_inla$unrestricted)
#>     parameter       VIF
#> 1 (Intercept)  1.076460
#> 2          X1  1.326317
#> 3          X2 15.300052
SVIF(weibull_inla$unrestricted, rsfm_inla$restricted)
#>     parameter      VIF
#> 1 (Intercept) 1.042351
#> 2          X1 1.239772
#> 3          X2 1.390637
```
