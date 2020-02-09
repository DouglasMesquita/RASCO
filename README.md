
# RASCO: An R package to Alleviate the Spatial Confounding

<!-- badges: start -->

<!-- [![CRAN status](https://www.r-pkg.org/badges/version/FLAMES)](https://cran.r-project.org/package=FLAMES) -->

<!-- [![Travis build status](https://travis-ci.org/DouglasMesquita/FLAMES.svg?branch=master)](https://travis-ci.org/DouglasMesquita/FLAMES) -->

<!-- [![Codecov test coverage](https://codecov.io/gh/DouglasMesquita/FLAMES/branch/master/graph/badge.svg)](https://codecov.io/gh/DouglasMesquita/FLAMES?branch=master) -->

<!-- badges: end -->

## Overview

RASCO is an R package that allow the practioners to fit restricted
spatial models for three class of models: + Generalized linear mixed
models + Shared component models + Spatial frailty models

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
#> alpha1  0.488012  0.488569 0.048372  0.389522  0.574772
#> alpha2  0.142272  0.141245 0.049397  0.044328  0.232808
#> X11_1  -0.448335 -0.447102 0.070483 -0.580809 -0.309785
#> X12_1  -0.530059 -0.507789 0.450775 -1.429520  0.334260
#> X21_2  -0.795016 -0.794650 0.056349 -0.904521 -0.688811
#> X12_2  -0.300217 -0.298368 0.219307 -0.710260  0.157520
rscm_inla$summary_fixed
#>             mean    median       sd     lower     upper
#> alpha1  0.475652  0.476027 0.052161  0.376391  0.575305
#> alpha2  0.141362  0.141404 0.048442  0.042682  0.231837
#> X11_1  -0.438496 -0.439143 0.076733 -0.585092 -0.287717
#> X12_1  -0.375071 -0.378957 0.106843 -0.580092 -0.167523
#> X21_2  -0.740850 -0.738411 0.055293 -0.856597 -0.644021
#> X12_2  -0.372616 -0.370667 0.059332 -0.496261 -0.263359

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

We can notice that the restricted model corrects the variance inflation

``` r
SVIF(rscm_inla, scm_inla)
#>   parameter        VIF
#> 1    alpha1  0.8599957
#> 2    alpha2  1.0398172
#> 3     X11_1  0.8437318
#> 4     X12_1 17.8003072
#> 5     X21_2  1.0385613
#> 6     X12_2 13.6624005
```

## Restricted spatial frailty models

``` r
set.seed(1)

##-- Spatial structure
data("neigh_RJ")

##-- Individuals and regions
n_reg <- length(neigh_RJ)
n_id <- sample(x = 3:5, size = n_reg, replace = T)
coefs <- c(0.3, -0.3)
tau <- 0.75 # Scale of spatial effect

##-- Data
data <- rsurv(n_id = n_id,
              coefs = coefs, cens = 0.5, scale = FALSE,
              cens_type = "right", hazard = "weibull",
              hazard_params = hazard_params <- list(weibull = list(alpha = 1.2, variant = 0)),
              spatial = "ICAR",
              neigh = neigh_RJ, tau = tau, confounding = "linear", proj = "none")

##-- Models
weibull_inla <- rsfm(data = data, time = "L", status = "status",
                     covariates = c("X1", "X2"), intercept = TRUE,
                     family = "weibull", proj = "rhz", nsamp = 1000, approach = "inla")
#> Warning in inla.surv.1(time, event, time2, truncation): 'time2' is ignored for
#> data that are not interval censored

rsfm_inla <- rsfm(data = data, time = "L", status = "status", area = "reg",
                  covariates = c("X1", "X2"), intercept = TRUE,
                  model = "restricted_besag", neigh = neigh_RJ,
                  family = "weibull", proj = "rhz", nsamp = 1000, approach = "inla")
#> Warning in inla.surv.1(time, event, time2, truncation): 'time2' is ignored for
#> data that are not interval censored

weibull_inla$unrestricted$summary_fixed
#>                  mean    median       sd     lower     upper
#> (Intercept) -0.615542 -0.614855 0.076730 -0.755362 -0.459974
#> X1           0.230099  0.228612 0.074603  0.087270  0.374713
#> X2          -0.942590 -0.938899 0.096096 -1.120827 -0.759549
rsfm_inla$unrestricted$summary_fixed
#>                  mean    median       sd     lower     upper
#> (Intercept) -0.539812 -0.537958 0.080770 -0.697389 -0.385100
#> X1           0.317940  0.317355 0.078439  0.168108  0.472043
#> X2          -1.072193 -1.075586 0.352117 -1.781310 -0.416856
rsfm_inla$restricted$summary_fixed
#>                  mean    median       sd     lower     upper
#> (Intercept) -0.508432 -0.508401 0.079677 -0.672593 -0.362172
#> X1           0.332966  0.330292 0.075201  0.183030  0.481675
#> X2          -1.288779 -1.289738 0.113797 -1.520087 -1.074102
```

``` r
SVIF(weibull_inla$unrestricted, rsfm_inla$unrestricted)
#>     parameter       VIF
#> 1 (Intercept)  1.108077
#> 2          X1  1.105482
#> 3          X2 13.426517
SVIF(weibull_inla$unrestricted, rsfm_inla$restricted)
#>     parameter      VIF
#> 1 (Intercept) 1.078290
#> 2          X1 1.016096
#> 3          X2 1.402333
```
