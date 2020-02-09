
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
#> alpha1  0.485808  0.486907 0.049308  0.393981  0.585722
#> alpha2  0.141786  0.141141 0.048648  0.043874  0.236382
#> X11_1  -0.444426 -0.443760 0.072838 -0.574318 -0.292077
#> X12_1  -0.546706 -0.547880 0.447061 -1.462246  0.234691
#> X21_2  -0.795418 -0.793623 0.056772 -0.909419 -0.684178
#> X12_2  -0.309800 -0.311300 0.214047 -0.717347  0.105813
rscm_inla$summary_fixed
#>             mean    median       sd     lower     upper
#> alpha1  0.474276  0.473831 0.051229  0.377303  0.574355
#> alpha2  0.140139  0.141570 0.049589  0.039709  0.231677
#> X11_1  -0.439000 -0.443620 0.080179 -0.602212 -0.294677
#> X12_1  -0.372206 -0.373594 0.105476 -0.585937 -0.172150
#> X21_2  -0.741666 -0.741454 0.055699 -0.843687 -0.625826
#> X12_2  -0.371649 -0.370018 0.059288 -0.485130 -0.256375

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
#> 1    alpha1  0.9264095
#> 2    alpha2  0.9624081
#> 3     X11_1  0.8252675
#> 4     X12_1 17.9649605
#> 5     X21_2  1.0388996
#> 6     X12_2 13.0342093
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

rsfm_inla <- rsfm(data = data, time = "L", status = "status", area = "reg",
                  covariates = c("X1", "X2"), intercept = TRUE,
                  model = "restricted_besag", neigh = neigh_RJ,
                  family = "weibull", proj = "rhz", nsamp = 1000, approach = "inla")

weibull_inla$unrestricted$summary_fixed
#>                  mean    median       sd     lower     upper
#> (Intercept) -0.611842 -0.608766 0.080801 -0.781361 -0.459660
#> X1           0.228842  0.231904 0.077503  0.079518  0.382685
#> X2          -0.942787 -0.947408 0.098957 -1.128245 -0.753624
rsfm_inla$unrestricted$summary_fixed
#>                  mean    median       sd     lower     upper
#> (Intercept) -0.539970 -0.537901 0.081628 -0.704039 -0.383428
#> X1           0.324066  0.326290 0.082192  0.170777  0.485429
#> X2          -1.066690 -1.056746 0.356526 -1.785210 -0.374065
rsfm_inla$restricted$summary_fixed
#>                  mean    median       sd     lower     upper
#> (Intercept) -0.508210 -0.505606 0.080053 -0.676535 -0.362712
#> X1           0.338540  0.338542 0.079803  0.179878  0.488082
#> X2          -1.288491 -1.294619 0.116528 -1.487207 -1.028042
```

``` r
SVIF(weibull_inla$unrestricted, rsfm_inla$unrestricted)
#>     parameter       VIF
#> 1 (Intercept)  1.020575
#> 2          X1  1.124662
#> 3          X2 12.980439
SVIF(weibull_inla$unrestricted, rsfm_inla$restricted)
#>     parameter       VIF
#> 1 (Intercept) 0.9815711
#> 2          X1 1.0602332
#> 3          X2 1.3866522
```
