
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
#> alpha1  0.488777  0.487693 0.049311  0.390932  0.579200
#> alpha2  0.142209  0.141845 0.048368  0.045552  0.233823
#> X11_1  -0.445002 -0.442557 0.072723 -0.578096 -0.295945
#> X12_1  -0.544089 -0.548066 0.460966 -1.391055  0.408615
#> X21_2  -0.794244 -0.793939 0.054391 -0.894604 -0.679798
#> X12_2  -0.310322 -0.314814 0.222780 -0.787786  0.082211
rscm_inla$summary_fixed
#>             mean    median       sd     lower     upper
#> alpha1  0.472144  0.472240 0.051107  0.366060  0.563171
#> alpha2  0.139478  0.138595 0.047901  0.047550  0.234562
#> X11_1  -0.443697 -0.443476 0.076266 -0.589910 -0.291900
#> X12_1  -0.372233 -0.373023 0.108394 -0.569993 -0.149101
#> X21_2  -0.740820 -0.739934 0.054355 -0.845104 -0.632694
#> X12_2  -0.369958 -0.370070 0.059483 -0.485250 -0.252443

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
#> 1    alpha1  0.9309510
#> 2    alpha2  1.0195936
#> 3     X11_1  0.9092465
#> 4     X12_1 18.0853660
#> 5     X21_2  1.0013251
#> 6     X12_2 14.0270605
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
#> (Intercept) -0.617524 -0.617350 0.078713 -0.761019 -0.458617
#> X1           0.227056  0.231913 0.071525  0.079568  0.357547
#> X2          -0.937586 -0.940347 0.101809 -1.112742 -0.711154
rsfm_inla$unrestricted$summary_fixed
#>                  mean    median       sd     lower     upper
#> (Intercept) -0.540431 -0.539675 0.083703 -0.701963 -0.368178
#> X1           0.319638  0.323012 0.081655  0.165675  0.478105
#> X2          -1.054561 -1.047469 0.359264 -1.786012 -0.409842
rsfm_inla$restricted$summary_fixed
#>                  mean    median       sd     lower     upper
#> (Intercept) -0.508318 -0.508820 0.082508 -0.681266 -0.350496
#> X1           0.333549  0.335730 0.077672  0.178660  0.476266
#> X2          -1.290920 -1.287474 0.115016 -1.517228 -1.071286
```

``` r
SVIF(weibull_inla$unrestricted, rsfm_inla$unrestricted)
#>     parameter       VIF
#> 1 (Intercept)  1.130809
#> 2          X1  1.303316
#> 3          X2 12.452457
SVIF(weibull_inla$unrestricted, rsfm_inla$restricted)
#>     parameter      VIF
#> 1 (Intercept) 1.098751
#> 2          X1 1.179270
#> 3          X2 1.276275
```
