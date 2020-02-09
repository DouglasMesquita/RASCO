
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

## Shared component models

``` r
library(RASCO)
library(spdep)
#> Loading required package: sp
#> Loading required package: spData
#> To access larger datasets in this package, install the spDataLarge
#> package with: `install.packages('spDataLarge',
#> repos='https://nowosad.github.io/drat/', type='source')`
#> Loading required package: sf
#> Linking to GEOS 3.6.2, GDAL 2.2.3, PROJ 4.9.3

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
#> alpha1  0.487311  0.490186 0.051162  0.390259  0.586618
#> alpha2  0.139328  0.138687 0.050037  0.047309  0.235870
#> X11_1  -0.444688 -0.442888 0.072485 -0.586889 -0.305352
#> X12_1  -0.563873 -0.572522 0.476051 -1.548699  0.327136
#> X21_2  -0.795867 -0.795942 0.056983 -0.901921 -0.680412
#> X12_2  -0.311263 -0.315459 0.222657 -0.728190  0.128891
rscm_inla$summary_fixed
#>             mean    median       sd     lower     upper
#> alpha1  0.474940  0.475680 0.048930  0.361073  0.555738
#> alpha2  0.140155  0.141711 0.050559  0.035038  0.230865
#> X11_1  -0.441265 -0.439722 0.077226 -0.594094 -0.287478
#> X12_1  -0.367949 -0.366770 0.102884 -0.575222 -0.157409
#> X21_2  -0.742207 -0.743272 0.055695 -0.846905 -0.631465
#> X12_2  -0.368564 -0.369143 0.059262 -0.480126 -0.248888

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
#> 1    alpha1  1.0933132
#> 2    alpha2  0.9794575
#> 3     X11_1  0.8809864
#> 4     X12_1 21.4097345
#> 5     X21_2  1.0467867
#> 6     X12_2 14.1162746
```
