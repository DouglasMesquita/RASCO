
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
library(spdep)

set.seed(1)
```

## Restricted generalized linear mixed models

``` r
##-- Spatial structure
data("neigh_RJ")

beta <- c(-0.5, -0.2)
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
#>                  mean    median       sd     lower     upper
#> (Intercept)  0.234140  0.233724 0.096053  0.051814  0.425890
#> X1          -0.363056 -0.359885 0.088340 -0.521581 -0.175938
#> X2          -0.039969 -0.041096 0.093569 -0.227393  0.140030
sglmm_mod$unrestricted$summary_fixed
#>                  mean    median       sd     lower     upper
#> (Intercept)  0.234726  0.233819 0.096558  0.046675  0.415938
#> X1          -0.360839 -0.359893 0.091211 -0.536051 -0.185585
#> X2          -0.037153 -0.036118 0.095402 -0.235432  0.135097
rglmm_mod$unrestricted$summary_fixed
#>                  mean    median       sd     lower     upper
#> (Intercept)  0.233733  0.230952 0.096278  0.043769  0.410482
#> X1          -0.361621 -0.360012 0.091000 -0.518872 -0.166869
#> X2          -0.038363 -0.038613 0.094865 -0.226378  0.152786
rglmm_mod$restricted$summary_fixed
#>                  mean    median       sd     lower     upper
#> (Intercept)  0.233733  0.230952 0.096278  0.043769  0.410482
#> X1          -0.361619 -0.360160 0.090986 -0.517856 -0.165391
#> X2          -0.038386 -0.039536 0.094249 -0.221058  0.161132

sglmm_mod$unrestricted$summary_hyperpar
#>                       mean   median       sd    lower    upper
#> Precision for reg 19877.01 13808.96 19544.57 25.59507 60062.62
```

The *SVIF* function provides the coefficients variance ratio between two
models. We can notice that the restricted model alleviates the variance
inflation

``` r
SVIF(rglmm_mod$restricted, rglmm_mod$unrestricted)
#>     parameter      VIF
#> 1 (Intercept) 1.000000
#> 2          X1 1.000308
#> 3          X2 1.013114
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
#> alpha1  0.444835  0.444290 0.049017  0.344603  0.534537
#> alpha2  0.114686  0.114271 0.051880  0.013314  0.214531
#> X11_1  -0.530009 -0.533021 0.078988 -0.683207 -0.381498
#> X12_1  -0.715365 -0.717143 0.473128 -1.650690  0.229479
#> X21_2  -0.784633 -0.785811 0.046222 -0.879345 -0.698487
#> X12_2  -0.790912 -0.792536 0.253973 -1.287329 -0.288363
rscm_mod$summary_fixed
#>             mean    median       sd     lower     upper
#> alpha1  0.418908  0.422525 0.052107  0.319380  0.518323
#> alpha2  0.126134  0.126813 0.050311  0.016541  0.215540
#> X11_1  -0.601610 -0.603692 0.085459 -0.764508 -0.444107
#> X12_1  -0.043488 -0.049713 0.122457 -0.273821  0.205812
#> X21_2  -0.764726 -0.765670 0.048473 -0.854432 -0.665676
#> X12_2  -0.381929 -0.381375 0.073681 -0.523138 -0.239222

scm_mod$summary_hyperpar
#>                           mean    median        sd    lower     upper
#> Precision for psi     2.047080  1.971893  0.579188 1.031142  3.210446
#> Precision for phi1   15.126391 11.275114 13.359959 0.801200 40.649718
#> Precision for phi2   19.102276 15.414731 13.359131 2.966504 45.169040
#> Beta for psi_gamma    2.049450  2.034844  0.258429 1.560730  2.566521
#> Delta                 1.429609  1.427436  0.090595 1.265242  1.613417
#> Precision for psi SC  0.990664  0.964977  0.213011 0.637999  1.444251
rscm_mod$summary_hyperpar
#>                           mean    median        sd    lower     upper
#> Precision for psi     0.693862  0.671965  0.165764 0.403133  1.030029
#> Precision for phi1   21.891229 16.672072 18.747190 0.903862 58.052288
#> Precision for phi2   22.920028 17.766901 18.460906 1.715648 58.374784
#> Beta for psi_gamma    1.990893  1.979407  0.193787 1.623614  2.378853
#> Delta                 1.411150  1.407817  0.066212 1.286138  1.530028
#> Precision for psi SC  0.346761  0.340935  0.066084 0.220676  0.471030
```

``` r
SVIF(rscm_mod, scm_mod)
#>   parameter        VIF
#> 1    alpha1  0.8849145
#> 2    alpha2  1.0633446
#> 3     X11_1  0.8542926
#> 4     X12_1 14.9276028
#> 5     X21_2  0.9092801
#> 6     X12_2 11.8813031
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
#> (Intercept) -0.568108 -0.572712 0.077053 -0.704936 -0.407553
#> X1           0.388676  0.388689 0.075055  0.251463  0.541067
#> X2          -0.025103 -0.023536 0.076707 -0.177105  0.120737
rsfm_mod$unrestricted$summary_fixed
#>                  mean    median       sd     lower     upper
#> (Intercept) -0.297759 -0.300678 0.078167 -0.445786 -0.143095
#> X1           0.450748  0.454387 0.100643  0.246986  0.648340
#> X2          -0.032407 -0.019305 0.542743 -1.034116  1.099458
rsfm_mod$restricted$summary_fixed
#>                  mean    median       sd     lower     upper
#> (Intercept) -0.336549 -0.336928 0.076812 -0.493728 -0.200647
#> X1           0.452597  0.453224 0.093384  0.272926  0.644396
#> X2          -0.248973 -0.247070 0.080098 -0.395440 -0.086972
```

``` r
SVIF(weibull_mod$unrestricted, rsfm_mod$unrestricted)
#>     parameter       VIF
#> 1 (Intercept)  1.029124
#> 2          X1  1.798075
#> 3          X2 50.063184
SVIF(weibull_mod$unrestricted, rsfm_mod$restricted)
#>     parameter       VIF
#> 1 (Intercept) 0.9937543
#> 2          X1 1.5480525
#> 3          X2 1.0903686
```
