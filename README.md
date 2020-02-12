
# RASCO: An R package to Alleviate the Spatial Confounding

<!-- badges: start -->

<!-- [![CRAN status](https://www.r-pkg.org/badges/version/FLAMES)](https://cran.r-project.org/package=FLAMES) -->

[![Travis build
status](https://travis-ci.org/DouglasMesquita/RASCO.svg?branch=master)](https://travis-ci.org/DouglasMesquita/RASCO)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/DouglasMesquita/RASCO?branch=master&svg=true)](https://ci.appveyor.com/project/DouglasMesquita/RASCO)
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
#> (Intercept)  0.287840  0.286986 0.088930  0.128059 0.466437
#> X1          -0.107546 -0.107553 0.087970 -0.270789 0.068594
#> X2           0.424493  0.421990 0.089728  0.243294 0.599016
sglmm_mod$unrestricted$summary_fixed
#>                  mean    median       sd     lower    upper
#> (Intercept)  0.195514  0.200821 0.107322 -0.026039 0.389527
#> X1          -0.103671 -0.099070 0.091310 -0.291783 0.058065
#> X2           0.550647  0.533660 0.263821  0.026987 1.104850
rglmm_mod$unrestricted$summary_fixed
#>                  mean    median       sd     lower    upper
#> (Intercept)  0.194461  0.197492 0.114990 -0.032533 0.420285
#> X1          -0.104303 -0.105881 0.093481 -0.299888 0.070544
#> X2           0.548703  0.519262 0.274365  0.049196 1.143555
rglmm_mod$restricted$summary_fixed
#>                  mean    median       sd     lower    upper
#> (Intercept)  0.194461  0.197492 0.114990 -0.032533 0.420285
#> X1          -0.120101 -0.120904 0.087624 -0.296421 0.048025
#> X2           0.428989  0.430657 0.101304  0.227754 0.618744

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
#> 2          X1 1.138153
#> 3          X2 7.335069
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
#> alpha1  0.491486  0.491392 0.052139  0.384812  0.587316
#> alpha2  0.110576  0.109214 0.052267  0.015981  0.219435
#> X11_1  -0.471922 -0.470421 0.094502 -0.646395 -0.277691
#> X12_1  -0.603204 -0.598649 0.506690 -1.496355  0.479589
#> X21_2  -0.899790 -0.898724 0.053795 -1.003211 -0.791118
#> X12_2  -0.361042 -0.356485 0.226431 -0.795899  0.105931
rscm_mod$summary_fixed
#>             mean    median       sd     lower     upper
#> alpha1  0.527391  0.527986 0.049138  0.429674  0.621838
#> alpha2  0.116545  0.115311 0.052402  0.007486  0.211752
#> X11_1  -0.556035 -0.557181 0.095442 -0.733661 -0.366424
#> X12_1  -0.256408 -0.257747 0.112656 -0.465997 -0.016378
#> X21_2  -0.933321 -0.932266 0.050189 -1.034959 -0.833360
#> X12_2  -0.138244 -0.139146 0.054773 -0.256340 -0.038614

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
#> 1    alpha1  1.1258757
#> 2    alpha2  0.9948542
#> 3     X11_1  0.9803992
#> 4     X12_1 20.2290744
#> 5     X21_2  1.1488590
#> 6     X12_2 17.0898681
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
#> (Intercept) -0.459192 -0.459147 0.074733 -0.614491 -0.320534
#> X1           0.369517  0.369172 0.072660  0.219651  0.505794
#> X2          -0.651012 -0.652905 0.085803 -0.815979 -0.482628
rsfm_mod$unrestricted$summary_fixed
#>                  mean    median       sd     lower     upper
#> (Intercept) -0.434588 -0.431110 0.077017 -0.585203 -0.280597
#> X1           0.369686  0.367476 0.074903  0.220007  0.514144
#> X2          -0.673414 -0.676832 0.236477 -1.184048 -0.239608
rsfm_mod$restricted$summary_fixed
#>                  mean    median       sd     lower     upper
#> (Intercept) -0.437905 -0.434835 0.076374 -0.595776 -0.295067
#> X1           0.375699  0.374915 0.073978  0.244355  0.535687
#> X2          -0.692499 -0.690008 0.084389 -0.862317 -0.538028
```

``` r
SVIF(weibull_mod$unrestricted, rsfm_mod$unrestricted)
#>     parameter      VIF
#> 1 (Intercept) 1.062058
#> 2          X1 1.062693
#> 3          X2 7.595789
SVIF(weibull_mod$unrestricted, rsfm_mod$restricted)
#>     parameter       VIF
#> 1 (Intercept) 1.0443985
#> 2          X1 1.0366076
#> 3          X2 0.9673124
```
