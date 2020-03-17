
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
                    area = "reg", model = "r_besag", neigh = neigh_RJ,
                    family = family,
                    proj = "rhz", nsamp = 1000)

sglm_mod$unrestricted$summary_fixed
#>                  mean    median       sd     lower    upper
#> (Intercept)  0.287226  0.288730 0.092605  0.118322 0.472360
#> X1          -0.106900 -0.106967 0.085543 -0.266927 0.055962
#> X2           0.427531  0.428168 0.090722  0.246916 0.587882
sglmm_mod$unrestricted$summary_fixed
#>                  mean    median       sd     lower    upper
#> (Intercept)  0.195157  0.191272 0.110204 -0.047175 0.380089
#> X1          -0.103344 -0.102715 0.093884 -0.279431 0.075227
#> X2           0.559350  0.530745 0.264203  0.063551 1.105293
rglmm_mod$unrestricted$summary_fixed
#>                  mean    median       sd     lower    upper
#> (Intercept)  0.192668  0.192736 0.116342 -0.053419 0.395670
#> X1          -0.101675 -0.104354 0.094498 -0.268667 0.100495
#> X2           0.543842  0.520862 0.272767 -0.013356 1.068388
rglmm_mod$restricted$summary_fixed
#>                  mean    median       sd     lower    upper
#> (Intercept)  0.192668  0.192736 0.116342 -0.053419 0.395670
#> X1          -0.115735 -0.113369 0.087372 -0.278361 0.064544
#> X2           0.436432  0.438264 0.098314  0.232467 0.619361

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
#> 2          X1 1.169771
#> 3          X2 7.697557
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
#> alpha1  0.491614  0.491195 0.054164  0.393997  0.609417
#> alpha2  0.106997  0.108629 0.051381 -0.000913  0.199187
#> X11_1  -0.476352 -0.478711 0.088874 -0.657020 -0.301415
#> X12_1  -0.617828 -0.612019 0.525844 -1.673994  0.401673
#> X21_2  -0.902219 -0.903144 0.053276 -0.998208 -0.790118
#> X12_2  -0.365149 -0.361666 0.226272 -0.802082  0.082989
rscm_mod$summary_fixed
#>             mean    median       sd     lower     upper
#> alpha1  0.527345  0.526665 0.051598  0.422537  0.627027
#> alpha2  0.110543  0.110613 0.048732  0.014524  0.203526
#> X11_1  -0.556610 -0.555476 0.099379 -0.749712 -0.362264
#> X12_1  -0.259005 -0.261003 0.108517 -0.474078 -0.052476
#> X21_2  -0.939981 -0.940507 0.051287 -1.038394 -0.845972
#> X12_2  -0.139469 -0.141918 0.054965 -0.246342 -0.031624

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
#>   parameter       VIF
#> 1    alpha1  1.101934
#> 2    alpha2  1.111672
#> 3     X11_1  0.799761
#> 4     X12_1 23.481092
#> 5     X21_2  1.079068
#> 6     X12_2 16.946857
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
                 model = "r_besag", neigh = neigh_RJ, family = "weibull",
                 proj = "rhz", nsamp = 1000, approach = "inla")

weibull_mod$unrestricted$summary_fixed
#>                  mean    median       sd     lower     upper
#> (Intercept) -0.458072 -0.459351 0.074184 -0.611280 -0.322398
#> X1           0.367851  0.366297 0.076857  0.229124  0.526923
#> X2          -0.649918 -0.647815 0.087821 -0.823269 -0.487666
rsfm_mod$unrestricted$summary_fixed
#>                  mean    median       sd     lower     upper
#> (Intercept) -0.433391 -0.436618 0.075309 -0.572904 -0.268757
#> X1           0.368563  0.365297 0.074889  0.233098  0.527644
#> X2          -0.665535 -0.664599 0.230197 -1.151567 -0.227115
rsfm_mod$restricted$summary_fixed
#>                  mean    median       sd     lower     upper
#> (Intercept) -0.436251 -0.439294 0.074865 -0.579336 -0.281124
#> X1           0.373529  0.371103 0.073830  0.230856  0.519180
#> X2          -0.691434 -0.690229 0.086667 -0.853196 -0.530846
```

``` r
SVIF(weibull_mod$unrestricted, rsfm_mod$unrestricted)
#>     parameter       VIF
#> 1 (Intercept) 1.0305600
#> 2          X1 0.9494437
#> 3          X2 6.8707249
SVIF(weibull_mod$unrestricted, rsfm_mod$restricted)
#>     parameter       VIF
#> 1 (Intercept) 1.0184440
#> 2          X1 0.9227815
#> 3          X2 0.9738919
```
