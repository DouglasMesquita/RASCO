
# RASCO: An R package to Alleviate the Spatial Confounding

<!-- badges: start -->

<!-- [![CRAN status](https://www.r-pkg.org/badges/version/FLAMES)](https://cran.r-project.org/package=FLAMES) -->

[![Travis build
status](https://travis-ci.org/DouglasMesquita/RASCO.svg?branch=master)](https://travis-ci.org/DouglasMesquita/RASCO)
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
#> (Intercept)  0.290505  0.288100 0.093411  0.129896 0.497969
#> X1          -0.106084 -0.103652 0.084485 -0.263661 0.063578
#> X2           0.423949  0.426745 0.093956  0.237774 0.596724
sglmm_mod$unrestricted$summary_fixed
#>                  mean    median       sd     lower    upper
#> (Intercept)  0.196164  0.198939 0.111103 -0.023160 0.411708
#> X1          -0.100264 -0.097855 0.096350 -0.290509 0.078318
#> X2           0.545708  0.523933 0.275050  0.016044 1.111274
rglmm_mod$unrestricted$summary_fixed
#>                  mean    median       sd     lower    upper
#> (Intercept)  0.197078  0.199832 0.112016 -0.017172 0.407581
#> X1          -0.097194 -0.095851 0.089325 -0.258042 0.097413
#> X2           0.566850  0.540171 0.270261  0.091666 1.204857
rglmm_mod$restricted$summary_fixed
#>                  mean    median       sd     lower    upper
#> (Intercept)  0.197078  0.199832 0.112016 -0.017172 0.407581
#> X1          -0.115910 -0.115383 0.082340 -0.293573 0.033782
#> X2           0.438248  0.440041 0.094824  0.232093 0.602719

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
#> 2          X1 1.176859
#> 3          X2 8.123257
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
#> alpha1  0.489924  0.490364 0.051656  0.399211  0.599754
#> alpha2  0.106686  0.106732 0.051504  0.006322  0.205167
#> X11_1  -0.474191 -0.475857 0.093691 -0.639875 -0.279506
#> X12_1  -0.571173 -0.567986 0.509287 -1.507656  0.453113
#> X21_2  -0.902860 -0.901992 0.052857 -1.017253 -0.812152
#> X12_2  -0.351371 -0.351164 0.215107 -0.785360  0.046951
rscm_mod$summary_fixed
#>             mean    median       sd     lower     upper
#> alpha1  0.525635  0.525204 0.051292  0.424096  0.618425
#> alpha2  0.113613  0.116184 0.050411  0.009800  0.209775
#> X11_1  -0.556932 -0.559174 0.098462 -0.730937 -0.339811
#> X12_1  -0.256256 -0.254786 0.114594 -0.501099 -0.049720
#> X21_2  -0.937573 -0.937245 0.050177 -1.044410 -0.842962
#> X12_2  -0.136420 -0.135977 0.055455 -0.248657 -0.028354

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
#> 1    alpha1  1.0142436
#> 2    alpha2  1.0438337
#> 3     X11_1  0.9054374
#> 4     X12_1 19.7515609
#> 5     X21_2  1.1096746
#> 6     X12_2 15.0462290
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
#> (Intercept) -0.460805 -0.461169 0.073425 -0.600938 -0.320121
#> X1           0.368488  0.369997 0.071175  0.229249  0.499219
#> X2          -0.647891 -0.646635 0.083315 -0.815971 -0.492603
rsfm_mod$unrestricted$summary_fixed
#>                  mean    median       sd     lower     upper
#> (Intercept) -0.426102 -0.426004 0.077148 -0.574345 -0.269414
#> X1           0.365454  0.364447 0.075787  0.226314  0.512000
#> X2          -0.662287 -0.657237 0.233928 -1.124656 -0.213650
rsfm_mod$restricted$summary_fixed
#>                  mean    median       sd     lower     upper
#> (Intercept) -0.428896 -0.427880 0.076157 -0.578690 -0.277932
#> X1           0.370961  0.371268 0.074667  0.229603  0.510738
#> X2          -0.687778 -0.689013 0.083782 -0.854702 -0.532128
```

``` r
SVIF(weibull_mod$unrestricted, rsfm_mod$unrestricted)
#>     parameter      VIF
#> 1 (Intercept) 1.103981
#> 2          X1 1.133795
#> 3          X2 7.883481
SVIF(weibull_mod$unrestricted, rsfm_mod$restricted)
#>     parameter      VIF
#> 1 (Intercept) 1.075801
#> 2          X1 1.100531
#> 3          X2 1.011242
```
