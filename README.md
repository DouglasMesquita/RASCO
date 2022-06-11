
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

-   Generalized linear mixed models

-   Shared component models

-   Spatial frailty models

## Installing LSB Core

You may need to install the LSB Core

lsb\_release is part of a software package called the LSB core, which is
not necessarily installed on your system by default. To install it, run
the command below that corresponds to your specific system:

### Ubuntu, Debian

sudo apt-get update && sudo apt-get install lsb-core

### CentOS

sudo yum update && sudo yum install redhat-lsb-core

### Fedora

sudo dnf update && sudo dnf install redhat-lsb-core

### OpenSUSE

sudo zypper update && sudo zypper install lsb-core

### Arch

pacman -Syu lsb-release

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
#>                  mean    median      mode       sd     lower    upper
#> (Intercept)  0.290180  0.293232  0.297064 0.096921  0.084517 0.464152
#> X1          -0.112440 -0.111390 -0.110643 0.083823 -0.270183 0.053553
#> X2           0.422433  0.422807  0.416497 0.095397  0.221108 0.604123
sglmm_mod$unrestricted$summary_fixed
#>                  mean    median      mode       sd     lower    upper
#> (Intercept)  0.193310  0.192793  0.195883 0.112536 -0.019229 0.408931
#> X1          -0.105333 -0.108120 -0.116720 0.097176 -0.299285 0.075365
#> X2           0.554684  0.524568  0.511294 0.254681  0.046469 1.043592
rglmm_mod$unrestricted$summary_fixed
#>                  mean    median      mode       sd     lower    upper
#> (Intercept)  0.192251  0.189350  0.183941 0.110390 -0.024940 0.390233
#> X1          -0.104634 -0.107476 -0.105396 0.096672 -0.302749 0.074206
#> X2           0.553324  0.518179  0.485339 0.268085  0.053777 1.120851
rglmm_mod$restricted$summary_fixed
#>                  mean    median      mode       sd     lower    upper
#> (Intercept)  0.192251  0.189350  0.183941 0.110390 -0.024940 0.390233
#> X1          -0.120125 -0.118708 -0.119207 0.086898 -0.287106 0.054653
#> X2           0.428399  0.427712  0.426570 0.101527  0.214849 0.608716

sglmm_mod$unrestricted$summary_hyperpar
#>                       mean   median     mode       sd   lower    upper
#> Precision for reg 93.31269 2.227233 1.487469 464.0461 0.40322 415.3316
```

The *SVIF* function provides the coefficients variance ratio between two
models. We can notice that the restricted model alleviates the variance
inflation

``` r
SVIF(rglmm_mod$restricted, rglmm_mod$unrestricted)
#>     parameter      VIF
#> 1 (Intercept) 1.000000
#> 2          X1 1.237604
#> 3          X2 6.972394
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
#>             mean    median      mode       sd     lower     upper
#> alpha1  0.468548  0.470915  0.470332 0.066425  0.342601  0.599288
#> alpha2  0.093245  0.093468  0.096582 0.054381 -0.014831  0.197089
#> X11_1  -0.458648 -0.455249 -0.455418 0.107509 -0.670253 -0.251454
#> X12_1  -0.556626 -0.565192 -0.556368 0.727899 -2.011461  0.815949
#> X21_2  -0.908133 -0.908408 -0.910508 0.057989 -1.015133 -0.790359
#> X12_2  -0.385854 -0.374461 -0.361140 0.315915 -0.969523  0.239076
rscm_mod$summary_fixed
#>             mean    median      mode       sd     lower     upper
#> alpha1  0.514514  0.514852  0.518343 0.062271  0.391794  0.628163
#> alpha2  0.107683  0.108888  0.108218 0.052797  0.002684  0.205958
#> X11_1  -0.556916 -0.554289 -0.549720 0.111126 -0.769681 -0.333068
#> X12_1  -0.191676 -0.196055 -0.199189 0.150284 -0.487269  0.095290
#> X21_2  -0.945358 -0.946199 -0.946929 0.052217 -1.039213 -0.834244
#> X12_2  -0.099804 -0.098899 -0.095106 0.069917 -0.228329  0.037779

scm_mod$summary_hyperpar
#>                           mean    median      mode        sd    lower     upper
#> Precision for psi     2.297820  2.213172  2.137552  0.612074 1.233218  3.534963
#> Precision for phi1   13.357326  9.728952  7.277940 12.234102 0.807289 36.489783
#> Precision for phi2   17.315173 14.488147 12.275385 10.875393 3.090969 38.617265
#> Beta for psi_gamma    2.592330  2.574339  2.561000  0.282854 2.060594  3.161575
#> Delta                 1.602811  1.596515  1.591446  0.087321 1.434889  1.768535
#> Precision for psi SC  0.878675  0.866814  0.852711  0.174781 0.585320  1.266731
rscm_mod$summary_hyperpar
#>                           mean    median      mode        sd    lower     upper
#> Precision for psi     0.802387  0.775819  0.750927  0.201922 0.445629  1.208257
#> Precision for phi1   13.995620  9.359445  6.070808 14.971263 0.257390 42.276127
#> Precision for phi2   22.518171 16.917707 12.112072 19.809224 0.855579 60.670636
#> Beta for psi_gamma    2.607235  2.591025  2.577442  0.254488 2.127917  3.118987
#> Delta                 1.613391  1.612117  1.606086  0.077287 1.459759  1.756817
#> Precision for psi SC  0.305762  0.301283  0.293250  0.060613 0.199598  0.427545
```

``` r
SVIF(rscm_mod, scm_mod)
#>   parameter        VIF
#> 1    alpha1  1.1378669
#> 2    alpha2  1.0609035
#> 3     X11_1  0.9359621
#> 4     X12_1 23.4593921
#> 5     X21_2  1.2332962
#> 6     X12_2 20.4162006
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
#>                  mean    median      mode       sd     lower     upper
#> (Intercept) -0.458939 -0.456202 -0.455831 0.077948 -0.609153 -0.308553
#> X1           0.366999  0.366651  0.368283 0.076112  0.217046  0.505295
#> X2          -0.652859 -0.656143 -0.660494 0.089487 -0.812906 -0.464613
rsfm_mod$unrestricted$summary_fixed
#>                  mean    median      mode       sd     lower     upper
#> (Intercept) -0.435600 -0.435694 -0.441499 0.076557 -0.581959 -0.288991
#> X1           0.369562  0.366782  0.366116 0.078220  0.221902  0.519001
#> X2          -0.641793 -0.642224 -0.649456 0.230732 -1.069026 -0.160841
rsfm_mod$restricted$summary_fixed
#>                  mean    median      mode       sd     lower     upper
#> (Intercept) -0.438975 -0.438133 -0.441816 0.075317 -0.573785 -0.280953
#> X1           0.374667  0.372933  0.375549 0.077295  0.234508  0.527128
#> X2          -0.691750 -0.689692 -0.687067 0.084866 -0.865797 -0.539728
```

``` r
SVIF(weibull_mod$unrestricted, rsfm_mod$unrestricted)
#>     parameter      VIF
#> 1 (Intercept) 0.964628
#> 2          X1 1.056159
#> 3          X2 6.648073
SVIF(weibull_mod$unrestricted, rsfm_mod$restricted)
#>     parameter       VIF
#> 1 (Intercept) 0.9336327
#> 2          X1 1.0313273
#> 3          X2 0.8993890
```
