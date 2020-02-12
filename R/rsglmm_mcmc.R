#' @title Restricted Spatial Generalized Linear Mixed model
#'
#' @description Fit a Restricted Spatial Generalized Linear Mixed model using ngspatial
#'
#' @param data data.frame containing, at least, \code{time}, \code{status}, \code{covariates}, \code{area} list
#' @param formula INLA formula ?inla.surv
#' @param family 'exponential', 'weibull', 'weibullcure', 'loglogistic', 'gamma', 'lognormal' or 'pwe'
#' @param W Adjacency matrix
#' @param area Areal variable name in data
#' @param proj 'none', 'rhz' or 'spock'
#' @param nsamp Sample size to use the projection approach
#' @param burnin Burnin period
#' @param thin Lag parameter
#' @param attractive The number of attractive Moran eigenvectors to use
#' @param ... Other parameters used in ?ngspatial::sparse.sglmm
#'
#' @return INLA object with corrected parameters
#'
#' @import INLA
#'
#' @export

rsglmm_mcmc <- function(data, formula, family,
                        W, area,
                        proj, nsamp = 1000, burnin = 5000, thin = 1,
                        attractive = round(0.5*(nrow(W)/2)),
                        ...) {
  ##-- Time
  time_start <- Sys.time()

  ##-- Model
  time_start_mcmc <- Sys.time()
  dimnames(W)[[1]] <- NULL
  x <- TRUE

  if(!is.null(data$offset)) {
    mod <- ngspatial::sparse.sglmm(data = data, formula = formula,
                                   family = family, offset = offset,
                                   method = "RSR", A = W,
                                   minit = burnin, maxit = burnin + nsamp*thin,
                                   attractive = attractive, x = x,
                                   ...)
  } else {
    mod <- ngspatial::sparse.sglmm(data = data, formula = formula,
                                   family = family,
                                   method = "RSR", A = W,
                                   minit = burnin, maxit = burnin + nsamp*thin,
                                   attractive = attractive, x = x,
                                   ...)
  }

  time_end_mcmc <- Sys.time()

  ##-- Samples
  pos_samp <- seq(1, (nsamp*thin), by = thin)

  tau_s <- mod$tau.s.sample[pos_samp]
  if(family == "gaussian") {
    tau_g <- mod$tau.h.sample[pos_samp]

    hyperpar_ast <- cbind.data.frame(tau_g, tau_s)
    names(hyperpar_ast) <- c("Precision for the Gaussian observations", sprintf("Precision for %s", area))
  }
  if(family %in% c("poisson", "binomial")) {
    tau_g <- mod$tau.h.sample

    hyperpar_ast <- cbind.data.frame(tau_s)
    names(hyperpar_ast) <- sprintf("Precision for %s", area)
  }

  W_ast <- data.frame(mod$gamma.sample[pos_samp,]%*%t(mod$M))
  names(W_ast) <- paste(area, 1:ncol(W_ast), sep = "_")

  beta_ast <- data.frame(mod$beta.sample[pos_samp, , drop = FALSE])
  names(beta_ast) <- names(mod$coefficients)

  sample_ast <- cbind.data.frame(hyperpar_ast, beta_ast, W_ast)

  ##-- Time
  time_tab <- data.frame(MCMC = as.numeric(difftime(time1 = time_end_mcmc, time2 = time_start_mcmc, units = "secs")),
                         total = as.numeric(difftime(time1 = Sys.time(), time2 = time_start, units = "secs")))

  ##-- Return
  out <- list()

  out$unrestricted <- list()

  out$restricted <- list()
  out$restricted$sample <- sample_ast
  out$restricted$summary_fixed <- chain_summary(obj = beta_ast)
  out$restricted$summary_hyperpar <- chain_summary(obj = hyperpar_ast)
  out$restricted$summary_random <- chain_summary(obj = W_ast)

  out$time <- time_tab

  return(out)
}
