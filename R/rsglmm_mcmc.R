#' @title Restricted Spatial Generalized Linear Mixed model
#'
#' @description Fit a Restricted Spatial Generalized Linear Mixed model using ngspatial
#'
#' @param data an data frame or list containing the variables in the model.
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param family allowed families are: 'gaussian', 'poisson' and 'binomial'.
#' @param E known component, in the mean for the Poisson likelihoods defined as E = exp(\eqn{\eta}), where \eqn{\eta} is the linear predictor. If not provided it is set to 1.
#' @param n a vector containing the number of trials for the binomial likelihood and variantes, or the number of required successes for the nbinomial2 likelihood. Default value is set to 1.
#' @param W adjacency matrix.
#' @param area areal variable name in \code{data}.
#' @param proj 'hh'
#' @param nsamp number of desired. samples Default = 1000.
#' @param burnin burnin size.
#' @param lag lag parameter.
#' @param attractive the number of attractive Moran eigenvectors to use. See ?ngspatial::sparse.sglmm for more information.
#' @param ... other parameters used in ?ngspatial::sparse.sglmm
#'
#' @return \item{$unrestricted}{A list containing
#'                                \itemize{
#'                                   \item $sample a sample of size nsamp for all parameters in the model
#'                                   \item $summary_fixed summary measures for the coefficients
#'                                   \item $summary_hyperpar summary measures for hyperparameters
#'                                   \item $summary_random summary measures for random quantities
#'                                 }
#'                              }
#' \item{$restricted}{A list containing
#'                                \itemize{
#'                                   \item $sample a sample of size nsamp for all parameters in the model
#'                                   \item $summary_fixed summary measures for the coefficients
#'                                   \item $summary_hyperpar summary measures for hyperparameters
#'                                   \item $summary_random summary measures for random quantities
#'                                 }
#'                              }
#'
#' \item{$out}{ngspatial output}
#' \item{$time}{time elapsed for fitting the model}
#'
#' @importFrom ngspatial sparse.sglmm
#' @importFrom stats offset
#'
#' @export

rsglmm_mcmc <- function(data, formula, family, E, n,
                        W, area,
                        proj, nsamp, burnin, lag,
                        attractive = round(0.5*(nrow(W)/2)),
                        ...) {
  ##-- Time
  time_start <- Sys.time()

  if(!is.null(E)) data$offset <- log(E)
  if(!is.null(n)) data$offset <- log(n)

  ##-- Model
  time_start_mcmc <- Sys.time()
  dimnames(W)[[1]] <- NULL
  x <- TRUE

  if(suppressWarnings(!is.null(data$offset))) {
    mod <- ngspatial::sparse.sglmm(data = data, formula = formula,
                                   family = family, offset = offset,
                                   method = "RSR", A = W,
                                   attractive = attractive, x = x,
                                   minit = burnin, maxit = burnin + lag*nsamp,
                                   ...)
  } else {
    mod <- ngspatial::sparse.sglmm(data = data, formula = formula,
                                   family = family,
                                   method = "RSR", A = W,
                                   attractive = attractive, x = x,
                                   minit = burnin, maxit = burnin + lag*nsamp,
                                   ...)
  }

  time_end_mcmc <- Sys.time()

  ##-- Samples
  pos_samp <- seq(burnin + lag,  burnin + lag*nsamp, by = lag)

  tau_s <- mod$tau.s.sample[pos_samp]
  if(family == "gaussian") {
    tau_g <- mod$tau.h.sample[pos_samp]

    hyperpar_ast <- cbind.data.frame(tau_g, tau_s)
    names(hyperpar_ast) <- c("Precision for the Gaussian observations", sprintf("Precision for %s", area))
  } else {
    hyperpar_ast <- cbind.data.frame(tau_s)
    names(hyperpar_ast) <- sprintf("Precision for %s", area)
  }

  W_ast <- data.frame(mod$gamma.sample[pos_samp, ]%*%t(mod$M))
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
  out$restricted$summary_fixed <- chain_summary(sample = beta_ast)
  out$restricted$summary_hyperpar <- chain_summary(sample = hyperpar_ast)
  out$restricted$summary_random <- chain_summary(sample = W_ast)

  out$out <- mod
  out$time <- time_tab

  return(out)
}
