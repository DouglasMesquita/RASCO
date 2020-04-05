#' @title Restricted Spatial Frailty Model in BUGS
#'
#' @description Fit a Restricted Spatial Frailty model using BUGS
#'
#' @usage rsfm_mcmc(model, data, inits, parameters, covariates, area,
#'                  proj = "none", fast = TRUE, nsamp = 1000, burnin = 5000, lag = 10, ...)
#'
#' @param model character with the model in text format.
#' @param data an data frame or list containing the variables in the model.
#' @param inits initial values for the chain.
#' @param parameters vector of parameters names to restore.
#' @param covariates vector of covariates names.
#' @param area areal variable name in \code{data}.
#' @param proj 'none' or 'rhz'.
#' @param fast TRUE to use the reduction operator.
#' @param nsamp number of desired. samples Default = 1000.
#' @param burnin burnin size.
#' @param lag lag parameter.
#' @param ... other parameters used in ?R2OpenBUGS::bugs
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
#' \item{$out}{BUGS output}
#' \item{$time}{time elapsed for fitting the model}
#'
#' @importFrom R2OpenBUGS bugs

rsfm_mcmc <- function(model, data, inits, parameters, covariates, area,
                      proj = "none", fast = TRUE, nsamp = 1000, burnin = 5000, lag = 10, ...){
  if(is.na(findOpenBUGS())) stop("We didn't find OpenBugs on your computer.")

  ##-- Time
  time_start <- Sys.time()

  ##-- Model
  time_start_mcmc <- Sys.time()
  mod <- R2OpenBUGS::bugs(data = data, inits = inits,
                          model.file = model,
                          parameters.to.save = parameters,
                          n.chains = 1, n.iter = burnin + nsamp*lag, n.burnin = burnin, n.thin = lag)
  time_end_mcmc <- Sys.time()

  sample <- as.data.frame(mod$sims.list)
  sample <- sample[, -ncol(sample)]

  if(!is.null(area)) {
    names(sample)[1:(3+length(covariates))] <- c("alpha parameter for weibullsurv", sprintf("Precision for %s", area), "(Intercept)", covariates)
    hyperpar_samp <- sample[, c("alpha parameter for weibullsurv", sprintf("Precision for %s", area))]
    W_samp <- sample[, -which(c("alpha parameter for weibullsurv", sprintf("Precision for %s", area), "(Intercept)", covariates) %in% names(sample))]
  } else{
    names(sample) <- c("alpha parameter for weibullsurv", "(Intercept)", covariates)
    hyperpar_samp <- sample[, "alpha parameter for weibullsurv", drop = FALSE]
    W_samp <- sample[, -which(c("alpha parameter for weibullsurv", "(Intercept)", covariates) %in% names(sample))]
  }

  beta_samp <- sample[, c("(Intercept)", covariates), drop = FALSE]

  if(ncol(W_samp) == 0) W_samp <- NULL

  ##-- Correcting the model
  if(proj != "none" & !is.null(area)){
    X <- as.matrix(data.frame(`(Intercept)` = 1, data[covariates], check.names = FALSE))
    reg_pos <- data[[area]]

    ##-- Projection matrices
    time_start_correction <- Sys.time()
    if(fast){
      proj_aux <- proj_mat(X = X, groups = data[["area"]], method = proj)
    } else{
      proj_aux <- proj_mat(X = X, groups = NULL, method = proj)
    }

    ##-- Fixed effects
    if(fast){
      W_aux <- t(proj_aux$Paux%*%t(W_samp))
    } else{
      W_aux <- t(proj_aux$Paux%*%t(W_samp[, reg_pos]))
    }

    beta_ast <- beta_samp + W_aux
    hyperpar_ast <- hyperpar_samp

    ##-- Random effects
    if(fast){
      W_ast <- proj_aux$Px_ort%*%t(W_samp)
      Z_ast <- W_samp[, reg_pos] - t(W_ast[reg_pos, ])
    } else{
      W_ast <- proj_aux$Px_ort%*%t(W_samp[, reg_pos])
      Z_ast <- W_samp[, reg_pos] - t(W_ast)

      W_ast <- apply(X = W_ast, MARGIN = 2, meang, g = reg_pos)
    }

    ##-- Time
    time_end_correction <- Sys.time()
  } else {
    beta_ast <- NULL
    hyperpar_ast <- NULL
    W_ast <- NULL
    Z_ast <- NULL

    time_start_correction <- time_end_correction <- Sys.time()
  }

  if(proj != "none" & !is.null(area)) {
    sample_ast <- cbind(hyperpar_samp, beta_ast, t(W_ast))
    colnames(sample_ast) <- c(colnames(hyperpar_samp), covariates, colnames(W_samp))
  } else {
    sample_ast <- NULL
  }

  time_tab <- data.frame(MCMC = as.numeric(difftime(time1 = time_end_mcmc, time2 = time_start_mcmc, units = "secs")),
                         correction = as.numeric(difftime(time1 = time_end_correction, time2 = time_start_correction, units = "secs")),
                         total = as.numeric(difftime(time1 = Sys.time(), time2 = time_start, units = "secs")))

  ##-- Return
  out <- list()

  out$unrestricted <- list()
  out$unrestricted$sample <- sample
  out$unrestricted$summary_fixed <- chain_summary(sample = beta_samp)
  out$unrestricted$summary_hyperpar <- chain_summary(sample = hyperpar_samp)
  out$unrestricted$summary_random <- chain_summary(sample = W_samp)

  out$restricted <- list()
  out$restricted$sample <- sample_ast
  out$restricted$summary_fixed <- chain_summary(sample = beta_ast)
  out$restricted$summary_hyperpar <- chain_summary(sample = hyperpar_ast)
  out$restricted$summary_random <- chain_summary(sample = W_ast)

  out$out <- mod
  out$time <- time_tab

  return(out)
}
