#' @title Restricted Spatial Frailty Model in INLA
#'
#' @description Fit a Restricted Spatial Frailty model using INLA
#'
#' @param f INLA formula ?inla.surv
#' @param data data.frame containing, at least, \code{time}, \code{status}, \code{covariates}, \code{area} list
#' @param family 'exponential', 'weibull', 'weibullcure', 'loglogistic', 'gamma', 'lognormal' or 'pwe'
#' @param W Adjacency matrix
#' @param proj 'none', 'rhz' or 'spock'
#' @param fast Should we use the reduction operator?
#' @param nsamp Sample size to use the projection approach
#' @param ... Other parameters used in ?inla
#'
#' @return INLA object with corrected parameters
#'
#' @import INLA
#' @importFrom stats density
#'
#' @export

rsfm_inla <- function(f, data, family, W = NULL,
                      proj = "none", fast = TRUE, nsamp = 1000, ...) {
  ##-- Time
  time_start <- Sys.time()

  family <- switch(family,
                   "exponential" = "exponential.surv",
                   "weibull" = "weibull.surv",
                   "weibullcure" = "weibullcure",
                   "loglogistic" = "loglogistic.surv",
                   "gamma" = "gamma.surv",
                   "lognormal" = "lognormal.surv",
                   "pwe" = "coxph",
                   stop(family, "family is not implemented."))

  ##-- Updating formula
  inla_formula <- update_inla_formula(formula = f)

  f <- inla_formula$formula
  reg_name <- inla_formula$var_restricted
  reg_size <- inla_formula$size_restricted

  ##-- Model
  time_start_inla <- Sys.time()
  mod <- inla(formula = f, data = data, family = family,
              control.compute = list(config = TRUE), ...)

  model_sample <- inla.posterior.sample(result = mod, n = nsamp, use.improved.mean = TRUE)
  hyperpar_samp <- inla.hyperpar.sample(result = mod, n = nsamp, improve.marginals = TRUE)
  time_end_inla <- Sys.time()

  X <- as.matrix(mod$model.matrix)
  fixed_vars <- mod$names.fixed
  id_fixed <- paste0(fixed_vars, ":", 1)

  if(length(reg_name) > 0) {
    reg_pos <- data[[reg_name]]
    id_latent <- paste0(reg_name, ":", 1:reg_size)
  } else {
    id_latent <- ""
  }

  beta_samp <- do.call(args = lapply(model_sample, select_marginal, ids = id_fixed), what = "rbind")
  W_samp <- do.call(args = lapply(model_sample, select_marginal, ids = id_latent), what = "rbind")
  hyperpar_samp <- hyperpar_samp

  if(length(reg_name) > 0){
    colnames(W_samp) <- paste0("S", 1:reg_size)
  } else {
    W_samp <- NULL
  }

  colnames(beta_samp) <- fixed_vars

  ##-- Correcting the model
  if(proj != "none" & length(reg_name) > 0) {
    ##-- Projection matrices
    time_start_correction <- Sys.time()
    if(fast){
      proj_aux <- proj_mat(X = X, groups = reg_pos, method = proj)
    } else{
      proj_aux <- proj_mat(X = X, groups = NULL, method = proj)
    }

    ##-- Fixed effects
    if(fast) {
      W_aux <- t(proj_aux$Paux%*%t(W_samp))
    } else{
      W_aux <- t(proj_aux$Paux%*%t(W_samp[, reg_pos]))
    }

    beta_ast <- matrix(beta_samp + W_aux, ncol = ncol(X))
    colnames(beta_ast) <- fixed_vars

    hyperpar_ast <- hyperpar_samp

    ##-- Random effects
    if(fast) {
      W_ast <- proj_aux$Px_ort%*%t(W_samp)
      Z_ast <- W_samp[, reg_pos] - t(W_ast[reg_pos, ])
    } else{
      W_ast <- proj_aux$Px_ort%*%t(W_samp[, reg_pos])
      Z_ast <- W_samp[, reg_pos] - t(W_ast)

      W_ast <- apply(X = W_ast, MARGIN = 2, meang, g = reg_pos)
    }

    time_end_correction <- Sys.time()
  } else {
    beta_ast <- NULL
    hyperpar_ast <- NULL
    W_ast <- NULL
    Z_ast <- NULL

    time_start_correction <- time_end_correction <- Sys.time()
  }

  sample <- cbind(hyperpar_samp, beta_samp, W_samp)
  sample <- as.data.frame(sample)

  names(sample) <- c(colnames(hyperpar_samp), fixed_vars, colnames(W_samp))

  if(proj != "none" & length(reg_name) > 0) {
    sample_ast <- cbind(hyperpar_samp, beta_ast, t(W_ast))
    sample_ast <- as.data.frame(sample_ast)

    names(sample_ast) <- c(colnames(hyperpar_samp), fixed_vars, colnames(W_samp))
  } else {
    sample_ast <- NULL
  }

  ##-- Time
  time_tab <- data.frame(INLA = as.numeric(difftime(time1 = time_end_inla, time2 = time_start_inla, units = "secs")),
                         correction = as.numeric(difftime(time1 = time_end_correction, time2 = time_start_correction, units = "secs")),
                         total = as.numeric(difftime(time1 = Sys.time(), time2 = time_start, units = "secs")))

  ##-- Return
  out <- list()

  out$unrestricted <- list()
  out$unrestricted$sample <- sample
  out$unrestricted$summary_fixed <- chain_summary(obj = beta_samp)
  out$unrestricted$summary_hyperpar <- chain_summary(obj = hyperpar_samp)
  out$unrestricted$summary_random <- chain_summary(obj = W_samp)
  out$unrestricted$out <- mod

  out$restricted <- list()
  out$restricted$sample <- sample_ast
  out$restricted$summary_fixed <- chain_summary(obj = beta_ast)
  out$restricted$summary_hyperpar <- chain_summary(obj = hyperpar_ast)
  out$restricted$summary_random <- chain_summary(obj = W_ast)

  out$time <- time_tab

  return(out)
}
