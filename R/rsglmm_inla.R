#' @title Restricted Spatial Generalized Linear Mixed model in INLA
#'
#' @description Fit a Restricted Spatial Generalized Linear Mixed model using INLA
#'
#' @param data an data frame or list containing the variables in the model.
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param family some allowed families are: 'gaussian', 'poisson' and 'binomial'. The family availability will depend on the approach.
#' @param E known component, in the mean for the Poisson likelihoods defined as E = exp(\eqn{\eta}), where \eqn{\eta} is the linear predictor. If not provided it is set to 1.
#' @param n a vector containing the number of trials for the binomial likelihood and variantes, or the number of required successes for the nbinomial2 likelihood. Default value is set to 1..
#' @param W adjacency matrix.
#' @param proj 'none', 'rhz' or 'spock'
#' @param nsamp number of desired. samples Default = 1000.
#' @param ... other parameters used in ?INLA::inla
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
#' \item{$out}{INLA output}
#' \item{$time}{time elapsed for fitting the model}
#'
#' @importFrom INLA inla inla.posterior.sample inla.hyperpar.sample
#'
#' @export

rsglmm_inla <- function(data, formula, family,
                        E = NULL, n = NULL,
                        W = NULL,
                        proj = "none", nsamp = 1000,
                        ...) {
  ##-- Time
  time_start <- Sys.time()

  ##-- Updating formula
  inla_formula <- update_inla_formula(formula = formula)

  formula <- inla_formula$formula

  reg_name_r <- inla_formula$var_restricted
  reg_size_r <- inla_formula$size_restricted

  reg_name_u <- inla_formula$vars_unrestricted
  reg_size_u <- inla_formula$size_unrestricted

  ##-- Model
  time_start_inla <- Sys.time()

  E_offset_inla <- E
  n_offset_inla <- n

  args <- list(...)
  args$control.compute$config <- TRUE

  inla_aux <- function(...) INLA::inla(formula = formula, data = data, family = family,
                                       E = E_offset_inla, Ntrials = n_offset_inla, ...)

  mod <- do.call(what = inla_aux, args = args)

  # mod <- inla(formula = formula, data = data, family = family,
  #             E = E_offset_inla, Ntrials = n_offset_inla, ...)

  model_sample <- INLA::inla.posterior.sample(result = mod, n = nsamp, use.improved.mean = TRUE)
  hyperpar_samp <- INLA::inla.hyperpar.sample(result = mod, n = nsamp, improve.marginals = TRUE)
  time_end_inla <- Sys.time()

  X <- as.matrix(mod$model.matrix)
  fixed_vars <- mod$names.fixed
  id_fixed <- paste0(fixed_vars, ":", 1)

  id_latent_r <- id_latent_u <- ""

  if(length(reg_name_r) > 0) {
    reg_pos_r <- data[[reg_name_r]]
    id_latent_r <- paste0(reg_name_r, ":", 1:reg_size_r)
  }

  if(length(reg_name_u) > 0) {
    reg_pos_u <- data[[reg_name_u]]
    id_latent_u <- paste0(reg_name_u, ":", 1:reg_size_u)
  }

  beta_samp <- do.call(args = lapply(model_sample, select_marginal, ids = id_fixed), what = "rbind")
  W_samp_r <- do.call(args = lapply(model_sample, select_marginal, ids = id_latent_r), what = "rbind")
  W_samp_u <- do.call(args = lapply(model_sample, select_marginal, ids = id_latent_u), what = "rbind")
  hyperpar_samp <- hyperpar_samp

  if(length(reg_name_r) > 0){
    colnames(W_samp_r) <- paste(reg_name_r, 1:reg_size_r, sep = "_")
  } else {
    W_samp_r <- NULL
  }

  if(length(reg_name_u) > 0){
    colnames(W_samp_u) <- paste(reg_name_u, 1:reg_size_u, sep = "_")
  } else {
    W_samp_u <- NULL
  }

  colnames(beta_samp) <- fixed_vars

  ##-- Correcting the model
  if(proj == "rhz" & length(reg_name_r) > 0) {
    ##-- Projection matrices
    time_start_correction <- Sys.time()
    proj_aux <- proj_mat(X = X, groups = reg_pos_r, method = proj)

    ##-- Fixed effects
    W_aux <- t(proj_aux$Paux%*%t(W_samp_r))

    beta_ast <- matrix(beta_samp + W_aux, ncol = ncol(X))
    colnames(beta_ast) <- fixed_vars

    hyperpar_ast <- hyperpar_samp

    ##-- Random effects
    W_ast <- proj_aux$Px_ort%*%t(W_samp_r)

    if(reg_size_r < nrow(X)) {
      Z_ast <- W_samp_r[, reg_pos_r] - t(W_ast[reg_pos_r, ])
    } else {
      Z_ast <- NULL
    }

    W_ast <- t(W_ast)
    W_ast_u <- W_samp_u

    colnames(W_ast) <- colnames(W_samp_r)
    if(reg_size_r < nrow(X)) colnames(Z_ast) <- paste("Z", 1:ncol(Z_ast), sep = "_")

    time_end_correction <- Sys.time()
  } else {
    if(proj == "spock" & length(reg_name_r) > 0) {
      beta_ast <- beta_samp
      hyperpar_ast <- hyperpar_samp
      W_ast <- W_samp_r
      W_ast_u <- W_samp_u
      Z_ast <- NULL
    } else {
      beta_ast <- NULL
      hyperpar_ast <- NULL
      W_ast <- NULL
      W_ast_u <- NULL
      Z_ast <- NULL
    }

    time_start_correction <- time_end_correction <- Sys.time()
  }

  if(proj != "spock") {
    sample <- cbind(hyperpar_samp, beta_samp, W_samp_r, W_samp_u)
    sample <- as.data.frame(sample)

    names(sample) <- c(colnames(hyperpar_samp), fixed_vars, colnames(W_samp_r), colnames(W_samp_u))
  } else {
    sample <- NULL
  }

  if(proj == "rhz" & length(reg_name_r) > 0) {
    sample_ast <- cbind(hyperpar_samp, beta_ast, W_ast, Z_ast)
    if(length(reg_name_u) > 0) sample_ast <- cbind(sample_ast, W_samp_u)

    sample_ast <- as.data.frame(sample_ast)

    if(length(reg_name_u) > 0) names(sample_ast) <- c(colnames(hyperpar_samp),
                                                      fixed_vars,
                                                      colnames(W_samp_r),
                                                      colnames(W_samp_u),
                                                      colnames(Z_ast))
    if(length(reg_name_u) == 0) names(sample_ast) <- c(colnames(hyperpar_samp),
                                                       fixed_vars,
                                                       colnames(W_samp_r),
                                                       colnames(Z_ast))
  } else {
    if(proj == "spock"){
      sample_ast <- cbind(hyperpar_samp, beta_samp, W_samp_r, W_samp_u)
      sample_ast <- as.data.frame(sample_ast)

      names(sample_ast) <- c(colnames(hyperpar_samp), fixed_vars, colnames(W_samp_r), colnames(W_samp_u))

      beta_samp <- NULL
      hyperpar_samp <- NULL
      W_samp_r <- NULL
      W_samp_u <- NULL
      Z_ast <- NULL
    } else {
      sample_ast <- NULL
    }
  }

  ##-- Time
  time_tab <- data.frame(INLA = as.numeric(difftime(time1 = time_end_inla, time2 = time_start_inla, units = "secs")),
                         correction = as.numeric(difftime(time1 = time_end_correction, time2 = time_start_correction, units = "secs")),
                         total = as.numeric(difftime(time1 = Sys.time(), time2 = time_start, units = "secs")))

  ##-- Return
  out <- list()

  out$unrestricted <- list()
  out$unrestricted$sample <- sample
  out$unrestricted$summary_fixed <- chain_summary(sample = beta_samp)
  out$unrestricted$summary_hyperpar <- chain_summary(sample = hyperpar_samp)
  out$unrestricted$summary_random <- chain_summary(sample = cbind(W_samp_r, W_samp_u))

  out$restricted <- list()
  out$restricted$sample <- sample_ast
  out$restricted$summary_fixed <- chain_summary(sample = beta_ast)
  out$restricted$summary_hyperpar <- chain_summary(sample = hyperpar_ast)
  out$restricted$summary_random <- rbind(chain_summary(sample = cbind(W_ast, W_ast_u)),
                                         chain_summary(sample =  Z_ast))

  out$out <- mod
  out$time <- time_tab

  return(out)
}
