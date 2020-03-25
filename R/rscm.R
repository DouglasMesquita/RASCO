#' @title Restricted Shared Component model
#'
#' @description Fit a Restricted Shared Component model for two diseases
#'
#' @usage rscm(data, formula1, formula2, family = c("poisson", "poisson"),
#'             E1 = NULL, E2 = NULL, area = NULL, neigh = NULL,
#'             proj = "none", nsamp = 1000,
#'             priors = list(prior_gamma = c(0, 0.1),
#'                           prior_prec = list(tau_s = c(0.5, 0.05),
#'                                             tau_1 = c(0.5, 0.05),
#'                                             tau_2 = c(0.5, 0.05))),
#'             random_effects = list(shared = TRUE,
#'                                   specific_1 = TRUE,
#'                                   specific_2 = TRUE),
#'             ...)
#'
#' @param data an data frame or list containing the variables in the model.
#' @param formula1 an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted for disease 1.
#' @param formula2 an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted for disease 2.
#' @param family a vector of size two with two families. Some allowed families are: poisson, nbinomial, zeroinflatedpoisson0, zeroinflatednbinomial0. See INLA::inla.list.models().
#' @param E1 known component, for disease 1, in the mean for the Poisson likelihoods defined as E = exp(\eqn{\eta}), where \eqn{\eta} is the linear predictor. If not provided it is set to 1.
#' @param E2 known component, for disease 2, in the mean for the Poisson likelihoods defined as E = exp(\eqn{\eta}), where \eqn{\eta} is the linear predictor. If not provided it is set to 1.
#' @param area areal variable name in \code{data}.
#' @param neigh neighborhood structure. A \code{SpatialPolygonsDataFrame} object.
#' @param proj 'none' or 'spock'.
#' @param nsamp number of desired. samples Default = 1000.
#' @param priors a list containing:
#'     \itemize{
#'        \item prior_gamma: a vector of size two containing mean and precision for the normal distribution applied for \eqn{\gamma}
#'        \item prior_prec: a list with:
#'        \itemize{
#'            \item tau_s: a vector of size two containing shape and scale for the gamma distribution applied for \eqn{\tau_s}
#'            \item tau_1: a vector of size two containing shape and scale for the gamma distribution applied for \eqn{\tau_1}
#'            \item tau_2: a vector of size two containing shape and scale for the gamma distribution applied for \eqn{\tau_2}
#'        }
#'     }
#' @param random_effects a list determining which effects should we include in the model. Default: list(shared = TRUE, specific_1 = TRUE, specific_2 = TRUE).
#' @param ... other parameters used in ?INLA::inla
#'
#' @details The fitted model is given by
#' \deqn{Y_1 ~ Poisson(E_1\theta_1),}
#' \deqn{Y_2 ~ Poisson(E_2\theta_2),}
#'
#' \deqn{log(\theta_1) = X\beta + \gamma\psi + \phi_1,}
#' \deqn{log(\theta_2) = X\beta + \psi + \phi_2,}
#'
#' \deqn{\psi ~ ICAR(\tau_s); \phi_1 ~ ICAR(\tau_1); \phi_2 ~ ICAR(\tau_2).}
#'
#' \deqn{\delta = \sqrt\gamma}
#'
#' @examples
#' library(spdep)
#'
#' set.seed(123456)
#'
#' ##-- Spatial structure
#' data("neigh_RJ")
#'
#' ##-- Parameters
#' alpha_1 <- 0.5
#' alpha_2 <- 0.1
#' beta_1 <- c(-0.5, -0.2)
#' beta_2 <- c(-0.8, -0.4)
#' tau_s <- 1
#' tau_1 <- tau_2 <- 10
#' delta <- 1.5
#'
#' ##-- Data
#' data <- rshared(alpha_1 = alpha_1, alpha_2 = alpha_2,
#'                 beta_1 = beta_1, beta_2 = beta_2,
#'                 delta = delta,
#'                 tau_1 = tau_1, tau_2 = tau_2, tau_s = tau_s,
#'                 confounding = "linear",
#'                 neigh = neigh_RJ)
#'
#' ##-- Models
#' scm_inla <- rscm(data = data,
#'                  formula1 = Y1 ~ X11 + X12,
#'                  formula2 = Y2 ~ X21 + X12,
#'                  family = c("nbinomial", "poisson"),
#'                  E1 = E1, E2 = E2,
#'                  area = "reg", neigh = neigh_RJ,
#'                  priors = list(prior_prec = list(tau_s = c(0.5, 0.05)), prior_gamma = c(0, 0.5)),
#'                  proj = "none", nsamp = 1000,
#'                  random_effects = list(shared = TRUE, specific_1 = TRUE, specific_2 = TRUE))
#'
#' rscm_inla <- rscm(data = data,
#'                   formula1 = Y1 ~ X11 + X12,
#'                   formula2 = Y2 ~ X21 + X12,
#'                   family = c("nbinomial", "poisson"),
#'                   E1 = E1, E2 = E2,
#'                   area = "reg", neigh = neigh_RJ,
#'                   priors = list(prior_prec = list(tau_s = c(0.5, 0.05)), prior_gamma = c(0, 0.5)),
#'                   proj = "spock", nsamp = 1000,
#'                   random_effects = list(shared = TRUE, specific_1 = TRUE, specific_2 = TRUE))
#'
#' ##-- Summary
#' scm_inla$summary_fixed
#' rscm_inla$summary_fixed
#'
#' scm_inla$summary_hyperpar
#' rscm_inla$summary_hyperpar
#'
#' @return \item{$sample}{a sample of size nsamp for all parameters in the model}
#' \item{$summary_fixed}{summary measures for the coefficients}
#' \item{$summary_hyperpar}{summary measures for hyperparameters}
#' \item{$summary_random}{summary measures for random quantities}
#' \item{$out}{INLA output}
#' \item{$time}{time elapsed for fitting the model}
#'
#' @importFrom stats as.formula terms model.frame update
#'
#' @export

rscm <- function(data, formula1, formula2, family = c("poisson", "poisson"),
                 E1 = NULL, E2 = NULL, area = NULL, neigh = NULL,
                 proj = "none", nsamp = 1000,
                 priors = list(prior_gamma = c(0, 0.1),
                               prior_prec = list(tau_s = c(0.5, 0.05),
                                                 tau_1 = c(0.5, 0.05),
                                                 tau_2 = c(0.5, 0.05))),
                 random_effects = list(shared = TRUE, specific_1 = TRUE, specific_2 = TRUE),
                 ...) {
  ##-- Time
  time_start <- Sys.time()

  ##-- Tests
  if(missing(formula1)) stop("You must provide the formula for disease 1")
  if(missing(formula2)) stop("You must provide the formula for disease 2")
  if(!proj %in% c("none", "spock")) stop("proj must be 'none' or 'spock'")

  random_effects <- append_list(list(shared = TRUE, specific_1 = TRUE, specific_2 = TRUE), random_effects)
  priors <- append_list(list(prior_gamma = c(0, 0.1),
                             prior_prec = list(tau_s = c(0.5, 0.05),
                                               tau_1 = c(0.5, 0.05),
                                               tau_2 = c(0.5, 0.05))),
                        priors)

  shared <- random_effects$shared
  specific_1 <- random_effects$specific_1
  specific_2 <- random_effects$specific_2

  prior_gamma <- priors$prior_gamma
  if(!is.list(priors$prior_prec) & length(priors$prior_prec) == 2) {
    message(sprintf("priors$prior_prec should be a list of tau_s, tau_1 and tau_2 entries.
                     Setting all priors as Gamma(%s, %s).", priors$prior_prec[1], priors$prior_prec[2]))
    prior_tau_s <- priors$prior_prec
    prior_tau_1 <- priors$prior_prec
    prior_tau_2 <- priors$prior_prec
  } else{
    prior_tau_s <- priors$prior_prec$tau_s
    prior_tau_1 <- priors$prior_prec$tau_1
    prior_tau_2 <- priors$prior_prec$tau_2
  }

  ##-- Setup INLA
  Y <- data[, c(deparse(formula1[[2]]), deparse(formula2[[2]]))]

  X1 <- all.vars(formula1[-2])
  X2 <- all.vars(formula2[-2])

  covs_s <- union(X1, X2)
  X1 <- data[, X1, drop = FALSE]
  X2 <- data[, X2, drop = FALSE]
  Xs <- data[, covs_s, drop = FALSE]

  n <- nrow(Y)
  n_covs1 <- ncol(X1)
  n_covs2 <- ncol(X2)
  n_covss <- ncol(Xs)

  if(!is.null(neigh)) {
    W <- nb2mat(neighbours = poly2nb(neigh), style = "B")
    Ws <- W1 <- W2 <- W
  }

  if(proj == "spock" & !is.null(area) & !is.null(neigh)) {
    time_start_correction <- Sys.time()

    if(specific_1) {
      neigh1 <- spock(X = as.matrix(cbind(1, X1)), map = neigh)
      W1 <- spdep::nb2mat(neigh1, style = "B", zero.policy = TRUE)
    }

    if(specific_2) {
      neigh2 <- spock(X = as.matrix(cbind(1, X2)), map = neigh)
      W2 <- spdep::nb2mat(neigh2, style = "B", zero.policy = TRUE)
    }

    if(shared) {
      neighs <- spock(X = as.matrix(cbind(1, Xs)), map = neigh)
      Ws <- spdep::nb2mat(neighs, style = "B", zero.policy = TRUE)
    }

    time_end_correction <- Sys.time()
  } else {
    time_end_correction <- time_start_correction <- Sys.time()
  }

  E1 <- substitute(E1)
  E2 <- substitute(E2)

  if(is.null(E1)) {
    E1 <- rep(1, n)
  } else {
    E1 <- data[, deparse(E1), drop = FALSE]
  }
  if(is.null(E2)) {
    E2 <- rep(1, n)
  } else {
    E2 <- data[, deparse(E2), drop = FALSE]
  }

  E <- data.frame(E1 = E1, E2 = E2)
  E <- c(E[, 1], E[, 2])

  inla_list <- list()
  inla_list$Y <- cbind(c(Y[, 1], rep(NA, n)),
                       c(rep(NA, n), Y[, 2]))

  inla_list$alpha1 <- rep(c(1, NA), each = n)          ##-- Intercept Y1
  inla_list$alpha2 <- rep(c(NA, 1), each = n)          ##-- Intercept Y2

  if(!is.null(area)) {
    inla_list$psi_gamma <- c(data[[area]], rep(NA, n))   ##-- gamma*psi
    inla_list$psi <- c(rep(NA, n), data[[area]])         ##-- psi

    inla_list$phi1 <- inla_list$psi_gamma                ##-- psi_1
    inla_list$phi2 <- inla_list$psi                      ##-- psi_2
  }

  inla_list$X <- matrix(NA, nrow = 2*n, ncol = n_covs1 + n_covs2)

  inla_list$X[1:n, 1:n_covs1] <- as.matrix(X1)
  inla_list$X[(n+1):(2*n), (n_covs1+1):(n_covs1 + n_covs2)] <- as.matrix(X2)
  colnames(inla_list$X) <- c(paste(colnames(X1), 1, sep = "_"), paste(colnames(X2), 2, sep = "_"))

  ##-- Fit shared
  f_s <- Y ~ -1 + alpha1 + alpha2 + X

  if(!is.null(area) & !is.null(neigh)) {
    if(shared) {
      f_s <- update(f_s, ~ . + f(psi,
                                 model = "besag",
                                 graph = Ws,
                                 hyper = list(prec =
                                                list(prior = "loggamma",
                                                     param = prior_tau_s))) +
                      f(psi_gamma,
                        copy = "psi",
                        range = c(0, Inf),
                        hyper = list(beta =
                                       list(fixed = FALSE,
                                            prior = "normal",
                                            param = prior_gamma)))
      )
    }

    if(specific_1) {
      f_s <- update(f_s, ~ . + f(phi1,
                                 model = "besag",
                                 graph = W1,
                                 hyper = list(prec =
                                                list(prior = "loggamma",
                                                     param = prior_tau_1)))
      )
    }

    if(specific_2) {
      f_s <- update(f_s, ~ . + f(phi2,
                                 model = "besag",
                                 graph = W2,
                                 hyper = list(prec =
                                                list(prior = "loggamma",
                                                     param = prior_tau_2)))
      )
    }
  }

  time_start_inla <- Sys.time()
  args <- list(...)
  args$control.compute$config <- TRUE
  args$control.inla$strategy <- "laplace"

  inla_aux <- function(...) inla(formula = f_s, family = family, data = inla_list, E = as.vector(E), ...)
  mod <- do.call(what = inla_aux, args = args)

  model_sample <- inla.posterior.sample(result = mod, n = nsamp, use.improved.mean = TRUE)
  hyperpar_samp <- inla.hyperpar.sample(result = mod, n = nsamp, improve.marginals = TRUE)
  time_end_inla <- Sys.time()

  id_fixed <- paste0(c("alpha1", "alpha2", colnames(inla_list$X)), ":", 1)

  if(!is.null(area) & !is.null(neigh)) {
    phi1 <- paste0("phi1:", 1:n)
    phi2 <- paste0("phi2:", 1:n)
    psi <- paste0("psi:", 1:n)
    psi_gamma <- paste0("psi_gamma:", 1:n)
    gamma <- paste0("gamma:", 1)
  }

  beta_samp <- do.call(args = lapply(model_sample, select_marginal, ids = id_fixed), what = "rbind")
  colnames(beta_samp) <- c("alpha1", "alpha2", colnames(inla_list$X))

  latent_sample <- c()

  if(!is.null(area) & !is.null(neigh)) {
    if(specific_1) {
      phi1_sample <- do.call(args = lapply(model_sample, select_marginal, ids = phi1), what = "rbind")
      colnames(phi1_sample) <- phi1

      latent_sample <- cbind(latent_sample, phi1_sample)
    }
    if(specific_2) {
      phi2_sample <- do.call(args = lapply(model_sample, select_marginal, ids = phi2), what = "rbind")
      colnames(phi2_sample) <- phi2

      latent_sample <- cbind(latent_sample, phi2_sample)
    }
    if(shared) {
      psi_sample <- do.call(args = lapply(model_sample, select_marginal, ids = psi), what = "rbind")
      colnames(psi_sample) <- psi
      psi_gamma_sample <- do.call(args = lapply(model_sample, select_marginal, ids = psi_gamma), what = "rbind")
      colnames(psi_gamma_sample) <- psi_gamma
      delta_sample <- matrix(get_trans_samp(fun = sqrt, marg = mod$marginals.hyperpar$`Beta for psi_gamma`, n = nsamp, trunc = TRUE), ncol = 1)
      colnames(delta_sample) <- "Delta"
      psi_precision_sample <- hyperpar_samp[, "Precision for psi"]/hyperpar_samp[, "Beta for psi_gamma"]

      hyperpar_samp <- cbind(hyperpar_samp, delta_sample, `Precision for psi SC` = psi_precision_sample)
      latent_sample <- cbind(latent_sample, psi_sample, psi_gamma_sample, psi_precision_sample)
    }

    sample <- cbind(hyperpar_samp, beta_samp, latent_sample)
  } else {
    sample <- cbind(hyperpar_samp, beta_samp)
  }

  ##-- Time
  time_tab <- data.frame(INLA = as.numeric(difftime(time1 = time_end_inla, time2 = time_start_inla, units = "secs")),
                         correction = as.numeric(difftime(time1 = time_end_correction, time2 = time_start_correction, units = "secs")),
                         total = as.numeric(difftime(time1 = Sys.time(), time2 = time_start, units = "secs")))

  ##-- Return
  out <- list()
  out$sample <- sample
  out$summary_fixed <- chain_summary(sample = beta_samp)
  out$summary_hyperpar <- chain_summary(sample = hyperpar_samp)
  if(is.null(area)) `<-`(out$summary_random, chain_summary(sample = latent_sample)) else `<-`(out$summary_random, NULL)
  out$out <- mod

  out$time <- time_tab

  return(out)
}
