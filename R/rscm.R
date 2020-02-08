#' @title Restricted Shared Component model
#'
#' @description Fit a Restricted Shared Component model for two diseases
#'
#' @usage rscm(data, Y1, Y2, X1, X2, E1 = NULL, E2 = NULL, W, neigh, area,
#'             proj = "none", nsamp = 1000,
#'             prior_gamma = c(0, 0.1), prior_prec = c(0.5, 0.05),
#'             ...)
#'
#' @param data data.frame containing Y1, Y2, X1, X2
#' @param Y1 Y1 name
#' @param Y2 Y2 name
#' @param X1 X1 names
#' @param X2 X2 names
#' @param E1 Expected counts for Y1
#' @param E2 Expected counts for Y2
#' @param W Adjacency matrix
#' @param neigh Neighborhood structure
#' @param area Areal variable name
#' @param proj 'none' or 'spock'
#' @param nsamp Number of samples to return
#' @param prior_gamma Prior (mean and precision) for the shared component coefficient. Default: N(0, 0.1).
#' @param prior_prec Prior (shape, scale) for the precision parameters
#' @param ... Other parameters used in ?inla
#'
#' @return \item{$sample}{A sample of size nsamp for all parameters in the model}
#' \item{$summary_fixed}{Summary measures for the coefficients}
#' \item{$summary_hyperpar}{Summary measures for hyperparameters}
#' \item{$summary_random}{Summary measures for random quantities}
#' \item{$out}{INLA output}
#' \item{$time}{Time elapsed for fitting the model}
#'
#' @export

rscm <- function(data, Y1, Y2, X1, X2, E1 = NULL, E2 = NULL, W, neigh, area,
                 proj = "none", nsamp = 1000,
                 prior_gamma = c(0, 0.1), prior_prec = c(0.5, 0.05), ...) {
  ##-- Time
  time_start <- Sys.time()

  ##-- Tests
  if(missing(Y1)) stop("You must provide the Y1 name")
  if(missing(Y2)) stop("You must provide the Y2 name")
  if(missing(X1)) stop("You must provide the X1 names")
  if(missing(X2)) stop("You must provide the X2 names")
  if(missing(area)) stop("You must provide the area name")
  if(!proj %in% c("none", "spock")) stop("proj must be 'none' or 'spock'")

  ##-- Setup INLA
  Y <- data[, c(Y1, Y2)]

  covs_s <- union(X1, X2)
  X1 <- data[, X1, drop = FALSE]
  X2 <- data[, X2, drop = FALSE]
  Xs <- data[, covs_s, drop = FALSE]

  n <- nrow(Y)
  n_covs1 <- ncol(X1)
  n_covs2 <- ncol(X2)
  n_covss <- ncol(Xs)

  Ws <- W1 <- W2 <- W

  if(proj == "spock") {
    time_start_correction <- Sys.time()

    neigh1 <- spock(X = as.matrix(cbind(1, X1)), map = neigh)
    W1 <- spdep::nb2mat(neigh1, style = "B", zero.policy = TRUE)

    neigh2 <- spock(X = as.matrix(cbind(1, X2)), map = neigh)
    W2 <- spdep::nb2mat(neigh2, style = "B", zero.policy = TRUE)

    if(ncol(Xs) > 0) {
      neighs <- spock(X = as.matrix(cbind(1, Xs)), map = neigh)
      Ws <- spdep::nb2mat(neighs, style = "B", zero.policy = TRUE)
    }

    time_end_correction <- Sys.time()
  } else {
    time_end_correction <- time_start_correction <- Sys.time()
  }

  if(is.null(E1)) {
    E1 <- rep(1, n)
  } else {
    E1 <- data[, E1, drop = FALSE]
  }
  if(is.null(E2)) {
    E2 <- rep(1, n)
  } else {
    E2 <- data[, E2, drop = FALSE]
  }

  E <- data.frame(E1 = E1, E2 = E2)
  E <- c(E[, 1], E[, 2])

  inla_list <- list()
  inla_list$Y <- cbind(c(Y[, 1], rep(NA, n)),
                       c(rep(NA, n), Y[, 2]))

  inla_list$alpha1 <- rep(c(1, NA), each = n)          ##-- Intercept Y1
  inla_list$alpha2 <- rep(c(NA, 1), each = n)          ##-- Intercept Y2

  inla_list$psi_gamma <- c(data[[area]], rep(NA, n))   ##-- gamma*psi
  inla_list$psi <- c(rep(NA, n), data[[area]])         ##-- psi

  inla_list$phi1 <- inla_list$psi_gamma                ##-- psi_1
  inla_list$phi2 <- inla_list$psi                      ##-- psi_2

  inla_list$X <- matrix(NA, nrow = 2*n, ncol = n_covs1 + n_covs2)

  inla_list$X[1:n, 1:n_covs1] <- as.matrix(X1)
  inla_list$X[(n+1):(2*n), (n_covs1+1):(n_covs1 + n_covs2)] <- as.matrix(X2)
  colnames(inla_list$X) <- c(paste(colnames(X1), 1, sep = "_"), paste(colnames(X2), 2, sep = "_"))

  ##-- Fit shared
  f_s <- Y ~ -1 + alpha1 + alpha2 + X +
    f(psi, model = "besag", graph = Ws, hyper = list(prec = list(prior = "loggamma",
                                                                 param = prior_prec))) +
    f(psi_gamma, copy = "psi", range = c(0, Inf), hyper = list(beta = list(fixed = FALSE,
                                                                           prior = "normal",
                                                                           param = prior_gamma))) +
    f(phi1, model = "besag", graph = W1, hyper = list(prec = list(prior = "loggamma",
                                                                  param = prior_prec))) +
    f(phi2, model = "besag", graph = W2, hyper = list(prec = list(prior = "loggamma",
                                                                  param = prior_prec)))

  if(missing(control.inla)) control.inla <- list(strategy = "laplace") else control.inla$strategy <- "laplace"
  if(missing(control.compute)) control.compute <- list(config = TRUE) else control.compute$config <- TRUE

  time_start_inla <- Sys.time()
  mod <- inla(formula = f_s,
              family = c("poisson", "poisson"),
              data = inla_list,
              E = as.vector(E), control.inla = control.inla,
              control.compute = control.compute, ...)
  model_sample <- inla.posterior.sample(result = mod, n = nsamp, use.improved.mean = TRUE)
  hyperpar_samp <- inla.hyperpar.sample(result = mod, n = nsamp, improve.marginals = TRUE)
  time_end_inla <- Sys.time()

  id_fixed <- paste0(c("alpha1", "alpha2", colnames(inla_list$X)), ":", 1)

  phi1 <- paste0("phi1:", 1:n)
  phi2 <- paste0("phi2:", 1:n)
  psi <- paste0("psi:", 1:n)
  psi_gamma <- paste0("psi_gamma:", 1:n)
  gamma <- paste0("gamma:", 1)

  beta_samp <- do.call(args = lapply(model_sample, select_marginal, ids = id_fixed), what = "rbind")
  colnames(beta_samp) <- c("alpha1", "alpha2", colnames(inla_list$X))
  phi1_sample <- do.call(args = lapply(model_sample, select_marginal, ids = phi1), what = "rbind")
  colnames(phi1_sample) <- phi1
  phi2_sample <- do.call(args = lapply(model_sample, select_marginal, ids = phi2), what = "rbind")
  colnames(phi2_sample) <- phi2
  psi_sample <- do.call(args = lapply(model_sample, select_marginal, ids = psi), what = "rbind")
  colnames(psi_sample) <- psi
  psi_gamma_sample <- do.call(args = lapply(model_sample, select_marginal, ids = psi_gamma), what = "rbind")
  colnames(psi_gamma_sample) <- psi_gamma
  delta_sample <- matrix(get_trans_samp(fun = sqrt, marg = mod$marginals.hyperpar$`Beta for psi_gamma`, n = nsamp, trunc = TRUE), ncol = 1)
  colnames(delta_sample) <- "Delta"
  psi_precision_sample <- hyperpar_samp[, "Precision for psi"]/hyperpar_samp[, "Beta for psi_gamma"]
  hyperpar_samp <- cbind(hyperpar_samp, delta_sample, `Precision for psi SC` = psi_precision_sample)

  latent_sample <- cbind(phi1_sample, phi2_sample, psi_sample, psi_gamma_sample)

  sample <- cbind(hyperpar_samp, beta_samp, latent_sample)

  ##-- Time
  time_tab <- data.frame(INLA = as.numeric(difftime(time1 = time_end_inla, time2 = time_start_inla, units = "secs")),
                         correction = as.numeric(difftime(time1 = time_end_correction, time2 = time_start_correction, units = "secs")),
                         total = as.numeric(difftime(time1 = Sys.time(), time2 = time_start, units = "secs")))

  ##-- Return
  out <- list()
  out$sample <- sample
  out$summary_fixed <- chain_summary(obj = beta_samp)
  out$summary_hyperpar <- chain_summary(obj = hyperpar_samp)
  out$summary_random <- chain_summary(obj = latent_sample)
  out$out <- mod

  out$time <- time_tab

  return(out)
}
