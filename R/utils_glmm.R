#' @title rglmm
#'
#' @description Generating data from Generalized Linear Mixed models
#'
#' @param beta coefficients parameters.
#' @param tau precision for ICAR.
#' @param family "gaussian", "binomial" or "poisson".
#' @param confounding "none", "linear", "quadratic" or "cubic".
#' @param neigh neighborhood structure. A \code{SpatialPolygonsDataFrame} object.
#' @param scale scale covariates? TRUE or FALSE. See ?scale for more information.
#' @param control_family a list with, at least, a \code{invlink}.
#'
#' @importFrom stats rnorm rpois rbinom
#' @importFrom sp coordinates
#'
#' @export

rglmm <- function(beta = c(0.1, -0.1), tau = 1, family = "gaussian",
                  confounding = 'none', neigh, scale = TRUE,
                  control_family = NULL){

  if(is.null(neigh)) stop("You must to define neigh (SpatialPolygonsDataFrame object).")
  if(!confounding %in% c("none", "linear", "quadratic", "cubic")) stop("It is a not valid confounding specification. Please try: 'none', 'linear', 'quadratic', 'cubic'.")

  control_family <- append_list(control_family_dft(family = family), control_family)

  ##-- Covariates
  n <- length(neigh)
  ncov <- length(beta)

  X <- matrix(rnorm(n = ncov*n, mean = 0, sd = 1), ncol = ncov)

  if(confounding != "none"){
    conf <- sp::coordinates(neigh)[, 1]

    X[, ncov] <- switch(confounding,
                        "linear" = conf,
                        "quadratic" = conf^2,
                        "cubic" = conf^3)
  }

  colnames(X) <- paste0("X", 1:ncol(X))

  if(scale) X <- scale(X)

  ##-- Spatial effects
  W <- nb2mat(neighbours = poly2nb(neigh), style = "B")
  eps <- ricar(W = W, sig = 1/tau)

  ##-- Outcome
  effect <- X%*%beta + eps

  if(family == "gaussian") {
    invlink <- control_family$invlink
    sd <- control_family$sd

    Y <- rnorm(n = n, mean = invlink(effect), sd = sd)

    data <- data.frame(reg = 1:n, Y = Y, X, eps = eps)
  }

  if(family == "poisson") {
    invlink <- control_family$invlink

    E_pois <- control_family$E
    if(length(E_pois) == 1) E_pois <- rep(E_pois, n)
    if(length(E_pois) != n) stop("control_family$E must to be a vector of size n")

    theta <- invlink(effect)

    Y <- rpois(n = n, lambda = E_pois*theta)

    data <- data.frame(reg = 1:n, Y = Y, E = E_pois, X, eps = eps)
  }

  if(family == "binomial") {
    invlink <- control_family$invlink

    n_bin <- control_family$n
    if(length(n_bin) == 1) n_bin <- rep(n_bin, n)
    if(length(n_bin) != n) stop("control_family$n must to be a vector of size n")

    Y <- rbinom(n = n, size = n_bin, prob = invlink(effect))

    data <- data.frame(reg = 1:n, Y = Y, n = n_bin, X, eps = eps)
  }

  return(data)
}

#' @title Default values to control GLMM families
#'
#' @param family "gaussian", "poisson", "binomial".
#'
#' @keywords internal

control_family_dft <- function(family){
  if(family == "gaussian") {
    l_out <- list(invlink = identity, sd = 1)
  }

  if(family == "poisson") {
    l_out <- list(invlink = exp, E = 1)
  }

  if(family == "binomial") {
    l_out <- list(invlink = function(x) exp(x)/(exp(x)+1), n = 1)
  }

  return(l_out)
}
