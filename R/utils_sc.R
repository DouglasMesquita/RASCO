#' @title rshared
#'
#' @description Generating data from Shared Component model
#'
#' @param alpha_1 intercept for disease 1.
#' @param alpha_2 intercept for disease 2.
#' @param beta_1 coefficients for disease 1.
#' @param beta_2 coefficients for disease 2.
#' @param delta dependence on the shared component.
#' @param tau_1 precision for ICAR specific for disease 1.
#' @param tau_2 precision for ICAR specific for disease 2.
#' @param tau_s precision for ICAR specific for shared component.
#' @param confounding 'none', 'linear', 'quadratic' or 'cubic'.
#' @param neigh neighborhood structure. A \code{SpatialPolygonsDataFrame} object.
#' @param scale scale covariates? TRUE or FALSE.
#'
#' @importFrom stats rnorm rpois
#' @importFrom sp coordinates
#'
#' @export

rshared <- function(alpha_1 = 0, alpha_2 = 0, beta_1 = c(0.1, -0.1), beta_2 = c(0.1, -0.1), delta = 1,
                    tau_1 = 1, tau_2 = 1, tau_s = 1,
                    confounding = 'none', neigh, scale = TRUE){

  if(is.null(neigh)) stop("You must to define neigh (SpatialPolygonsDataFrame object).")
  if(!confounding %in% c("none", "linear", "quadratic", "cubic")) stop("It is a not valid confounding specification. Please try: 'none', 'linear', 'quadratic', 'cubic'.")

  ##-- Covariates
  n <- length(neigh)

  ncov1 <- length(beta_1)
  ncov2 <- length(beta_2)

  X1 <- matrix(rnorm(n = ncov1*n, mean = 0, sd = 1), ncol = ncov1)
  X2 <- matrix(rnorm(n = ncov2*n, mean = 0, sd = 1), ncol = ncov2)

  if(confounding != "none"){
    conf <- sp::coordinates(neigh)[, 1]

    X1[, ncov1] <- switch(confounding,
                          "linear" = conf,
                          "quadratic" = conf^2,
                          "cubic" = conf^3)

    X2[, ncov2] <- switch(confounding,
                          "linear" = conf,
                          "quadratic" = conf^2,
                          "cubic" = conf^3)
  }

  colnames(X1) <- paste0("X1", 1:ncol(X1))
  colnames(X2) <- paste0("X2", 1:ncol(X2))

  if(scale) X1 <- scale(X1)
  if(scale) X2 <- scale(X2)

  ##-- Spatial effects
  W <- nb2mat(neighbours = poly2nb(neigh), style = "B")

  sc <- ricar(W = W, sig = 1/tau_s)
  s1 <- ricar(W = W, sig = 1/tau_1)
  s2 <- ricar(W = W, sig = 1/tau_2)

  spatial_1 <- sc*delta + s1
  spatial_2 <- sc/delta + s2

  ##-- SMR, expected value and counts
  srm1 <- exp(alpha_1 + X1%*%beta_1 + spatial_1)
  E1 <- sample(1:10, n, replace = TRUE)
  Y1 <- rpois(n = n, lambda = E1*srm1)

  srm2 <- exp(alpha_2 + X2%*%beta_2 + spatial_2)
  E2 <- sample(1:10, n, replace = TRUE)
  Y2 <- rpois(n = n, lambda = E2*srm2)

  ##-- Dataset ----
  data <- data.frame(reg = 1:n, Y1 = Y1, Y2 = Y2, E1 = E1, E2 = E2, X1, X2, s1 = s1, s2 = s2, sc = sc, srm1 = srm1, srm2 = srm2)

  return(data)
}

#' @title get_trans_samp
#'
#' @description Generating a sample from the transformed marginal
#'
#' @param marg a marginal from a INLA model.
#' @param fun a transformation function.
#' @param n a number of samples desired.
#' @param trunc TRUE: x > 0.
#' @param method ?inla.tmarginal method.
#'
#' @importFrom INLA inla.tmarginal inla.rmarginal

get_trans_samp <- function(marg, fun, n = 1000, trunc = FALSE, method = "linear") {

  if(trunc) marg <- marg[marg[, 1] > 0, ]

  trans_marg <- INLA::inla.tmarginal(fun = fun, marginal = marg, method = method)
  trans_samp <- INLA::inla.rmarginal(n = n, marginal = trans_marg)

  return(trans_samp)
}

#' @title gamma_prior
#'
#' @description Given mean and variance return a and b
#'
#' @param mean mean of a desired gamma distribution (restore a and b).
#' @param var variance of a desired gamma distribution (restore a and b).
#' @param a shape of a desired gamma distribution (restore mean and variance).
#' @param b scale of a desired gamma distribution (restore mean and variance).

gamma_prior <- function(mean, var, a = NULL, b = NULL){
  if(all(!is.null(c(a, b)))){
    mean <- a/b
    var <- a/b^2

    return(c(mean = mean, var = var))
  } else{
    a <- (mean^2)/var
    b <- mean/var

    return(c(a = a, b = b))
  }
}
