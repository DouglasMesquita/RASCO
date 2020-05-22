#' @title rsurv
#'
#' @description Generating data with spatial confounding
#'
#' @param n_id vector with the number of individuals in each region to be generated.
#' @param coefs vector of coefficients.
#' @param cens proportion of censoring.
#' @param cens_type censoring scheme: "none", "left", "right" or "interval".
#' @param hazard hazard model: "exponential", "weibull" or "pwe" for piecewise exponential.
#' @param hazard_params named list with parameters for the hazard model: hazard_dft().
#' @param spatial spatial model: "none" for the conventional Cox model,
#' "gamma" for an independent gamma frailty,
#' "lognormal" for an independent lognormal frailty,
#' "ICAR" or "BYM" for spatial structured models.
#' @param neigh neighborhood structure. A \code{SpatialPolygonsDataFrame} object.
#' @param tau precision for ICAR and BYM models.
#' @param confounding "none", "linear", "quadratic" or "cubic".
#' @param proj "none", "rhz", "hh" or "spock".
#' @param sd_x standard deviation to generating confounding.
#' @param X matrix of covariates. Default = NULL.
#' @param scale scale X. TRUE or FALSE.
#'
#' @importFrom stats rnorm rgamma
#' @importFrom sp coordinates
#'
#' @export

rsurv <- function(n_id, coefs = c(0.1, 0.4, -0.3),
                  cens = 0, cens_type = "interval",
                  hazard = "weibull",
                  hazard_params = hazard_dft(),
                  spatial = "ICAR", neigh = NULL, tau = 1,
                  confounding = "none", proj = "none", sd_x = 0,
                  X = NULL, scale = TRUE) {

  ##-- Initial checking ----
  if(!(spatial %in% c("none", "gamma", "lognormal", "ICAR", "BYM")))
    stop("It is a not valid spatial model. Please try: 'none', 'gamma', 'lognormal', 'ICAR' or 'BYM'")
  if(!(hazard %in% c("weibull", "exponential", "pwe")))
    stop("It is a not valid hazard model. Please try: 'weibull', 'exponential' or 'pwe'")
  if(!(cens_type %in% c("none", "right", "left", "interval")))
    stop("It is a not valid censoring scheme. Please try: 'right', 'left' or 'interval'")
  if(is.null(neigh)) stop("You must to define neigh (SpatialPolygonsDataFrame object).")
  if(!confounding %in% c("none", "linear", "quadratic", "cubic")) stop("It is a not valid confounding specification. Please try: 'none', 'linear', 'quadratic', 'cubic'.")

  ##-- Appending lists ----
  hazard_params <- append_list(hazard_dft(), hazard_params)

  ##-- General objects ----
  frailty <- eps <- eps_eff <- eps_ort <- NA
  W <- nb2mat(neighbours = poly2nb(neigh), style = "B")

  ##-- Individuals and regions ----
  n_reg <- nrow(W)
  N <- sum(n_id)
  pos_reg <- rep(1:n_reg, n_id)

  ##-- Spatial effects ----
  if((spatial %in% c("ICAR", "BYM") | confounding != "none")){
    eps <- ricar(W = W, sig = 1/tau)

    if(is.null(X) & confounding != "none"){
      if(!is.null(neigh)){
        conf_var <- scale(rnorm(n = N, rowSums(sp::coordinates(neigh)[pos_reg, ]), sd = sd_x))
      } else{
        conf_var <- rnorm(n = N, mean = eps[pos_reg], sd = sd_x)
      }

      conf_var <- switch(confounding,
                         "linear" = conf_var,
                         "quadratic" = conf_var^2,
                         "cubic" = conf_var^3)
    }
  }

  ##-- Frailty effects ----
  if(spatial %in% c("gamma", "lognormal", "BYM")){
    frailty <- hazard_params$frailty$frailty
    frailty <- switch(frailty,
                      "gamma" = log(rgamma(n = N, shape = hazard_params$frailty$params$gamma$shape, rate = hazard_params$frailty$params$gamma$rate)),
                      "lognormal" = rnorm(N, mean = hazard_params$frailty$params$lognormal$meanlog, sd = hazard_params$frailty$params$lognormal$sdlog))

    frailty <- scale(frailty[pos_reg], scale = FALSE)
  }

  ##-- Covariates ----
  if(is.null(X)){
    P <- length(coefs)

    if(confounding != "none"){
      X <- matrix(rnorm(n = N*(P-1), mean = rep(0, P-1)), nrow = N, ncol = P-1)
      X <- cbind(X, conf_var)
      rownames(X) <- NULL
    } else{
      X <- matrix(rnorm(n = N*P, mean = rep(0, P)), nrow = N, ncol = P)
      rownames(X) <- NULL
    }
  }

  if(scale) X <- scale(X)

  colnames(X) <- paste0("X", 1:ncol(X))

  ##-- Projecting effects ----
  if(spatial %in% c("ICAR", "BYM")){
    if(proj != "none"){
      mat_proj <- proj_mat(X = X, groups = pos_reg, method = proj)
      Px <- mat_proj$Px
      Px_ort <- mat_proj$Px_ort

      eps_ort <- as.numeric(Px_ort%*%eps)
      eps_eff <- eps_ort[pos_reg]
    } else{
      eps_eff <- eps[pos_reg]
    }

    eps_eff <- scale(eps_eff, scale = FALSE)
    eps <- tapply(eps_eff, pos_reg, `[`, 1)
  }

  ##-- Fixed effects + random effects ----
  effects <- switch(spatial,
                    "BYM" = exp(as.numeric(coefs%*%t(X)) + eps_eff + frailty),
                    "ICAR" = exp(as.numeric(coefs%*%t(X)) + eps_eff),
                    "gamma" = exp(as.numeric(coefs%*%t(X)) + frailty),
                    "lognormal" = exp(as.numeric(coefs%*%t(X)) + frailty),
                    "none" = exp(as.numeric(coefs%*%t(X))))

  ##-- Failure times ----
  times <- switch(hazard,
                  "exponential" = rexpsurv(N = N, rate = effects*hazard_params$exponential$rate),
                  "weibull" = rweibullsurv(N = N, alpha = hazard_params$weibull$alpha, lambda = effects*hazard_params$weibull$lambda, variant = hazard_params$weibull$variant),
                  "pwe" = rpwesurv(N = N, effects = effects, rates = hazard_params$pwe$rates, tgrid = hazard_params$pwe$tgrid))

  ##-- Censoring times ----
  times_status <- switch(cens_type,
                         "none" = rsurv_none(times = times),
                         "right" = rsurv_right(times = times, cens = cens),
                         "left" = rsurv_left(times = times, cens = cens),
                         "interval" = rsurv_interval(times = times, cens = cens))

  ##-- Data ----
  data <- data.frame(reg = pos_reg,
                     id = as.vector(unlist(sapply(n_id, function(x) 1:x))),
                     L = times_status$L,
                     t = times_status$t,
                     R = times_status$R,
                     status = times_status$status,
                     X,
                     eps = eps[pos_reg],
                     eps_ort = eps_ort[pos_reg],
                     frailty = frailty[pos_reg],
                     check.names = F)

  return(data)
}

#' @title surv
#'
#' @description Auxiliar function for survival models
#'
#' @param time time to event (right censoring), or lower limit (interval censoring).
#' @param time2 upper limit (interval censoring).
#' @param event event indicator; 1 = observed event, 0 = right censored event, 2 = left censored event, 3 = interval censored event.
#'
#' @keywords internal

surv <- function(time, time2 = NULL, event) {
  time <- deparse(substitute(time))
  time2 <- ifelse(!is.null(time2), deparse(substitute(time2)), character(0))
  event <- deparse(substitute(event))

  return(list(time = time, time2 = time2, event = event))
}

#' @title rsurv_none
#'
#' @description Generating data without censoring
#'
#' @param times observed time vector.
#'
#' @keywords internal

rsurv_none <- function(times){
  N <- length(times)

  L <- R <- times

  censure_value <- rep(x = 1, times = N)

  times_status <- data.frame(L = L,
                             t = times,
                             R = R,
                             status = censure_value)

  return(times_status)
}

#' @title rsurv_right
#'
#' @description Generating data with right censoring
#'
#' @param times observed time vector.
#' @param cens censored proportion.
#'
#' @importFrom stats runif
#'
#' @keywords internal

rsurv_right <- function(times, cens){
  N <- length(times)

  L <- R <- times

  ##-- Censure ----
  if(cens == 0) {
    cens_times <- rep(max(times) + 1, N)
  } else {
    cens_times <- runif(n = N, min = 0, max = times/cens)
  }

  cens_ids <- cens_times < times
  cens_val <- rep(x = 1, times = N)

  censure <- cens_ids

  L[cens_ids] <- cens_times[cens_ids]
  R[censure] <- Inf

  cens_val[censure] <- 0

  times_status <- data.frame(L = L,
                             t = times,
                             R = R,
                             status = cens_val)

  return(times_status)
}

#' @title rsurv_left
#'
#' @description Generating data with left censoring
#'
#' @param times observed time vector.
#' @param cens censored proportion.
#'
#' @importFrom stats runif
#'
#' @keywords internal

rsurv_left <- function(times, cens){
  N <- length(times)

  L <- R <- times

  ##-- Censure ----
  if(cens == 0) {
    cens_times <- rep(0, N)
  } else {
    cens_times <- runif(n = N, min = 0, max = times/(1-cens))
  }

  cens_ids <- cens_times > times
  cens_val <- rep(x = 1, times = N)
  censure <- cens_ids

  L[censure] <- 0
  R[cens_ids] <- cens_times[cens_ids]

  cens_val[censure] <- 2

  times_status <- data.frame(L = L,
                             t = times,
                             R = R,
                             status = cens_val)

  return(times_status)
}

#' @title rsurv_interval
#'
#' @description Generating data with interval censoring
#'
#' @param times observed time vector.
#' @param cens censored proportion.
#'
#' @importFrom stats rbinom median runif
#'
#' @keywords internal

rsurv_interval <- function(times, cens){
  L <- ceiling(times)
  R <- L + 1

  cens_val <- 3

  times_status <- data.frame(L = L,
                             t = times,
                             R = R,
                             status = cens_val)

  return(times_status)
}

#' @title rexpsurv
#'
#' @description Generate data from Exponential hazard model
#'
#' @param N number of observations.
#' @param rate rate vector.
#'
#' @importFrom stats rexp
#'
#' @keywords internal

rexpsurv <- function(N, rate = rep(1, N)){
  times <- rexp(n = N, rate = rate)

  return(times)
}

#' @title hexpsurv
#'
#' @description Exponential hazards
#'
#' @param t time vector.
#' @param rate rate vector.
#'
#' @keywords internal

hexpsurv <- function(t, rate){
  hazard <- rate

  return(hazard)
}

#' @title cum_hexpsurv
#'
#' @description Exponential cumulative hazards
#'
#' @param t time vector.
#' @param rate rate vector.
#'
#' @keywords internal

cum_hexpsurv <- function(t, rate){
  cum_hazard <- t*rate

  return(cum_hazard)
}

#' @title rweibullsurv
#'
#' @description Generate data from Weibull hazard model
#'
#' @param N number of observations.
#' @param alpha shape parameter.
#' @param lambda scale parameter.
#' @param variant variant (0 or 1).
#'
#' @importFrom stats rweibull
#'
#' @keywords internal

rweibullsurv <- function(N, alpha, lambda = rep(1, N), variant = 1){
  if(variant == 1){
    times <- rweibull(n = N, shape = alpha, scale = 1/lambda)
  } else{
    times <- rweibull(n = N, shape = alpha, scale = lambda^(-1/alpha))
  }

  return(times)
}

#' @title hweibullsurv
#'
#' @description Weibull hazards
#'
#' @param t time vector.
#' @param alpha alpha parameter.
#' @param lambda lambda parameter.
#'
#' @keywords internal

hweibullsurv <- function(t, alpha, lambda){
  hazard <- alpha*lambda*t^(alpha - 1)

  return(hazard)
}

#' @title cum_hweibullsurv
#'
#' @description Weibull cumulative hazards
#'
#' @param t time vector.
#' @param alpha alpha parameter.
#' @param lambda lambda parameter.
#'
#' @keywords internal

cum_hweibullsurv <- function(t, alpha, lambda){
  cum_hazard <- lambda*t^alpha

  return(cum_hazard)
}

#' @title rpwesurv
#'
#' @description Generate data from Piecewise Exponential hazard model
#'
#' @param N number of observations.
#' @param effects the linear predictor
#' @param rates exponential rates.
#' @param tgrid time grid.
#'
#' @importFrom stats runif
#'
#' @keywords internal

rpwesurv <- function(N, effects = rep(1, N), rates, tgrid){
  u <- runif(n = N, min = 0, max = 1)

  h <- diff(tgrid)
  areas <- c(0, cumsum(h*rates))
  acum <- 1 - exp(-areas)
  id <- as.numeric(cut(u, acum, include.lowest = TRUE))

  times <- ifelse(id == 1,
                  -log(1-u)/(effects*rates[1]),
                  tgrid[id] - (areas[id] + log(1-u))/(effects*rates[id]))

  return(times)
}

#' @title hpwesurv
#'
#' @description Piecewise Eponential hazards
#'
#' @param t time vector.
#' @param rates exponential rates.
#' @param tgrid time grid.
#'
#' @keywords internal

hpwesurv <- function(t, rates, tgrid){

  pos <- as.numeric(cut(t, tgrid, include.lowest = TRUE))
  hazard <- rates[pos]

  return(hazard)
}

#' @title cum_hpwesurv
#'
#' @description Piecewise Exponential cumulative hazards
#'
#' @param t time vector.
#' @param rates exponential rates.
#' @param tgrid time grid.
#'
#' @keywords internal

cum_hpwesurv <- function(t, tgrid, rates){

  h <- diff(tgrid)
  n <- length(t)
  pos <- as.numeric(cut(t, tgrid, include.lowest = TRUE))
  area <- cumsum(h*rates)

  if(n == 1){
    if(pos == 1){
      cum_hazard <- t*rates[1]
    } else{
      cum_hazard <- (t-tgrid[pos])*rates[pos] + area[pos-1]
    }
  } else{
    cum_hazard  <- rep(0, n)
    cum_hazard[which(pos == 1)] <- t[which(pos == 1)]*rates[1]
    aux <- which(pos > 1)
    cum_hazard[aux] <- (t[aux] - tgrid[pos[aux]])*rates[pos[aux]] + area[pos[aux] - 1]
  }

  return(cum_hazard)
}

#' @title Default values for hazard models
#'
#' @keywords internal

hazard_dft <- function(){
  l_out <- list(frailty = list(frailty = "gamma",
                               params = list(gamma = list(shape = 50, rate = 50),
                                             lognormal = list(meanlog = 0, sdlog = 0.1))),
                weibull = list(alpha = 1, lambda = 1, variant = 0),
                exponential = list(rate = 1),
                pwe = list(rates = c(0.1, 0.3, 0.5),
                           tgrid = c(0, 3, 6, Inf)))

  return(l_out)
}
