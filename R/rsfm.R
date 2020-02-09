#' @title Restricted Spatial Frailty Model
#'
#' @description Fit a Restricted Spatial Frailty model
#'
#' @usage rsfm(data, time, status, covariates, area, model, family, proj = "none", ...)
#'
#' @param data data.frame containing, at least, \code{time}, \code{status}, \code{covariates}, \code{area} list
#' @param time For right censored data, this is the follow up time. For interval data, this is the starting time for the interval.
#' @param time2 Ending time for the interval censured data.
#' @param status The status indicator, 1 = observed event, 0 = right censored event, 2 = left censored event, 3 = interval censored event.
#' @param covariates Vector of covariates names
#' @param intercept TRUE or FALSE
#' @param area Column representing the region of each individual
#' @param model Spatial model adopted
#' @param W Adjacency matrix
#' @param family 'exponential', 'weibull', 'weibullcure', 'loglogistic', 'gamma', 'lognormal' or 'pwe'
#' @param proj 'none', 'rhz', 'hh' or 'spock'
#' @param fast To use the reduction operator
#' @param nsamp Sample size to use the projection approach
#' @param approach 'inla' or 'bugs'
#' @param ... Other parameters used in ?inla or ?bugs
#'
#' @examples
#' set.seed(1)
#'
#' ##-- Spatial structure
#' data("neigh_RJ")
#'
#' ##-- Individuals and regions
#' n_reg <- length(neigh_RJ)
#' n_id <- sample(x = 3:5, size = n_reg, replace = T)
#' coefs <- c(0.3, -0.3)
#' tau <- 0.75 # Scale of spatial effect
#'
#' ##-- Data
#' data <- rsurv(n_id = n_id,
#'               coefs = coefs, cens = 0.5, scale = FALSE,
#'               cens_type = "right", hazard = "weibull",
#'               hazard_params = hazard_params <- list(weibull = list(alpha = 1.2, variant = 0)),
#'               spatial = "ICAR",
#'               neigh = neigh_RJ, tau = tau, confounding = "linear", proj = "none")
#'
#' ##-- Models
#' weibull_inla <- rsfm(data = data, time = "L", status = "status",
#'                      covariates = c("X1", "X2"), intercept = TRUE,
#'                      family = "weibull", proj = "rhz", nsamp = 1000, approach = "inla")
#'
#' rsfm_inla <- rsfm(data = data, time = "L", status = "status", area = "reg",
#'                   covariates = c("X1", "X2"), intercept = TRUE,
#'                   model = "restricted_besag", neigh = neigh_RJ,
#'                   family = "weibull", proj = "rhz", nsamp = 1000, approach = "inla")
#'
#' weibull_inla$unrestricted$summary_fixed
#' rsfm_inla$unrestricted$summary_fixed
#' rsfm_inla$restricted$summary_fixed
#'
#' @return Restricted model
#'
#' @export

rsfm <- function(data, time, time2 = NULL, status, covariates, intercept = TRUE, area = NULL,
                 model = NULL, neigh = NULL,
                 family, proj = "none", fast = TRUE, nsamp = 1000,
                 approach = "inla", ...) {

  if(is.null(covariates)) stop("You must provide at least one covariate")

  if(!is.null(area)) {
    W <- nb2mat(neighbours = poly2nb(neigh), style = "B")
  } else {
    W <- NULL
  }

  ##-- INLA
  if(approach == "inla") {
    f_fixed <- paste(covariates, collapse = " + ")
    if(!intercept) f_fixed <- paste("-1 +", f_fixed)

    if(!is.null(area)) {
      f_random <- sprintf("f(%s, model = '%s', graph = %s)", area, model, "W")
      f_pred <- paste(f_fixed, f_random, sep = " + ")
    } else{
      f_pred <- f_fixed
    }

    if(is.null(time2)) {
      f <- sprintf("INLA::inla.surv(time = %s, event = %s) ~ %s", time, status, f_pred)
    } else {
      f <- sprintf("INLA::inla.surv(time = %s, time2 = %s, event = %s) ~ %s", time, time2, status, f_pred)
    }

    f <- as.formula(f)

    out <- rsfm_inla(f, data, W = W, family, proj = proj, fast = fast, nsamp = nsamp, ...)
  }

  ##-- BUGS
  if(approach == "bugs") {
    data[[status]] <- ifelse(data[[status]] == 0, data[[time]], 0)
    data[[time]][data[[status]] > 0] <- NA_real_

    bugs_data <- list()

    ##-- Time and status
    bugs_data$t <- data[[time]]
    bugs_data$status <- data[[status]]

    ##-- Fixed effects
    bugs_data$N <- nrow(data)

    f_fixed <- paste('~', paste(covariates, collapse = " + "))
    if(!intercept) f_fixed <- paste("-1 +", f_fixed)
    df_covariates <- model.frame(formula = f_fixed, data = data)

    beta_names <- names(df_covariates)

    for(i in 1:length(beta_names)){
      bugs_data[[beta_names[i]]] <- df_covariates[[beta_names[i]]]
    }

    ##-- Random effects
    if(!is.null(area)){
      bugs_data$n <- length(unique(data[[area]]))
      bugs_data$area <- data[[area]]

      bugs_data$adj <- unlist(apply(X = W, MARGIN = 1, which))
      bugs_data$num <- colSums(W)
      bugs_data$weights <- 1 + 0*bugs_data$adj
    }

    ##-- Make model
    model_file <- tempfile(fileext = ".txt")
    unlink(model_file)

    if(!is.null(area)) {
      ##-- Model
      linear_pred <- paste0("beta_intercept + ", paste0(paste0("beta_", beta_names, " * ", beta_names, "[i]", collapse = " + ")), " + S[area[i]]")
      random_effect <- "S[1 : n] ~ car.normal(adj[], weights[], num[], tau);"

      ##-- Prior
      priors_fixed <- paste0(c("beta_intercept", paste0("beta_", beta_names)), " ~ dnorm(0.0, 0.0001);", collapse = " \n ")
      priors_hyper <- "tau ~ dgamma(0.001, 0.001); \n alpha ~ dgamma(1.0, 0.001);"

      ##-- Initials
      l_init <- vector(mode = "list", length = 3 + length(beta_names))
      names(l_init) <- c("beta_intercept", paste0("beta_", beta_names), "alpha", "tau")
      l_init[1:length(l_init)] <- c(rep(0, length(beta_names) + 1), 1, 1)

      inits <- function() l_init

      ##-- Parameters to save
      parameters <- c("alpha", "tau", "beta_intercept", paste0("beta_", beta_names), "S")
    } else {
      ##-- Model
      linear_pred <- paste0("beta_intercept + ", paste0(paste0("beta_", beta_names, " * ", beta_names, "[i]", collapse = " + ")))
      random_effect <- ""

      ##-- Prior
      priors_fixed <- paste0(c("beta_intercept", paste0("beta_", beta_names)), " ~ dnorm(0.0, 0.0001);", collapse = " \n ")
      priors_hyper <- "alpha ~ dgamma(1.0, 0.001);"

      ##-- Initials
      l_init <- vector(mode = "list", length = 2 + length(beta_names))
      names(l_init) <- c("beta_intercept", paste0("beta_", beta_names), "alpha")
      l_init[1:length(l_init)] <- c(rep(0, length(beta_names) + 1), 1)

      inits <- function() l_init

      ##-- Parameters to save
      parameters <- c("alpha", "beta_intercept", paste0("beta_", beta_names))
    }

    text_model <- sprintf("model {
                           # Model:
                           ## Fixed effects:
                           for (i in 1 : N) {
                               t[i] ~ dweib(alpha, mu[i]) I(status[i], );
                               log(mu[i]) <- %s;
                           }

                           ## Random effects:
                           %s

                           # Priors:
                           ## Fixed effects:
                           %s

                           ## Hyperparameters:
                           %s
                           }", linear_pred, random_effect, priors_fixed, priors_hyper)

    writeLines(text = text_model, con = model_file)

    out <- rsfm_bugs(model = model_file, data = bugs_data, inits = inits, parameters = parameters, covariates = covariates, area = area,
                     proj = proj, fast = fast, nsamp = nsamp, ...)
  }

  return(out)
}
