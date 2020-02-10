#' @title Restricted Spatial Frailty Model
#'
#' @description Fit a Restricted Spatial Frailty model
#'
#' @usage rsfm(data, formula, area = NULL,
#'             model = NULL, neigh = NULL,
#'             family, proj = "none", fast = TRUE, nsamp = 1000,
#'             approach = "inla", ...)
#'
#' @param data data.frame containing, at least, \code{time}, \code{event}, \code{covariates}, \code{area} list
#' @param formula A formula in the format surv(time, time2, event) ~ X1 + X2
#' @param area Column of data specifying the region of each individual
#' @param model Spatial model adopted: "besag" or  "restricted_besag"
#' @param neigh Neighborhood structure. A \code{SpatialPolygonsDataFrame} object
#' @param family 'exponential', 'weibull', 'weibullcure', 'loglogistic', 'gamma', 'lognormal' or 'pwe'
#' @param proj 'none', 'rhz' or 'spock'
#' @param fast To use the reduction operator
#' @param nsamp Sample size to use the projection approach
#' @param approach 'inla' or 'mcmc'
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
#' n_id <- sample(x = 3:5, size = n_reg, replace = TRUE)
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
#' weibull_inla <- rsfm(data = data,
#'                      formula = surv(time = L, event = status) ~ X1 + X2,
#'                      model = "none", family = "weibull",
#'                      proj = "rhz", nsamp = 1000, approach = "inla")
#'
#' rsfm_inla <- rsfm(data = data, area = "reg",
#'                   formula = surv(time = L, event = status) ~ X1 + X2,
#'                   model = "restricted_besag", neigh = neigh_RJ, family = "weibull",
#'                   proj = "rhz", nsamp = 1000, approach = "inla")
#'
#' weibull_inla$unrestricted$summary_fixed
#' rsfm_inla$unrestricted$summary_fixed
#' rsfm_inla$restricted$summary_fixed
#'
#' @return \item{$unrestricted}{A list containing $sample, $summary_fixed, $summary_hyperpar, $summary_random, $out and $time}
#' \item{$restricted}{A list containing $sample, $summary_fixed, $summary_hyperpar, $summary_random, $out and $time}
#'
#' @export

rsfm <- function(data, formula, area = NULL,
                 model = NULL, neigh = NULL,
                 family, proj = "none", fast = TRUE, nsamp = 1000,
                 approach = "inla", ...) {

  if(is.null(formula)) stop("You must provide a formula for the fixed effects")

  if(!is.null(area)) {
    W <- nb2mat(neighbours = poly2nb(neigh), style = "B")
  } else {
    W <- NULL
  }

  f_fixed <- format(formula)

  ##-- INLA
  if(approach == "inla") {
    if(!is.null(area)) {
      f_random <- sprintf("f(%s, model = '%s', graph = %s)", area, model, "W")
      f_pred <- paste(f_fixed, f_random, sep = " + ")
    } else{
      f_pred <- f_fixed
    }

    f <- gsub(x = f_pred, pattern = "^surv\\(", replacement = "INLA::inla.surv(")
    f <- as.formula(f)

    out <- rsfm_inla(f, data, W = W, family, proj = proj, fast = fast, nsamp = nsamp, ...)
  }

  ##-- BUGS
  if(approach == "mcmc") {
    surv_elements <- eval(terms(formula)[[2]])

    event <- surv_elements$event
    time <- surv_elements$time

    data[[event]] <- ifelse(data[[event]] == 0, data[[time]], 0)
    data[[time]][data[[event]] > 0] <- NA_real_

    bugs_data <- list()

    ##-- Time and event
    bugs_data$t <- data[[time]]
    bugs_data$event <- data[[event]]

    ##-- Fixed effects
    bugs_data$N <- nrow(data)

    df_covariates <- model.frame(formula = formula[-2], data = data)
    beta_names <- names(df_covariates)

    for(i in 1:length(beta_names)){
      bugs_data[[beta_names[i]]] <- df_covariates[[beta_names[i]]]
    }

    ##-- Random effects
    if(!is.null(area)){
      bugs_data$n <- length(unique(data[[area]]))
      bugs_data$area <- data[[area]]

      bugs_data$adj <- unlist(apply(X = W == 1, MARGIN = 1, which))
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
                               t[i] ~ dweib(alpha, mu[i]) I(event[i], );
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

    covariates <- all.vars(formula[-2])
    out <- rsfm_bugs(model = model_file, data = bugs_data, inits = inits, parameters = parameters, covariates = covariates, area = area,
                     proj = proj, fast = fast, nsamp = nsamp, ...)
  }

  return(out)
}
