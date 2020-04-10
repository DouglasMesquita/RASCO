#' @title Restricted Spatial Frailty Model
#'
#' @description Fit a Restricted Spatial Frailty model
#'
#' @usage rsfm(data, formula, family, area = NULL,
#'             model = NULL, neigh = NULL,
#'             proj = "none", nsamp = 1000,
#'             fast = TRUE, approach = "inla", priors,
#'             ...)
#'
#' @param data a data frame or list containing the variables in the model.
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param family "exponential", "weibull", "weibullcure", "loglogistic", "gamma", "lognormal" or "pwe".
#' @param area areal variable name in \code{data}.
#' @param model spatial model adopted. Examples: "besag", "besag2" or "r_besag". See INLA::inla.list.models() for other models.
#' @param neigh neighborhood structure. A \code{SpatialPolygonsDataFrame} object.
#' @param proj "none" or "rhz".
#' @param nsamp number of samples. Default = 1000.
#' @param fast TRUE to use the reduction operator.
#' @param approach "inla" or "mcmc". "mcmc" has less model options.
#' @param priors a list containing:
#'     \itemize{
#'        \item prior_prec: a vector of size two, representing shape and scale for the gamma distribution applied for model precision
#'     }
#' @param ... other parameters used in ?INLA::inla or ?R2OpenBUGS::bugs
#'
#' @examples
#' set.seed(123456)
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
#'               hazard_params = list(weibull = list(alpha = 1.2, variant = 0)),
#'               spatial = "ICAR",
#'               neigh = neigh_RJ, tau = tau, confounding = "linear", proj = "none")
#'
#' ##-- Models
#' weibull_inla <- rsfm(data = data,
#'                      formula = surv(time = L, event = status) ~ X1 + X2,
#'                      family = "weibull", model = "none",
#'                      proj = "rhz", nsamp = 1000, approach = "inla")
#'
#' rsfm_inla <- rsfm(data = data,
#'                   formula = surv(time = L, event = status) ~ X1 + X2,
#'                   family = "weibull", area = "reg",
#'                   model = "r_besag", neigh = neigh_RJ,
#'                   proj = "rhz", nsamp = 1000, approach = "inla")
#'
#' weibull_inla$unrestricted$summary_fixed
#' rsfm_inla$unrestricted$summary_fixed
#' rsfm_inla$restricted$summary_fixed
#'
#' @return \item{$unrestricted}{A list containing
#'                                \itemize{
#'                                   \item $sample: a sample of size nsamp for all parameters in the model
#'                                   \item $summary_fixed: summary measures for the coefficients
#'                                   \item $summary_hyperpar: summary measures for hyperparameters
#'                                   \item $summary_random: summary measures for random quantities
#'                                 }
#'                              }
#' \item{$restricted}{A list containing
#'                                \itemize{
#'                                   \item $sample: a sample of size nsamp for all parameters in the model
#'                                   \item $summary_fixed: summary measures for the coefficients
#'                                   \item $summary_hyperpar: summary measures for hyperparameters
#'                                   \item $summary_random: summary measures for random quantities
#'                                 }
#'                              }
#'
#' \item{$out}{INLA or BUGS output}
#' \item{$time}{time elapsed for fitting the model}
#'
#' @export

rsfm <- function(data, formula, family,
                 area = NULL, model = NULL, neigh = NULL,
                 proj = "none", nsamp = 1000,
                 fast = TRUE, approach = "inla",
                 priors = list(prior_prec = c(0.5, 0.05)),
                 ...) {

  if(is.null(formula)) stop("You must provide a formula for the fixed effects")
  if(!proj %in% c("none", "rhz")) stop("proj must be 'none' or 'rhz'")

  if(!is.null(area)) {
    W <- nb2mat(neighbours = poly2nb(neigh), style = "B")
  } else {
    W <- NULL
  }

  priors <- append_list(list(prior_prec = c(0.5, 0.05)), priors)
  prior_prec <- priors$prior_prec

  f_fixed <- format(formula)

  ##-- INLA
  if(approach == "inla") {
    if(!is.null(area)) {
      f_random <- sprintf("f(%s,
                             model = '%s',
                             graph = %s,
                             hyper = list(prec = list(prior = 'loggamma',
                                                      param = c(%s, %s))))",
                          area, model, "W", prior_prec[1], prior_prec[2])
      f_pred <- paste(f_fixed, f_random, sep = " + ")
    } else{
      f_pred <- f_fixed
    }

    formula <- gsub(x = f_pred, pattern = "^surv\\(", replacement = "INLA::inla.surv(")
    formula <- as.formula(formula)

    out <- rsfm_inla(data = data, formula = formula, W = W, family, proj = proj, fast = fast, nsamp = nsamp, ...)
  }

  ##-- BUGS
  if(approach == "mcmc") {
    surv_elements <- eval(terms(formula)[[2]])

    event <- surv_elements$event
    time <- surv_elements$time

    data[[event]] <- ifelse(data[[event]] == 0, data[[time]], 0)
    data[[time]][data[[event]] > 0] <- NA_real_

    mcmc_data <- list()

    ##-- Time and event
    mcmc_data$t <- data[[time]]
    mcmc_data$event <- data[[event]]

    ##-- Fixed effects
    mcmc_data$N <- nrow(data)

    df_covariates <- model.frame(formula = formula[-2], data = data)
    beta_names <- names(df_covariates)

    for(i in 1:length(beta_names)){
      mcmc_data[[beta_names[i]]] <- df_covariates[[beta_names[i]]]
    }

    ##-- Random effects
    if(!is.null(area)){
      mcmc_data$n <- length(unique(data[[area]]))
      mcmc_data$area <- data[[area]]

      mcmc_data$adj <- unlist(apply(X = W == 1, MARGIN = 1, which))
      mcmc_data$num <- colSums(W)
      mcmc_data$weights <- 1 + 0*mcmc_data$adj
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
      priors_hyper <- sprintf("tau ~ dgamma(%s, %s); \n alpha ~ dgamma(1.0, 0.001);", prior_prec[1], prior_prec[2])

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
    out <- rsfm_mcmc(model = model_file, data = mcmc_data,
                     inits = inits, parameters = parameters,
                     covariates = covariates, area = area,
                     proj = proj, fast = fast, nsamp = nsamp, ...)
  }

  return(out)
}
