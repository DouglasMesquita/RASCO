#' @title Restricted Spatial Generalized Linear Mixed model
#'
#' @description Fit a Restricted Spatial Generalized Linear Mixed model
#'
#' @usage rsglmm(data, formula, family,
#'               E = NULL, n = NULL,
#'               area = NULL, model = NULL, neigh = NULL,
#'               proj = "none", nsamp = 1000, burnin = 5000, lag = 1,
#'               priors = list(prior_prec = list(tau = c(0.5, 0.0005))),
#'               approach = "inla",
#'               ...)
#'
#' @param data a data frame or list containing the variables in the model.
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param family some allowed families are: 'gaussian', 'poisson' and 'binomial'. The family availability will depend on the approach.
#' @param E known component, in the mean for the Poisson likelihoods defined as E = exp(\eqn{\eta}), where \eqn{\eta} is the linear predictor. Default = 1.
#' @param n a vector containing the number of trials for the binomial likelihood, or the number of required successes for the nbinomial2 likelihood. Default value is set to 1.
#' @param area areal variable name in \code{data}.
#' @param model spatial model adopted. Examples: "besag", "besag2" or "r_besag". See INLA::inla.list.models() for other models.
#' @param neigh neighborhood structure. A \code{SpatialPolygonsDataFrame} object.
#' @param proj "none", "rhz", "hh" or "spock"
#' @param nsamp number of samples. Default = 1000.
#' @param burnin burn-in size (just for hh).
#' @param lag lag parameter (just for hh).
#' @param priors a list containing (for now):
#'     \itemize{
#'        \item prior_prec: a list with:
#'        \itemize{
#'            \item tau: a vector of size two containing shape and scale for the gamma distribution applied for \eqn{\tau}
#'        }
#'     }
#'
#' @param approach 'inla' or 'mcmc'
#' @param ... other parameters used in ?INLA::inla or ?ngspatial::sparse.sglmm
#'
#' @details The fitted model is given by
#' \deqn{Y ~ Poisson(E\theta),}
#'
#' \deqn{log(\theta) = X\beta + \psi,}
#'
#' \deqn{\psi ~ ICAR(\tau).}
#'
#' @examples
#' set.seed(123456)
#'
#' ##-- Spatial structure
#' data("neigh_RJ")
#'
#' beta <- c(-0.5, -0.2)
#' tau <- 1
#'
#' ##-- Data ----
#' family <- "poisson"
#' data <- rglmm(beta = beta, tau = tau, family = family,
#'               confounding = "none", neigh = neigh_RJ,
#'               scale = TRUE)
#'
#' ##-- Models ----
#' sglm_mod <- rsglmm(data = data, formula = Y ~ X1 + X2,
#'                    family = family,
#'                    proj = "none", nsamp = 1000)
#'
#' sglmm_mod <- rsglmm(data = data, formula = Y ~ X1 + X2,
#'                     family = family,
#'                     area = "reg", model = "besag", neigh = neigh_RJ,
#'                     proj = "none", nsamp = 1000)
#'
#' rglmm_rhz <- rsglmm(data = data, formula = Y ~ X1 + X2,
#'                     family = family,
#'                     area = "reg", model = "r_besag", neigh = neigh_RJ,
#'                     proj = "rhz", nsamp = 1000)
#'
#' rglmm_spock <- rsglmm(data = data, formula = Y ~ X1 + X2,
#'                       family = family,
#'                       area = "reg", model = "r_besag", neigh = neigh_RJ,
#'                       proj = "spock", nsamp = 1000)
#'
#' rglmm_hh <- rsglmm(data = data, formula = Y ~ X1 + X2,
#'                    family = family,
#'                    area = "reg", model = "r_besag", neigh = neigh_RJ,
#'                    approach = "mcmc",
#'                    proj = "hh", burnin = 5000, nsamp = 1000, lag = 10)
#'
#' sglm_mod$unrestricted$summary_fixed
#' sglmm_mod$unrestricted$summary_fixed
#' rglmm_rhz$unrestricted$summary_fixed
#' rglmm_rhz$restricted$summary_fixed
#' rglmm_spock$restricted$summary_fixed
#' rglmm_hh$restricted$summary_fixed
#'
#' @importFrom methods as
#' @importFrom ngspatial sparse.sglmm
#' @importFrom stats model.matrix update.formula
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
#' \item{$out}{INLA or ngspatial output}
#' \item{$time}{time elapsed for fitting the model}
#'
#' @export

rsglmm <- function(data, formula, family,
                   E = NULL, n = NULL,
                   area = NULL, model = NULL, neigh = NULL,
                   proj = "none", nsamp = 1000, burnin = 5000, lag = 1,
                   priors = list(prior_prec = list(tau = c(0.5, 0.0005))),
                   approach = "inla",
                   ...) {

  if(missing(formula)) stop("You must provide the formula")
  if(!proj %in% c("none", "rhz", "hh", "spock")) stop("proj must be 'none', 'rhz', 'hh' or 'spock'")
  if(proj == "hh" & approach == "inla") {
    message("hh is only implemented in MCMC. Changing approach to 'mcmc'")
    approach <- 'mcmc'
  }

  if("sf" %in% class(neigh)) neigh <- as(neigh, "Spatial")

  f_fixed <- paste0(format(formula), collapse = "")

  if(!is.null(area)) {
    if(proj == "spock" & grepl(x = model, pattern = "restricted_|r_")) {
      X <- model.matrix(object = formula[-2], data = data)

      if(nrow(X) > length(neigh)) {
        message(sprintf("SPOCK still can't deal with different lengths of X and %s. Setting proj = 'rhz' instead.", area))
        proj <- 'rhz'
        W <- nb2mat(neighbours = poly2nb(neigh), style = "B")
      } else{
        neigh <- spock(X = X, map = neigh)
        W <- nb2mat(neighbours = neigh, style = "B")
      }
    } else {
      W <- nb2mat(neighbours = poly2nb(neigh), style = "B", zero.policy = TRUE)
    }
  } else {
    W <- NULL
  }

  if(!is.null(area) && (nrow(W) < nrow(data) & proj == 'hh')) {
    message(sprintf("hh still can't deal with different lengths of X and %s. Setting proj = 'rhz' and approach = 'inla' instead.", area))
    proj <- "rhz"
    approach <- 'inla'
  }

  prior_tau <- priors$prior_prec

  ##-- INLA
  if(approach == "inla") {
    if(!is.null(area)) {
      f_random <- sprintf("f(%s,
                             model = '%s',
                             graph = %s,
                             hyper = list(prec = list(prior = 'loggamma',
                                                      param = c(%s, %s))))",
                          area, model, "W", prior_tau[1], prior_tau[2])
      f_pred <- paste(f_fixed, f_random, sep = " + ")
    } else{
      f_pred <- f_fixed
    }

    formula <- gsub(x = f_pred, pattern = "^surv\\(", replacement = "INLA::inla.surv(")
    formula <- as.formula(formula)

    out <- rsglmm_inla(data = data, formula = formula, W = W, family, proj = proj, nsamp = nsamp, E = E, n = n, ...)
  }

  ##-- BUGS
  if(approach == "mcmc") {
    if(proj != "hh") message("Just hh method is available for now. Changing proj to 'hh'")
    proj <- 'hh'

    if(is.character(family) && family == "nbinomial") family <- "negbinomial"

    if(is.null(list(...)$attractive)) {
      attractive <- round(0.5*(nrow(W)/2))
      message(sprintf("'attractive' parameter not defined. Trying attractive = %s. See ?ngspatial::sparse.sglmm", attractive))
    }

    out <- rsglmm_mcmc(data = data, E = E, n = n, formula = formula,
                       family = family, W = W, area = area,
                       proj = proj, burnin = burnin, nsamp = nsamp, lag = lag,
                       ...)

    return(out)
  }

  return(out)
}

