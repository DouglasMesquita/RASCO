#' @title Restricted Spatial Generalized Linear Mixed model
#'
#' @description Fit a Restricted Spatial Generalized Linear Mixed model
#'
#' @usage rsglmm(data, formula, family,
#'               E, n,
#'               area = NULL, model = NULL, neigh = NULL,
#'               proj = "none", nsamp = 1000,
#'               approach = "inla",
#'               ...)
#'
#' @param data an data frame or list containing the variables in the model.
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param family some allowed families are: 'gaussian', 'poisson' and 'binomial'. The family availability will depend on the approach.
#' @param E known component, in the mean for the Poisson likelihoods defined as E = exp(\eqn{\eta}), where \eqn{\eta} is the linear predictor. If not provided it is set to 1.
#' @param n a vector containing the number of trials for the binomial likelihood and variantes, or the number of required successes for the nbinomial2 likelihood. Default value is set to 1..
#' @param area areal variable name in \code{data}.
#' @param model spatial model adopted. Examples: "besag", "besag2" or "restricted_besag". See INLA::inla.list.models() for other models.
#' @param neigh neighborhood structure. A \code{SpatialPolygonsDataFrame} object.
#' @param proj 'none', 'rhz', 'hh' or 'spock'
#' @param nsamp number of desired. samples Default = 1000.
#' @param approach 'inla' or 'mcmc'
#' @param ... other parameters used in ?INLA::inla or ?ngspatial::sparse.sglmm
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
#'                     area = "reg", model = "restricted_besag", neigh = neigh_RJ,
#'                     proj = "rhz", nsamp = 1000)
#'
#' rglmm_spock <- rsglmm(data = data, formula = Y ~ X1 + X2,
#'                       family = family,
#'                       area = "reg", model = "restricted_besag", neigh = neigh_RJ,
#'                       proj = "spock", nsamp = 1000)
#'
#' rglmm_hh <- rsglmm(data = data, formula = Y ~ X1 + X2,
#'                    family = family,
#'                    area = "reg", model = "restricted_besag", neigh = neigh_RJ,
#'                    approach = "mcmc", proj = "hh",
#'                    nsamp = 1000)
#'
#' sglm_mod$unrestricted$summary_fixed
#' sglmm_mod$unrestricted$summary_fixed
#' rglmm_rhz$unrestricted$summary_fixed
#' rglmm_rhz$restricted$summary_fixed
#' rglmm_spock$restricted$summary_fixed
#' rglmm_hh$restricted$summary_fixed
#'
#' @importFrom ngspatial sparse.sglmm
#' @importFrom stats model.matrix update.formula
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
#' \item{$out}{INLA or ngspatial output}
#' \item{$time}{time elapsed for fitting the model}
#'
#' @export

rsglmm <- function(data, formula, family,
                   E = NULL, n = NULL,
                   area = NULL, model = NULL, neigh = NULL,
                   proj = "none", nsamp = 1000,
                   approach = "inla",
                   ...) {

  if(missing(formula)) stop("You must provide the formula")
  if(!proj %in% c("none", "rhz", "hh", "spock")) stop("proj must be 'none', 'rhz', 'hh' or 'spock'")
  if(proj == "hh" & approach == "inla") {
    message("hh is only implemented in MCMC. Changing approach to 'mcmc'")
    approach <- 'mcmc'
  }

  f_fixed <- format(formula)

  if(!is.null(area)) {
    if(proj == "spock" & grepl(x = model, pattern = "restricted")) {
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

  ##-- INLA
  if(approach == "inla") {
    if(!is.null(area)) {
      f_random <- sprintf("f(%s,
                             model = '%s',
                             graph = %s,
                             hyper = list(prec = list(prior = 'loggamma',
                                                      param = c(0.5, 0.0005))))",
                          area, model, "W")
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

    if(family == "nbinomial") family <- "negbinomial"

    if(is.null(list(...)$attractive)) {
      attractive <- round(0.5*(nrow(W)/2))
      message(sprintf("'attractive' parameter not defined. Trying attractive = %s. See ?ngspatial::sparse.sglmm", attractive))
    }

    E <- substitute(E)
    n <- substitute(n)
    if(!is.null(E)) formula <- update.formula(formula, paste0("~ . + offset(log(", deparse(E), "))"))
    if(!is.null(n)) formula <- update.formula(formula, paste0("~ . + offset(log(", deparse(n), "))"))

    out <- rsglmm_mcmc(data = data, formula = formula,
                       family = family, W = W, area = area,
                       proj = proj, ...)

    return(out)
  }

  return(out)
}

