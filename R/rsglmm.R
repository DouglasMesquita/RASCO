#' @title Restricted Spatial Generalized Linear Mixed model
#'
#' @description Fit a Restricted Spatial Generalized Linear Mixed model
#'
#' @usage rsglmm(data, formula, E, n,
#'               area = NULL, neigh = NULL, model = NULL,
#'               family, proj = "none", nsamp = 1000,
#'               approach = "inla", ...)
#'
#' @param data data.frame containing the covariates in formula and E, n and area if necessary
#' @param formula Formula for the fixed effects
#' @param E Expected counts for poisson data. Default = 1 for all sample units
#' @param n N trails for binomial data. Default = 1 for all sample units
#' @param area Areal variable name in data
#' @param neigh Neighborhood structure. A \code{SpatialPolygonsDataFrame} object
#' @param model Spatial model adopted: "besag" or "restricted_besag", for example. The model availability will depend on the approach
#' @param family Some allowed families are: 'gaussian', 'poisson' and 'binomial'. The family availability will depend on the approach
#' @param proj 'none', 'rhz', 'hh' or 'spock'
#' @param nsamp Number of samples desired. Default = 1000
#' @param approach 'inla' or 'mcmc'
#' @param ... Other parameters used in ?INLA::inla or ?ngspatial::sparse.sglmm
#'
#' @examples
#' set.seed(1)
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
#'                     area = "reg", model = "besag", neigh = neigh_RJ,
#'                     family = family,
#'                     proj = "none", nsamp = 1000)
#'
#' rglmm_rhz <- rsglmm(data = data, formula = Y ~ X1 + X2,
#'                     area = "reg", model = "restricted_besag", neigh = neigh_RJ,
#'                     family = family,
#'                     proj = "rhz", nsamp = 1000)
#'
#' rglmm_spock <- rsglmm(data = data, formula = Y ~ X1 + X2,
#'                       area = "reg", model = "restricted_besag", neigh = neigh_RJ,
#'                       family = family,
#'                       proj = "spock", nsamp = 1000)
#'
#' rglmm_hh <- rsglmm(data = data, formula = Y ~ X1 + X2,
#'                    area = "reg", model = "restricted_besag", neigh = neigh_RJ,
#'                    family = family,
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
#' @return Restricted model
#'
#' @export

rsglmm <- function(data, formula,
                   E = NULL, n = NULL,
                   area = NULL, neigh = NULL, model = NULL,
                   family, proj = "none", nsamp = 1000,
                   approach = "inla",
                   ...) {

  if(missing(formula)) stop("You must provide the formula")
  if(!proj %in% c("none", "rhz", "hh", "spock")) stop("proj must be 'none', 'rhz', 'hh' or 'spock'")
  if(proj == "hh" & approach == "inla") message("hh is only implemented in INLA. Changing approach to 'inla'\n")

  f_fixed <- format(formula)

  if(!is.null(area)) {
    if(proj == "spock" & grepl(x = model, pattern = "restricted")) {
      X <- model.matrix(object = formula[-2], data = data)

      if(nrow(X) > length(neigh)) {
        message(sprintf("SPOCK still con't deal with different lengths of X and %s. Setting proj = 'rhz' instead. \n", area))
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
    if(proj != "hh") message("Just hh method is available for now. Changing proj to 'hh' \n")
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

