#' @title Restricted Spatial Generalized Linear Mixed model
#'
#' @description Fit a Restricted Spatial Generalized Linear Mixed model
#'
#' @usage rsglmm(data, formula,
#'               area = NULL, model = NULL, neigh = NULL,
#'               family, proj = "none", nsamp = 1000,
#'               approach = "inla", ...)
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
#' family <- "binomial"
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
#' rglmm_mod <- rsglmm(data = data, formula = Y ~ X1 + X2,
#'                     area = "reg", model = "restricted_besag", neigh = neigh_RJ,
#'                     family = family,
#'                     proj = "rhz", nsamp = 1000)
#'
#' sglm_mod$unrestricted$summary_fixed
#' sglmm_mod$unrestricted$summary_fixed
#' rglmm_mod$unrestricted$summary_fixed
#' rglmm_mod$restricted$summary_fixed
#'
#' sglm_mod$unrestricted$summary_hyperpar
#' sglmm_mod$unrestricted$summary_hyperpar
#'
#' @return Restricted model
#'
#' @export

rsglmm <- function(data, formula,
                   area = NULL, model = NULL, neigh = NULL,
                   family, proj = "none", nsamp = 1000,
                   approach = "inla", ...) {

  if(missing(formula)) stop("You must provide the formula")
  if(!proj %in% c("none", "rhz", "spock")) stop("proj must be 'none', 'rhz' or 'spock'")

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

    out <- rsglmm_inla(f, data, W = W, family, proj = proj, nsamp = nsamp, ...)
  }

  ##-- BUGS
  if(approach == "mcmc") {
    stop("Not implemented yet")
  }

  return(out)
}

