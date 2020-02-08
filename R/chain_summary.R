#' @title Summary measures for a(n) (un)restricted model
#'
#' @description Summary measures for a(n) (un)restricted model
#'
#' @param obj A matrix containing parameter's chains
#'
#' @return Measures
#'
#' @export

chain_summary <- function(obj, level = 0.95) {
  if(!is.null(obj)) {
    means <- colMeans(obj)
    medians <- apply(X = obj, MARGIN = 2, FUN = median)
    sds <- apply(X = obj, MARGIN = 2, FUN = sd)
    hpds <- t(apply(X = obj, MARGIN = 2, FUN = function(x) coda::HPDinterval(coda::as.mcmc(x), prob = level)))

    out <- cbind.data.frame(round(means, 6), round(medians, 6), round(sds, 6), round(hpds, 6))
    names(out) <- c("mean", "median", "sd", "lower", "upper")
  } else {
    out <- NULL
  }

  return(out)
}
