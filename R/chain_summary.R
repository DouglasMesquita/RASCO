#' @title Summary measures for a sample
#'
#' @description Summary measures for a sample
#'
#' @usage chain_summary(sample, level)
#'
#' @param sample a matrix containing parameter's chains.
#' @param level credibility interval level.
#'
#' @importFrom coda HPDinterval as.mcmc
#' @importFrom stats median sd
#'
#' @return a matrix with \code{mean}, \code{median}, \code{sd}, \code{lower} and \code{upper} measures
#'
#' @export

chain_summary <- function(sample, level = 0.95) {
  if(!is.null(sample)) {
    means <- colMeans(sample)
    medians <- apply(X = sample, MARGIN = 2, FUN = median)
    sds <- apply(X = sample, MARGIN = 2, FUN = sd)
    hpds <- t(apply(X = sample, MARGIN = 2, FUN = function(x) coda::HPDinterval(coda::as.mcmc(x), prob = level)))

    out <- cbind.data.frame(round(means, 6), round(medians, 6), round(sds, 6), round(hpds, 6))
    names(out) <- c("mean", "median", "sd", "lower", "upper")
  } else {
    out <- NULL
  }

  return(out)
}
