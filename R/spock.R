#' @title SPOCK
#'
#' @description SPOCK (SPatial Orthogonal Centroid "K"orrection) approach to transform the neighborhood structure
#'
#' @references Prates, M. O., Assunção, R. M., & Rodrigues, E. C. (2019). \emph{"Alleviating Spatial Confounding for Areal Data Problems by Displacing the Geographical Centroids"}. Bayesian Analysis, 14(2), 623-647.
#'
#' @param X covariate matrix.
#' @param map object of class \code{nb}.
#'
#' @importFrom sp coordinates
#'
#' @return map_new object of class \code{nb}
#'
#' @importFrom stats dist
#'
#' @export

spock <- function(X, map){

  if("sf" %in% class(map)) map <- as(map, "Spatial")

  H <- X%*%solve(t(X)%*%X)%*%t(X)
  n <- ncol(H)

  Px <- cbind((diag(n) - H)%*%sp::coordinates(map))
  map_nb <- spdep::poly2nb(map)
  n_neigh <- sapply(map_nb, length)

  nearest_nb <- apply(as.matrix(dist(Px)), 2, function(x) order(x, decreasing = FALSE))
  k_nearest_nb <- sapply(1:length(map_nb), function(x) c(as.integer(nearest_nb[1:(n_neigh[x] + 1), x]))[-1])
  class(k_nearest_nb) <- "nb"

  A <- spdep::nb2mat(k_nearest_nb, style = "B", zero.policy = TRUE)
  A <- A + t(A)
  A[which(A>0, arr.ind = TRUE)] <- 1

  map_nb <- apply(A, 1, function(x) which(x > 0))
  class(map_nb) <- "nb"

  return(map_nb)
}
