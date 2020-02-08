#' @title SPOCK nb
#'
#' @description Apply the SPOCK approach to transform the neighborhood matrix
#'
#' @param X covariate matrix
#' @param map neighborhood pbject
#'
#' @importFrom sp coordinates
#' @importFrom spdep poly2nb nb2mat
#'
#' @return map_new new neighborhood object
#'
#' @export

spock <- function(X, map){
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
