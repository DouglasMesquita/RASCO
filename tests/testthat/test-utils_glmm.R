test_that("utils glmmm", {
  set.seed(123456)

  ##-- Spatial structure
  data("neigh_RJ")

  beta <- c(-0.5, -0.2)
  tau <- 1

  ##-- Data ----
  family <- "poisson"

  testthat::skip_on_appveyor()
  testthat::expect_error(
    rglmm(beta = beta, tau = tau, family = family,
          confounding = "none", neigh = NULL,
          scale = TRUE)
  )
  testthat::expect_error(
    rglmm(beta = beta, tau = tau, family = family,
          confounding = TRUE, neigh = neigh_RJ,
          scale = TRUE)
  )
  testthat::expect_error(
    rglmm(beta = beta, tau = tau, family = family,
          confounding = 'none', neigh = neigh_RJ,
          scale = TRUE, control_family = list(E = c(1, 2)))
  )
  testthat::expect_error(
    rglmm(beta = beta, tau = tau, family = "binomial",
          confounding = 'none', neigh = neigh_RJ,
          scale = TRUE, control_family = list(n = c(1, 2)))
  )


})
