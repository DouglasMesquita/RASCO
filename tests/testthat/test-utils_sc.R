test_that("utils scm", {
  library(spdep)

  set.seed(123456)

  ##-- Spatial structure
  data("neigh_RJ")

  ##-- Parameters
  alpha_1 <- 0.5
  alpha_2 <- 0.1
  beta_1 <- c(-0.5, -0.2)
  beta_2 <- c(-0.8, -0.4)
  tau_s <- 1
  tau_1 <- tau_2 <- 10
  delta <- 1.5

  testthat::skip_on_appveyor()
  testthat::expect_error(
    rshared(alpha_1 = alpha_1, alpha_2 = alpha_2,
            beta_1 = beta_1, beta_2 = beta_2,
            delta = delta,
            tau_1 = tau_1, tau_2 = tau_2, tau_s = tau_s,
            confounding = "linear",
            neigh = NULL)
  )
  testthat::expect_error(
    rshared(alpha_1 = alpha_1, alpha_2 = alpha_2,
            beta_1 = beta_1, beta_2 = beta_2,
            delta = delta,
            tau_1 = tau_1, tau_2 = tau_2, tau_s = tau_s,
            confounding = TRUE,
            neigh = neigh_RJ)
  )
  testthat::expect_warning(
    rshared(alpha_1 = NA, alpha_2 = alpha_2,
            beta_1 = beta_1, beta_2 = beta_2,
            delta = delta,
            tau_1 = tau_1, tau_2 = tau_2, tau_s = tau_s,
            confounding = 'none',
            neigh = neigh_RJ)
  )
})
