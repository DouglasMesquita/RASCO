test_that("utils sfm", {
  set.seed(123456)

  ##-- Spatial structure
  data("neigh_RJ")

  ##-- Individuals and regions
  n_reg <- length(neigh_RJ)
  n_id <- sample(x = 3:5, size = n_reg, replace = TRUE)
  coefs <- c(0.3, -0.3)
  tau <- 0.75 # Scale of spatial effect

  testthat::skip_on_appveyor()
  testthat::expect_error(
    rsurv(n_id = n_id,
          coefs = coefs, cens = 0.5, scale = FALSE,
          cens_type = "right", hazard = "weibull",
          hazard_params = hazard_params <- list(weibull = list(alpha = 1.2, variant = 0)),
          spatial = "besag",
          neigh = neigh_RJ, tau = tau, confounding = "linear", proj = "none")
  )
  testthat::expect_error(
    rsurv(n_id = n_id,
          coefs = coefs, cens = 0.5, scale = FALSE,
          cens_type = "right", hazard = "gamma",
          hazard_params = hazard_params <- list(weibull = list(alpha = 1.2, variant = 0)),
          spatial = "ICAR",
          neigh = neigh_RJ, tau = tau, confounding = "linear", proj = "none")
  )
  testthat::expect_error(
    rsurv(n_id = n_id,
          coefs = coefs, cens = 0.5, scale = FALSE,
          cens_type = TRUE, hazard = "weibull",
          hazard_params = hazard_params <- list(weibull = list(alpha = 1.2, variant = 0)),
          spatial = "ICAR",
          neigh = neigh_RJ, tau = tau, confounding = "linear", proj = "none")
  )
  testthat::expect_error(
    rsurv(n_id = n_id,
          coefs = coefs, cens = 0.5, scale = FALSE,
          cens_type = "right", hazard = "weibull",
          hazard_params = hazard_params <- list(weibull = list(alpha = 1.2, variant = 0)),
          spatial = "ICAR",
          neigh = NULL, tau = tau, confounding = "linear", proj = "none")
  )
  testthat::expect_error(
    rsurv(n_id = n_id,
          coefs = coefs, cens = 0.5, scale = FALSE,
          cens_type = "right", hazard = "weibull",
          hazard_params = hazard_params <- list(weibull = list(alpha = 1.2, variant = 0)),
          spatial = "ICAR",
          neigh = neigh_RJ, tau = tau, confounding = TRUE, proj = "none")
  )
  testthat::expect_length(
    rsurv(n_id = n_id,
          coefs = coefs, cens = 0.5, scale = FALSE,
          cens_type = "right", hazard = "weibull",
          hazard_params = hazard_params <- list(weibull = list(alpha = 1.2, variant = 0)),
          spatial = "ICAR",
          neigh = neigh_RJ, tau = tau, confounding = "linear", proj = "none"), n = 11
  )
  testthat::expect_length(
    rsurv(n_id = n_id,
          coefs = coefs, cens = 0.5, scale = FALSE,
          cens_type = "left", hazard = "weibull",
          hazard_params = hazard_params <- list(weibull = list(alpha = 1.2, variant = 0)),
          spatial = "ICAR",
          neigh = neigh_RJ, tau = tau, confounding = "linear", proj = "none"), n = 11
  )
  testthat::expect_length(
    rsurv(n_id = n_id,
          coefs = coefs, cens = 0.5, scale = FALSE,
          cens_type = "interval", hazard = "pwe",
          hazard_params = hazard_params <- list(weibull = list(alpha = 1.2, variant = 0)),
          spatial = "ICAR",
          neigh = neigh_RJ, tau = tau, confounding = "quadratic", proj = "rhz"), n = 11
  )
  testthat::expect_length(
    rsurv(n_id = n_id,
          coefs = coefs, cens = 0.5, scale = FALSE,
          cens_type = "interval", hazard = "exponential",
          hazard_params = hazard_params <- list(weibull = list(alpha = 1.2, variant = 0)),
          spatial = "ICAR",
          neigh = neigh_RJ, tau = tau, confounding = "cubic", proj = "none"), n = 11
  )
})
