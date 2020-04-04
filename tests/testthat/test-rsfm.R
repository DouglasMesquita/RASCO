test_that("rsfm options", {
  set.seed(123456)

  ##-- Spatial structure
  data("neigh_RJ")

  ##-- Individuals and regions
  n_reg <- length(neigh_RJ)
  n_id <- sample(x = 3:5, size = n_reg, replace = TRUE)
  coefs <- c(0.3, -0.3)
  tau <- 0.75 # Scale of spatial effect

  ##-- Data
  data <- rsurv(n_id = n_id,
                coefs = coefs, cens = 0.5, scale = FALSE,
                cens_type = "right", hazard = "weibull",
                hazard_params = hazard_params <- list(weibull = list(alpha = 1.2, variant = 0)),
                spatial = "ICAR",
                neigh = neigh_RJ, tau = tau, confounding = "linear", proj = "none")

  ##-- Models
  weibull_inla <- rsfm(data = data,
                       formula = surv(time = L, event = status) ~ X1 + X2,
                       family = "weibull", model = "none",
                       proj = "rhz", nsamp = 1000, approach = "inla")

  rsfm_inla <- rsfm(data = data,
                    formula = surv(time = L, event = status) ~ X1 + X2,
                    family = "weibull", area = "reg",
                    model = "r_besag", neigh = neigh_RJ,
                    proj = "rhz", nsamp = 1000, approach = "inla")

  # weibull_mcmc <- rsfm(data = data,
  #                      formula = surv(time = L, event = status) ~ X1 + X2,
  #                      family = "weibull", model = "none",
  #                      proj = "rhz", nsamp = 10, burnin = 0, lag = 1,
  #                      approach = "mcmc")
  #
  # rsfm_mcmc <- rsfm(data = data,
  #                   formula = surv(time = L, event = status) ~ X1 + X2,
  #                   family = "weibull", area = "reg",
  #                   model = "r_besag", neigh = neigh_RJ,
  #                   proj = "rhz", nsamp = 10, burnin = 0, lag = 1,
  #                   approach = "mcmc")

  testthat::skip_on_appveyor()
  testthat::expect_equal(object = length(weibull_inla$restricted), 0)
  testthat::expect_equal(object = length(weibull_mcmc$restricted), 0)
  # testthat::expect_equal(object = length(rsfm_inla$restricted), 4)
  # testthat::expect_equal(object = length(rsfm_mcmc$restricted), 4)
  testthat::expect_error(
    suppressWarnings(rsfm(data = data,
                          formula = L ~ X1 + X2,
                          family = "weibull", model = "none",
                          proj = "rhz", nsamp = 1000, approach = "inla"))
  )
  testthat::expect_error(
    rsfm(data = data,
         family = "weibull", model = "none",
         proj = "rhz", nsamp = 1000, approach = "inla")
  )
  testthat::expect_error(
    rsfm(data = data,
         formula = surv(time = L, event = status) ~ X1 + X2,
         family = "exponential.cure", model = "none",
         proj = "rhz", nsamp = 1000, approach = "inla")
  )
  testthat::expect_error(
    rsfm(data = data,
         formula = surv(time = L, event = status) ~ X1 + X2,
         family = "weibull", model = "none",
         proj = "hh", nsamp = 1000, approach = "inla")
  )
})
