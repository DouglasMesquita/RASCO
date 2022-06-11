test_that("rglmm options", {
  set.seed(123456)

  ##-- Spatial structure
  data("neigh_RJ")
  neigh_RJ_sf <- sf::st_as_sf(neigh_RJ)

  beta <- c(-0.5, -0.2)
  tau <- 1

  ##-- Data ----
  family <- "poisson"
  data <- rglmm(beta = beta, tau = tau, family = family,
                confounding = "none", neigh = neigh_RJ,
                scale = TRUE)

  sglm_mod <- rsglmm(data = data, formula = Y ~ X1 + X2,
                     family = family,
                     proj = "none", nsamp = 1000)

  sglmm_mod <- rsglmm(data = data, formula = Y ~ X1 + X2,
                      family = family,
                      area = "reg", model = "besag", neigh = neigh_RJ,
                      proj = "none", nsamp = 1000)

  sglmm_mod_sf <- rsglmm(data = data, formula = Y ~ X1 + X2,
                         family = family,
                         area = "reg", model = "besag", neigh = neigh_RJ_sf,
                         proj = "none", nsamp = 1000)

  rglmm_rhz <- rsglmm(data = data, formula = Y ~ X1 + X2,
                      family = family,
                      area = "reg", model = "r_besag", neigh = neigh_RJ,
                      proj = "rhz", nsamp = 1000)

  rglmm_rhz_sf <- rsglmm(data = data, formula = Y ~ X1 + X2,
                         family = family,
                         area = "reg", model = "r_besag", neigh = neigh_RJ_sf,
                         proj = "rhz", nsamp = 1000)

  testthat::skip_on_appveyor()
  testthat::skip_on_cran()

  testthat::expect_equal(object = length(sglm_mod$restricted), 0)
  testthat::expect_equal(object = length(sglmm_mod$restricted), 0)
  testthat::expect_equal(object = length(sglmm_mod_sf$restricted), 0)
  testthat::expect_equal(object = length(rglmm_rhz$restricted), 4)
  testthat::expect_equal(object = length(rglmm_rhz_sf$restricted), 4)
  testthat::expect_error(
    rsglmm(data = data,
           family = family,
           area = "reg", model = "besag", neigh = neigh_RJ,
           proj = "none", nsamp = 1000)
  )
  testthat::expect_error(
    rsglmm(data = data, formula = Y ~ X1 + X2,
           family = family,
           area = "reg", model = "r_besag", neigh = neigh_RJ,
           proj = "hughes", nsamp = 1000)
  )
  testthat::expect_message(
    rsglmm(data = data, formula = Y ~ X1 + X2,
           family = family,
           area = "reg", model = "r_besag", neigh = neigh_RJ,
           proj = "hh", approach = "inla", nsamp = 1000)
  )

  data_tst <- rbind(data, data)
  testthat::expect_message(
    rsglmm(data = data_tst, formula = Y ~ X1 + X2,
           family = family,
           area = "reg", model = "r_besag", neigh = neigh_RJ,
           proj = "hh", approach = "inla", nsamp = 1000)
  )
  testthat::expect_message(
    rsglmm(data = data, formula = Y ~ X1 + X2,
           family = family,
           area = "reg", model = "r_besag", neigh = neigh_RJ,
           proj = "hh", approach = "mcmc", nsamp = 1000)
  )
  testthat::expect_message(
    rsglmm(data = data, formula = Y ~ X1 + X2,
           family = family,
           area = "reg", model = "r_besag", neigh = neigh_RJ,
           proj = "rhz", approach = "mcmc", nsamp = 1000, attractive = 10)
  )
})
