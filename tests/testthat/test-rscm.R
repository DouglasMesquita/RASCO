test_that("rscm options", {
  set.seed(123456)

  ##-- Spatial structure
  data("neigh_RJ")
  neigh_RJ_sf <- st_as_sf(neigh_RJ)

  ##-- Parameters
  alpha_1 <- 0.5
  alpha_2 <- 0.1
  beta_1 <- c(-0.5, -0.2)
  beta_2 <- c(-0.8, -0.4)
  tau_s <- 1
  tau_1 <- tau_2 <- 10
  delta <- 1.5

  testthat::skip_on_appveyor()
  testthat::skip_on_cran()

  ##-- Data
  data <- rshared(alpha_1 = alpha_1, alpha_2 = alpha_2,
                  beta_1 = beta_1, beta_2 = beta_2,
                  delta = delta,
                  tau_1 = tau_1, tau_2 = tau_2, tau_s = tau_s,
                  confounding = "linear",
                  neigh = neigh_RJ)

  ##-- Models
  simple_model <- rscm(data = data,
                       formula1 = Y1 ~ X11 + X12,
                       formula2 = Y2 ~ X21 + X12,
                       family = c("poisson", "poisson"),
                       E1 = E1, E2 = E2)

  zip_model <- rscm(data = data,
                       formula1 = Y1 ~ X11 + X12,
                       formula2 = Y2 ~ X21 + X12,
                       family = c("poisson", "zeroinflatedpoisson0"),
                       E1 = E1, E2 = E2)

  spoi_model <- rscm(data = data,
                     formula1 = Y1 ~ X11 + X12,
                     formula2 = Y2 ~ X21 + X12,
                     family = c("poisson", "poisson"),
                     E1 = E1, E2 = E2,
                     area = "reg", neigh = neigh_RJ)

  spoi_model_sf <- rscm(data = data,
                        formula1 = Y1 ~ X11 + X12,
                        formula2 = Y2 ~ X21 + X12,
                        family = c("poisson", "poisson"),
                        E1 = E1, E2 = E2,
                        area = "reg", neigh = neigh_RJ_sf)

  rspoi_model <- rscm(data = data,
                      formula1 = Y1 ~ X11 + X12,
                      formula2 = Y2 ~ X21 + X12,
                      family = c("poisson", "poisson"),
                      E1 = E1, E2 = E2,
                      area = "reg", neigh = neigh_RJ,
                      proj = "spock")

  rspoi_model_sf <- rscm(data = data,
                         formula1 = Y1 ~ X11 + X12,
                         formula2 = Y2 ~ X21 + X12,
                         family = c("poisson", "poisson"),
                         E1 = E1, E2 = E2,
                         area = "reg", neigh = neigh_RJ_sf,
                         proj = "spock")

  testthat::expect_equal(object = simple_model$summary_hyperpar, NULL)
  testthat::expect_equal(object = rownames(zip_model$summary_random), NULL)
  testthat::expect_equal(object = rownames(zip_model$summary_hyperpar), "zero-probability parameter for zero-inflated poisson_0[2]")
  testthat::expect_equal(object = nrow(spoi_model$summary_hyperpar), 6)
  testthat::expect_message(
    rscm(data = data,
         formula1 = Y1 ~ X11 + X12,
         formula2 = Y2 ~ X21 + X12,
         family = c("poisson", "poisson"),
         E1 = E1, E2 = E2,
         area = "reg", neigh = neigh_RJ,
         proj = "spock", priors = list(prior_prec = c(0.5, 0.005)))
  )
  testthat::expect_error(
    rscm(data = data,
         formula1 = Y1 ~ X11 + X12,
         family = c("poisson", "poisson"),
         E1 = E1, E2 = E2,
         area = "reg", neigh = neigh_RJ,
         proj = "spock", priors = list(prior_prec = c(0.5, 0.005)))
  )
  testthat::expect_error(
    rscm(data = data,
         formula2 = Y2 ~ X21 + X12,
         family = c("poisson", "poisson"),
         E1 = E1, E2 = E2,
         area = "reg", neigh = neigh_RJ,
         proj = "spock", priors = list(prior_prec = c(0.5, 0.005)))
  )
  testthat::expect_error(
    rscm(data = data,
         formula1 = Y1 ~ X11 + X12,
         formula2 = Y2 ~ X21 + X12,
         family = c("poisson", "poisson"),
         E1 = E1, E2 = E2,
         area = "reg", neigh = neigh_RJ,
         proj = "hh", priors = list(prior_prec = c(0.5, 0.005)))
  )
  testthat::expect_error(
    rscm(data = data,
         formula1 = Y1 ~ X11 + X12,
         formula2 = Y2 ~ X21 + X12,
         family = c("poisson", "poisson"),
         E1 = E1, E2 = E2,
         area = "reg", neigh = neigh_RJ,
         proj = "hh", priors = list(prior_prec = c(0.5, 0.005, 0.01)))
  )
  testthat::expect_true(
    {
      out <- rscm(data = data,
                  formula1 = Y1 ~ X11 + X12,
                  formula2 = Y2 ~ X21 + X12,
                  family = c("poisson", "poisson"),
                  E1 = E1, E2 = E2,
                  area = "reg", neigh = neigh_RJ,
                  random_effects = list(shared = FALSE))

      is.list(out)
    }
  )
  testthat::expect_true(
    {
      out <- rscm(data = data,
         formula1 = Y1 ~ X11 + X12,
         formula2 = Y2 ~ X21 + X12,
         family = c("poisson", "poisson"),
         E1 = E1, E2 = E2,
         area = "reg", neigh = neigh_RJ,
         random_effects = list(specific_1 = FALSE))

      is.list(out)
    }
  )
  testthat::expect_true(
    {
      out <- rscm(data = data,
         formula1 = Y1 ~ X11 + X12,
         formula2 = Y2 ~ X21 + X12,
         family = c("poisson", "poisson"),
         E1 = E1, E2 = E2,
         area = "reg", neigh = neigh_RJ,
         random_effects = list(specific_2 = FALSE))

      is.list(out)
    }
  )
})
