test_that("utils", {
  set.seed(123456)

  ##-- Spatial structure
  data("neigh_RJ")

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

  rglmm_rhz <- rsglmm(data = data, formula = Y ~ X1 + X2,
                      family = family,
                      area = "reg", model = "r_besag", neigh = neigh_RJ,
                      proj = "rhz", nsamp = 1000,
                      control.compute = list(waic = TRUE, dic = TRUE))

  svif <- SVIF(base_model = sglm_mod$unrestricted, model = sglmm_mod$unrestricted)
  svrf <- SVRF(base_model = sglmm_mod$unrestricted, model = rglmm_rhz$restricted)
  waic <- WAIC(rglmm_rhz)
  dic <- DIC(rglmm_rhz)

  mat_aux <- matrix(c(1:10), ncol = 1)
  colnames(mat_aux) <- "X1"
  rownames(mat_aux) <- 1:10

  enl <- mat_aux %e% rep(1:10, each = 10)

  testthat::expect_length(svif, 2)
  testthat::expect_length(svrf, 2)
  testthat::expect_length(waic, 1)
  testthat::expect_length(enl, 100)
})
