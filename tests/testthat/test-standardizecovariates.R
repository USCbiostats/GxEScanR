context("test-standardizecovariates")

test_that("Covariates standardize", {
  covariates <- matrix(c(1., 1., 1., 1., 0., 2., 0., 2.), 4, 2)
  expect_error(StandardizeCovariates(covariates),
               "At least one of the covariates is constant")
  covariates[1,1] <- 0.
  covariates[2,1] <- 0.
  sigmaE <- sqrt(c(1./3., 4./3.))
  x1 <- 0.5 / sigmaE[1]
  x <- matrix(c(-x1, -x1, x1, x1, -x1, x1, -x1, x1), 4, 2)
  ans <- StandardizeCovariates(covariates)
  expect_equal(ans$sigmaE, sigmaE)
  expect_equal(ans$x, x)
})
