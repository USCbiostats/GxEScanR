context("test-standardizecovariates")

test_that("Covariates standardize", {
  geneticData <- simGeneticData
  geneticData$filename <- system.file("extdata", simGeneticData$filename, package = "GxEScanR")
  subjectInfo <- ProcessSubjectData(subjectData = simSubjectData)
  subjectInfo <- MatchSubjectsAndGeneticData(subjectInfo, geneticData)

  subjectInfo2 <- subjectInfo;
  subjectInfo2$covariates[,1] <- 1.
  expect_error(StandardizeCovariates(subjectInfo2),
               "At least one of the covariates is constant")
  
  result <- StandardizeCovariates(subjectInfo)
  expect_equal(mean(result$x[,1]), 1.)
  expect_equal(var(result$x[,1]), 0.)
  expect_true(abs(mean(result$x[,2])) < 1e-15)
  expect_equal(var(result$x[,2]), 1.)
  expect_true(abs(mean(result$x[,3])) < 1e-15)
  expect_equal(var(result$x[,3]), 1.)
  expect_true(all(result$eg == result$x[,1:2]))
  expect_true(all(result$cases == result$eg[subjectInfo$phenotype == 1,]))
  expect_true(all(result$controls == result$eg[subjectInfo$phenotype == 0,]))
})
