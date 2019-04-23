context("test-genetic-data")

# I need to think about adding tests for binary dosage files and vcf files

test_that("Test genetic data is of proper type", {
  expect_error(GxEScan(simSubjectData), "No genetic data specified")
  geneticData <- simGeneticData
  class(geneticData) <- c("badGeneticData")
  expect_error(GxEScan(simSubjectData, geneticData), "class of genetic data in not \"genetic-file-info\"")
})
