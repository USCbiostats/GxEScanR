context("test-gxescan-input")

test_that("Valid GxEScan input", {
  outFile = "junk.out"
  skipFile = "junk.skipped"
  popminMaf = 0.01
  sampleminMaf = 0.05
  binCov = TRUE
  snps = 1L:5L
  
  expect_error(GxEScan(subjectData = simSubjectData,
                       geneticData = simGeneticData,
                       skipFile = skipFile,
                       popminMaf = popminMaf,
                       sampleminMaf = sampleminMaf,
                       binCov = binCov,
                       snps = snps
                       ), "No output file specified")
  expect_error(GxEScan(subjectData = simSubjectData,
                       geneticData = simGeneticData,
                       outFile = 1,
                       skipFile = skipFile,
                       popminMaf = popminMaf,
                       sampleminMaf = sampleminMaf,
                       binCov = binCov,
                       snps = snps
                       ), "outFile must be of type character")
  
  expect_error(GxEScan(subjectData = simSubjectData,
                       geneticData = simGeneticData,
                       outFile = outFile,
                       skipFile = 1,
                       popminMaf = popminMaf,
                       sampleminMaf = sampleminMaf,
                       binCov = binCov,
                       snps = snps
                       ), "skipFile must be of type character")
  
  expect_error(GxEScan(subjectData = simSubjectData,
                       geneticData = simGeneticData,
                       outFile = outFile,
                       skipFile = skipFile,
                       popminMaf = "A",
                       sampleminMaf = sampleminMaf,
                       binCov = binCov,
                       snps = snps
                       ), "popminMaf must be a numeric value")
  expect_error(GxEScan(subjectData = simSubjectData,
                       geneticData = simGeneticData,
                       outFile = outFile,
                       skipFile = skipFile,
                       popminMaf = c(0.1, 0.2),
                       sampleminMaf = sampleminMaf,
                       binCov = binCov,
                       snps = snps
                       ), "popminMaf must be a single value")
  expect_error(GxEScan(subjectData = simSubjectData,
                       geneticData = simGeneticData,
                       outFile = outFile,
                       skipFile = skipFile,
                       popminMaf = matrix(0., 2, 2),
                       sampleminMaf = sampleminMaf,
                       binCov = binCov,
                       snps = snps
                       ), "popminMaf must be a single value")
  expect_error(GxEScan(subjectData = simSubjectData,
                       geneticData = simGeneticData,
                       outFile = outFile,
                       skipFile = skipFile,
                       popminMaf = 0.5,
                       sampleminMaf = sampleminMaf,
                       binCov = binCov,
                       snps = snps
                       ), "popminMaf must be between 0 and 0.25 inclusive")
  expect_error(GxEScan(subjectData = simSubjectData,
                       geneticData = simGeneticData,
                       outFile = outFile,
                       skipFile = skipFile,
                       popminMaf = -0.1,
                       sampleminMaf = sampleminMaf,
                       binCov = binCov,
                       snps = snps
                       ), "popminMaf must be between 0 and 0.25 inclusive")
  
  expect_error(GxEScan(subjectData = simSubjectData,
                       geneticData = simGeneticData,
                       outFile = outFile,
                       skipFile = 1,
                       popminMaf = popminMaf,
                       sampleminMaf = sampleminMaf,
                       binCov = binCov,
                       snps = snps
                       ), "skipFile must be of type character")
  
  expect_error(GxEScan(subjectData = simSubjectData,
                       geneticData = simGeneticData,
                       outFile = outFile,
                       skipFile = skipFile,
                       popminMaf = popminMaf,
                       sampleminMaf = "A",
                       binCov = binCov,
                       snps = snps
                       ), "sampleminMaf must be a numeric value")
  expect_error(GxEScan(subjectData = simSubjectData,
                       geneticData = simGeneticData,
                       outFile = outFile,
                       skipFile = skipFile,
                       popminMaf = popminMaf,
                       sampleminMaf = c(0.1, 0.2),
                       binCov = binCov,
                       snps = snps
                       ), "sampleminMaf must be a single value")
  expect_error(GxEScan(subjectData = simSubjectData,
                       geneticData = simGeneticData,
                       outFile = outFile,
                       skipFile = skipFile,
                       popminMaf = popminMaf,
                       sampleminMaf = matrix(0., 2, 2),
                       binCov = binCov,
                       snps = snps
                       ), "sampleminMaf must be a single value")
  expect_error(GxEScan(subjectData = simSubjectData,
                       geneticData = simGeneticData,
                       outFile = outFile,
                       skipFile = skipFile,
                       popminMaf = popminMaf,
                       sampleminMaf = 0.5,
                       binCov = binCov,
                       snps = snps
                       ), "sampleminMaf must be between 0 and 0.25 inclusive")
  expect_error(GxEScan(subjectData = simSubjectData,
                       geneticData = simGeneticData,
                       outFile = outFile,
                       skipFile = skipFile,
                       popminMaf = popminMaf,
                       sampleminMaf = -0.1,
                       binCov = binCov,
                       snps = snps
                       ), "sampleminMaf must be between 0 and 0.25 inclusive")
  expect_error(GxEScan(subjectData = simSubjectData,
                       geneticData = simGeneticData,
                       outFile = outFile,
                       skipFile = skipFile,
                       popminMaf = popminMaf,
                       sampleminMaf = sampleminMaf,
                       binCov = 1,
                       snps = snps
                       ), "binaryCovariate must be a logical value")
  
})
