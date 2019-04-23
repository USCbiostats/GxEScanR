context("test-gxescan-input")

test_that("Valid GxEScan input", {
  outFile = "junk.out"
  skippedFile = "junk.skipped"
  popminMaf = 0.01
  sampleminMaf = 0.05
  snps = 1L:5L
  
  expect_error(GxEScan(subjectData = simSubjectData,
                       geneticData = simGeneticData,
                       skippedFilename = skippedFile,
                       popminMaf = popminMaf,
                       sampleminMaf = sampleminMaf,
                       snps = snps
                       ), "No output file specified")
  expect_error(GxEScan(subjectData = simSubjectData,
                       geneticData = simGeneticData,
                       outputFile = 1,
                       skippedFilename = skippedFile,
                       popminMaf = popminMaf,
                       sampleminMaf = sampleminMaf,
                       snps = snps
                       ), "outputFile must be of type character")
  
  expect_error(GxEScan(subjectData = simSubjectData,
                       geneticData = simGeneticData,
                       outputFile = outFile,
                       skippedFilename = 1,
                       popminMaf = popminMaf,
                       sampleminMaf = sampleminMaf,
                       snps = snps
                       ), "skippedFilename must be of type character")
  
  expect_error(GxEScan(subjectData = simSubjectData,
                       geneticData = simGeneticData,
                       outputFile = outFile,
                       skippedFilename = skippedFile,
                       popminMaf = "A",
                       sampleminMaf = sampleminMaf,
                       snps = snps
                       ), "popminMaf must be a numeric value")
  expect_error(GxEScan(subjectData = simSubjectData,
                       geneticData = simGeneticData,
                       outputFile = outFile,
                       skippedFilename = skippedFile,
                       popminMaf = c(0.1, 0.2),
                       sampleminMaf = sampleminMaf,
                       snps = snps
                       ), "popminMaf must be a single value")
  expect_error(GxEScan(subjectData = simSubjectData,
                       geneticData = simGeneticData,
                       outputFile = outFile,
                       skippedFilename = skippedFile,
                       popminMaf = matrix(0., 2, 2),
                       sampleminMaf = sampleminMaf,
                       snps = snps
                       ), "popminMaf must be a single value")
  expect_error(GxEScan(subjectData = simSubjectData,
                       geneticData = simGeneticData,
                       outputFile = outFile,
                       skippedFilename = skippedFile,
                       popminMaf = 0.5,
                       sampleminMaf = sampleminMaf,
                       snps = snps
                       ), "popminMaf must be between 0 and 0.25 inclusive")
  expect_error(GxEScan(subjectData = simSubjectData,
                       geneticData = simGeneticData,
                       outputFile = outFile,
                       skippedFilename = skippedFile,
                       popminMaf = -0.1,
                       sampleminMaf = sampleminMaf,
                       snps = snps
                       ), "popminMaf must be between 0 and 0.25 inclusive")
  
  expect_error(GxEScan(subjectData = simSubjectData,
                       geneticData = simGeneticData,
                       outputFile = outFile,
                       skippedFilename = 1,
                       popminMaf = popminMaf,
                       sampleminMaf = sampleminMaf,
                       snps = snps
                       ), "skippedFilename must be of type character")
  
  expect_error(GxEScan(subjectData = simSubjectData,
                       geneticData = simGeneticData,
                       outputFile = outFile,
                       skippedFilename = skippedFile,
                       popminMaf = popminMaf,
                       sampleminMaf = "A",
                       snps = snps
                       ), "sampleminMaf must be a numeric value")
  expect_error(GxEScan(subjectData = simSubjectData,
                       geneticData = simGeneticData,
                       outputFile = outFile,
                       skippedFilename = skippedFile,
                       popminMaf = popminMaf,
                       sampleminMaf = c(0.1, 0.2),
                       snps = snps
                       ), "sampleminMaf must be a single value")
  expect_error(GxEScan(subjectData = simSubjectData,
                       geneticData = simGeneticData,
                       outputFile = outFile,
                       skippedFilename = skippedFile,
                       popminMaf = popminMaf,
                       sampleminMaf = matrix(0., 2, 2),
                       snps = snps
                       ), "sampleminMaf must be a single value")
  expect_error(GxEScan(subjectData = simSubjectData,
                       geneticData = simGeneticData,
                       outputFile = outFile,
                       skippedFilename = skippedFile,
                       popminMaf = popminMaf,
                       sampleminMaf = 0.5,
                       snps = snps
                       ), "sampleminMaf must be between 0 and 0.25 inclusive")
  expect_error(GxEScan(subjectData = simSubjectData,
                       geneticData = simGeneticData,
                       outputFile = outFile,
                       skippedFilename = skippedFile,
                       popminMaf = popminMaf,
                       sampleminMaf = -0.1,
                       snps = snps
                       ), "sampleminMaf must be between 0 and 0.25 inclusive")

})
