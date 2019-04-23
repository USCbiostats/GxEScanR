context("test-matchsubjectsandgeneticdata")

test_that("Matching subjects with genetic data", {
  simSubjectData2 <- data.frame(fid = paste("f", c(1:49,51:101), sep = ""), simSubjectData, stringsAsFactors = FALSE)
  simGeneticData2 <- simGeneticData
  simGeneticData2$Samples$FID <- paste("f", 1L:100L, sep = "")
  simGeneticData2$usesFID <- TRUE
  
  processedSubjects <- ProcessSubjectData(simSubjectData)
  expect_error(MatchSubjectsAndGeneticData(subjectData = processedSubjects,
                                           geneticData = simGeneticData2),
               "subjectData and geneticData must either both use family ID or neither use family ID")
  
  processedSubjects$usesFID <- TRUE
  expect_error(MatchSubjectsAndGeneticData(subjectData = processedSubjects,
                                           geneticData = simGeneticData2),
               "No subjects in subject data line up with subjects in genetic data")
  processedSubjects$usesFID <- FALSE

  matchedSubjects <- MatchSubjectsAndGeneticData(subjectData = processedSubjects, geneticData = simGeneticData)
  expect_equal(all(matchedSubjects$indices == c(1:49,51:100)), TRUE)
  expect_equal(all(matchedSubjects$phenotype == simSubjectData[1:99,2]), TRUE)
  expect_equal(all(matchedSubjects$covariates == simSubjectData[1:99,3:4]), TRUE)
  
  processedSubjects <- ProcessSubjectData(simSubjectData2)
  matchedSubjects <- MatchSubjectsAndGeneticData(subjectData = processedSubjects, geneticData = simGeneticData2)
  expect_equal(all(matchedSubjects$indices == c(1:49,51:100)), TRUE)
  expect_equal(all(matchedSubjects$phenotype == simSubjectData2[1:99,3]), TRUE)
  expect_equal(all(matchedSubjects$covariates == simSubjectData2[1:99,4:5]), TRUE)
})
