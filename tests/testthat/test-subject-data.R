context("test-subject-data")

test_that("Tests of subjectData", {
  badSubjectData <- readRDS("BadSubjectData")

  expect_error(GxEScan(), "No subject data specified")
  
  subjectData <- 1L
  expect_error(ProcessSubjectData(subjectData), "subjectData is not a data frame")
  
  subjectData <- data.frame()
  expect_error(ProcessSubjectData(subjectData), "subjectData has no data")
  
  subjectData <- badSubjectData[,c(3,1)]
  expect_error(ProcessSubjectData(subjectData), "first column of subject must be of type \"character\"")
  
  subjectData <- data.frame(badSubjectData[,1], stringsAsFactors = FALSE)
  expect_error(ProcessSubjectData(subjectData), "subjectData has no phenotype data")

  subjectData <- badSubjectData[,1:2]
  expect_error(ProcessSubjectData(subjectData), "subjectData has no phenotype data")

  subjectData <- badSubjectData[,c(1,3)]
  expect_error(ProcessSubjectData(subjectData), "subjectData has no covariate data")
  
  subjectData <- badSubjectData[,1:3]
  expect_error(ProcessSubjectData(subjectData), "subjectData has no covariate data")
  
  subjectData <- badSubjectData[,c(1:3, 5)]
  expect_error(ProcessSubjectData(subjectData), "Phenotype and covariates must be either integer or numeric values")
  
  subjectData <- badSubjectData[,c(1:4)]
  expect_error(ProcessSubjectData(subjectData), "No subject in subjectData have complete data")
})
