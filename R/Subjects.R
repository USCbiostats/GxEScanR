# Function to test format of subject data table passed to GxEScan
# 
# Function to test format of subject data table passed to GxEScan
# 
# @param subjectData
# Data frame with subject data
# First column must be a character string
# If second column is a character string, the first column is the
# is the family ID and the second column is the subject ID otherwise
# the first column is the subject
# The column following the subject ID is the subject's phenotype
# this column must be a numeric value and not an integer
# The remian columns are the covariates values. These must be numeric
# values. The covariate in the last column is the covariate interacting
# with the gene.
# @return
# res - List contain the following values
# success -- TRUE - data frame is in a valid format, FALSE otherwise
# hasFamilyID -- TRUE data frame has family ID, FALSE otherwise
# subject -- string matrix with family and subject IDs
# phenotype -- numeric vector containing subject phenotypes
# covariates -- numberic matrix with covariate values
# @export
TestSubjectDataR <- function(subjectData) {
  phenoCol = 2;
  firstCov = 3;
  res <- list(success = FALSE, hasFamilyID = FALSE, subjects = NULL, phenotype = NULL, covariates = NULL)
  if (is.data.frame(subjectData) == FALSE) {
    print("SubjectData is not a data frame")
    return (res)
  }
  if (ncol(subjectData) < 3) {
    print("Not enough columns in subject data frame")
    return (res)
  }
  if (is.character(subjectData[,1]) == FALSE) {
    print("First column of SubjectData must be a vector of strings")
    return (res)
  }
  if (is.character(subjectData[,2]) == TRUE) {
    res[["hasFamilyID"]] = TRUE
    if (ncol(subjectData) < 4) {
      print("Not enough columns in subject data frame")
      return (res)
    }
    phenoCol = 3;
    firstCov = 4;
  }
  if (is.numeric(subjectData[,phenoCol]) == FALSE) {
    print("Phenotype column is not numeric")
    return (res)
  }
  for (i in firstCov:ncol(subjectData)) {
    if (is.numeric(subjectData[,i]) == FALSE) {
      print("Not all covariates are numeric")
      return (res)
    }
  }
  res[["success"]] = TRUE
  if (res[["hasFamilyID"]] == TRUE) {
    res[["subjects"]] = as.matrix(subjectData[,c(1:2)])
  } else {
    res[["subjects"]] = as.matrix(subjectData[,1])
  }
  res[["phenotype"]] = subjectData[,phenoCol]
  res[["covariates"]] = as.matrix(subjectData[,firstCov:ncol(subjectData)])
  return (res)
#  df <- read.table("c:/mplink/test/Gold/Original/testfam.fam", stringsAsFactors = FALSE)
#  df2 <- df[,c(2:6)]
#  x <- as.matrix(df[,1:2])
#  return (x)
}