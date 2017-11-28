#' Function to test format of subject data table passed to GxEScan
#' 
#' Function to test format of subject data table passed to GxEScan
# 
#' @param subjectData
#' Data frame with subject data
#' First column must be a character string
#' If second column is a character string, the first column is the
#' is the family ID and the second column is the subject ID otherwise
#' the first column is the subject
#' The column following the subject ID is the subject's phenotype
#' this column must be a numeric value and not an integer
#' The remian columns are the covariates values. These must be numeric
#' values. The covariate in the last column is the covariate interacting
#' with the gene.
#' @return
#' 0 success
#' 1 failure
#' @export
TestSubjectDataR <- function() {
  print("Testing subjects")
  return (0)
#  df <- read.table("c:/mplink/test/Gold/Original/testfam.fam", stringsAsFactors = FALSE)
#  x <- as.matrix(df[,1:2])
#  return (x)
}