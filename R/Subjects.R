#' Function to test format of subject data table passed to GxEScan
#' 
#' Function to test format of subject data table passed to GxEScan
#' 
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
#' res - List contain the following values
#' usesFID -- TRUE data frame has family ID, FALSE otherwise
#' subject -- string matrix with family and subject IDs
#' phenotype -- numeric vector containing subject phenotypes
#' covariates -- numberic matrix with covariate values
#' @export
TestSubjectDataR <- function(subjectData) {
  phenoCol <- 2;
  firstCov <- 3;
  subjects <- data.frame(FID = character(), IID = character(), stringsAsFactors = FALSE)
  # must be a better way of initializing a blank matrix ??? or not
  res <- list(usesFID = FALSE, subjects = subjects, phenotype = numeric(), covariates = matrix(0,0,0))
  if (is.data.frame(subjectData) == FALSE)
    stop("SubjectData is not a data frame")
  if (ncol(subjectData) < 3)
    stop("Not enough columns in subject data frame")
  if (is.character(subjectData[,1]) == FALSE)
    stop("First column of SubjectData must be a character vector")
  if (is.character(subjectData[,2]) == TRUE) {
    res$usesFID = TRUE
    if (ncol(subjectData) < 4)
      stop("Not enough columns in subject data frame")
    phenoCol <- 3;
    firstCov <- 4;
  }
  if (is.numeric(subjectData[,phenoCol]) == FALSE)
    stop("Phenotype column is not numeric")
  for (i in firstCov:ncol(subjectData)) {
    if (is.numeric(subjectData[,i]) == FALSE)
      stop("Not all covariates are numeric")
  }
  
  if (res$usesFID == TRUE) {
    subjects <- subjectData[,1:2]
  } else {
    subjects <- data.frame(rep("", nrow(subjectData)), subjectData[,1], stringsAsFactors = FALSE)
  }
  colnames(subjects) <- c("FID", "IID")
  res$subjects <- subjects
  res$phenotype <- subjectData[,phenoCol]
  res$covariates <- as.matrix(subjectData[,firstCov:ncol(subjectData)])
  return (res)
}