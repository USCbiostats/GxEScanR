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
#' @param geneticData
#' List with information on reading genetic data
#' This is returned from GetGeneticInfo
#' @param outputFile
#' Name of file to store results in
#' @return
#' 0 - success
#' 1 - failure
#' @export
GxEScan <- function(subjectData, geneticData, outputFile) {
  if (missing(subjectData) == TRUE)
    stop("No subject data specified")
  if (missing(geneticData) == TRUE)
    stop("No genetic data specified")
  if (missing(outputFile) == TRUE)
    stop("No output file specified")
  subjectSubset <- SubsetSubjects(subjectData, geneticData)
  if (is.list(subjectSubset) == FALSE)
    return (1)
  return (GxEScanC(subjectSubset, geneticData, outputFile))
}

# Subset the subjects with complete phenotype and covariate data
# and find their indices in the genetic file
SubsetSubjects <- function(subjectData, geneticData) {
  subjectData <- subjectData[complete.cases(subjectData),]
  subjectTest <- TestSubjectDataR(subjectData)
  
  if (subjectTest$usesFID != geneticData$usesFID)
    stop("Subject data and genetic data must both use family ID on not use family ID")
  
  subIDs <- paste(subjectTest$subjects$FID, subjectTest$subjects$IID, sep = "_gxe_") 
  geneIDs <- paste(geneticData$subjects$FID, geneticData$subjects$IID, sep = "_gxe_")
  indices <- match(subIDs, geneIDs)
  
  return (list(phenotype = subjectTest$phenotype[!is.na(indices)],
               covariates = subjectTest$covariates[!is.na(indices),],
               geneIndex = indices[complete.cases(indices)]))
}  
