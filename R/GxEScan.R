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
  subjectSubset <- SubsetSubjects(subjectData, geneticData)
  if (is.list(subjectSubset) == FALSE)
    return (1)
  return (GxEScanC(subjectSubset, geneticData, outputFile))
#  return (subjectSubset)
}

SubsetSubjects <- function(subjectData, geneticData) {
  subjectData <- na.omit(subjectData)
  subjectTest <- TestSubjectDataR(subjectData)
  if (subjectTest$success == FALSE)
    return (1)
  
  if (subjectTest$hasFamilyID == TRUE) {
    if (geneticData$subjects$useFID == FALSE)
      return (1)
    colnames(subjectData)[3] <- "Phenotype"
    subjectData <- subjectData[order(subjectData$Phenotype),]
    m1 <- data.frame(1:nrow(subjectData), subjectData[,c(1,2)], stringsAsFactors = FALSE)
    colnames(m1) <- c("subNum", "FID", "IID")
    m2 <- data.frame(1:nrow(geneticData$subjects$Info), geneticData$subjects$Info)
    colnames(m2) <- c("gLoc", "FID", "IID")
    msg <- merge(m1, m2, by.x = c("FID", "IID"), by.y = c("FID", "IID"))
    msg <- msg[order(msg$subNum),]
    #    colnames(msg)[3] <- "Phenotype"
    phenotypes <- as.numeric(subjectData[c(msg[,"subNum"]),3])
    # covariates <- as.matrix(as.numeric(subjectData[c(msg[,"subNum"]),c(4:ncol(subjectData))]))
    covariates <- as.matrix(subjectData[c(msg[,"subNum"]),c(4:ncol(subjectData))])
  } else {
    colnames(subjectData)[2] <- "Phenotype"
    subjectData <- subjectData[order(subjectData$Phenotype),]
    m1 <- data.frame(1:nrow(subjectData), subjectData[,1], stringsAsFactors = FALSE)
    colnames(m1) <- c("subNum", "IID")
    m2 <- data.frame(1:nrow(geneticData$subjects$Info), geneticData$subjects$Info)
    if (ncol(m2) == 3) {
      colnames(m2) <- c("gLoc", "FID", "IID")
    } else {
      colnames(m2) <- c("gLoc", "IID")
    }
    msg <- merge(m1, m2, by.x = c("IID"), by.y = c("IID"))
    msg <- msg[order(msg$subNum),]

    phenotypes <- as.numeric(subjectData[c(msg[,"subNum"]),2])
#    covariates <- as.matrix(as.numeric(subjectData[c(msg[,"subNum"]),c(3:ncol(subjectData))]))
    covariates <- as.matrix(subjectData[c(msg[,"subNum"]),c(3:ncol(subjectData))])
  }
  gloc <- msg[,c("gLoc")]
  return (list(gLoc = gloc, phenotypes = phenotypes, covariates = covariates))
}  
