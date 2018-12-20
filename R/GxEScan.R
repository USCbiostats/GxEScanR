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
#' Name of file to write the results to
#' @param skippedFilename
#' Name of file to write info about SNPs that were skipped. If this is blank
#' no file is written. If this is the same as outputFilename, the skipped SNPs
#' are written to the output file along with NA for all tests.
#' @param minMaf
#' Minimum minor allele frequency. Has to be a value between 0.0001 and 0.25
#' @return
#' 0 - success
#' 1 - failure
#' @export
GxEScan <- function(subjectData, geneticData, outputFile, skippedFilename = "", minMaf = 0.05) {
  if (missing(subjectData) == TRUE)
    stop("No subject data specified")
  if (missing(geneticData) == TRUE)
    stop("No genetic data specified")
  if (missing(outputFile) == TRUE)
    stop("No output file specified")
  if (minMaf < 0.0001 | minMaf > 0.25)
    stop("Minimum minor allele frequency must be between 0.0001 and 0.25")
  subjectSubset <- SubsetSubjects(subjectData, geneticData)
  if (is.list(subjectSubset) == FALSE)
    return (1)
#  return (subjectSubset)
  return (GxEScanC(subjectSubset$phenotype, subjectSubset$covariates, subjectSubset$mu, subjectSubset$sigma2, subjectSubset$beta0, subjectSubset$geneIndex, geneticData, outputFile, skippedFilename, minMaf))
}

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
# usesFID -- TRUE data frame has family ID, FALSE otherwise
# subject -- string matrix with family and subject IDs
# phenotype -- numeric vector containing subject phenotypes
# covariates -- numberic matrix with covariate values
# @export
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
  res$phenotype <- as.numeric(subjectData[,phenoCol])
  res$covariates <- as.matrix(subjectData[,firstCov:ncol(subjectData)])
  return (res)
}

# Subset the subjects with complete phenotype and covariate data
# and find their indices in the genetic file
SubsetSubjects <- function(subjectData, geneticData) {
  subjectData <- subjectData[complete.cases(subjectData),]
  subjectTest <- TestSubjectDataR(subjectData)

  if (subjectTest$usesFID != geneticData$usesFID)
    stop("Subject data and genetic data must both use family ID or both must not use family ID")
  
  subIDs <- paste(subjectTest$subjects$FID, subjectTest$subjects$IID, sep = "_gxe_") 
  geneIDs <- paste(geneticData$Samples$FID, geneticData$Samples$SID, sep = "_gxe_")
  indices <- match(subIDs, geneIDs)
  covariates <- subjectTest$covariates[!is.na(indices),]
  mu <- colMeans(covariates)
  sigma2 <- apply(covariates, 2, var)
  covariates <- scale(covariates)
  df <- data.frame(subjectTest$phenotype[!is.na(indices)],covariates)
  colnames(df)[1] <- "Phenotype"
  lr <- glm(Phenotype ~ ., data = df, family = "binomial")
  return (list(phenotype = subjectTest$phenotype[!is.na(indices)],
               covariates = cbind((rep(1.,length(indices))), covariates),
               beta0 = lr$coefficients,
               mu = c(0., mu),
               sigma2 = c(1., sigma2),
               geneIndex = indices[complete.cases(indices)]))
}

AllocateLRModConstantMemory <- function(x, y) {
  n <- nrow(x)
  p <- ncol(x)
  beta <- numeric(p)
  score <- numeric(p)
  w <- numeric(n)
  wInv <- numeric(n)
  yp <- numeric(n)
  zt <- numeric(p)
  k <- numeric(p)
  ql <- matrix(0., nrow = n, ncol = p)
  rtl <- matrix(0., nrow = p, ncol = p)
  logLikelihood <- 0.
  return (list(n = n,
               p = p,
               beta = beta,
               score = score,
               w = w,
               wInv = wInv,
               yp = yp,
               zt = zt,
               k = k,
               ql = ql,
               rtl = rtl,
               logLikelihood = logLikelihood))
}

AllocateLRModFixedMemory <- function(n, p) {
  abx <- numeric(n)
  expabx <- numeric(n)
  expabxp1 <- numeric(n)
  expitabx <- numeric(n)
  yp <- numeric(n)
  zt <- numeric(p)
  k <- numeric(p)
  bt <- numeric(p)
  return (list(abx = abx,
               expabx = expabx,
               expabxp1 = expabxp1,
               expitabx = expitabx,
               yp = yp,
               zt = zt,
               k = k,
               bt = bt))
}

AllocateLRModNotFixedMemory <- function(n, p, q) {
  xrw <- matrix(0., nrow = n, ncol = q)
  beta <- numeric(p + q)
  score <- numeric(p + q)
  zb <- vector(mode = "numeric", length = q)
  bb <- vector(mode = "numeric", length = q)
  h <- matrix(0, nrow = p, ncol = q)
  rtr <- matrix(0, nrow = p, ncol = q)
  t <- matrix(0., nrow = n, ncol = q)
  qr <- matrix(0., nrow = n, ncol = q)
  rbr <- matrix(0., nrow = q, ncol = q)
  logLikelihood <- numeric(1)
  return (list(xrw = xrw,
               beta = beta,
               score = score,
               zb = zb,
               bb = bb,
               h = h,
               rtr = rtr,
               t = t,
               qr = qr,
               rbr = rbr,
               logLikelihood = logLikelihood))
}

