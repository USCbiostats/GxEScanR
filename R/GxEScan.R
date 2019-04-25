#' @useDynLib GxEScanR, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats complete.cases glm var logLik
NULL

# Function to process the subjectData data frame passed to GxEScan
# Checks are performed to make sure data frame is in the proper format.
# Once the data is verified to be in the proper format, the data is
# is broken into 3 groups, fid and iid, phenotype, and covariates.
# Only the rows with complete data are kept. An indicator variable is
# also returned indicated if family id, fid, is used. If fid is not
# used it is still returned but set to ""
ProcessSubjectData <- function(subjectData) {
  if (is.data.frame(subjectData) == FALSE)
    stop("subjectData is not a data frame")
  
  usesfid <- FALSE
  nc <- ncol(subjectData)
  nr <- nrow(subjectData)
  pc <- 2
  
  if (nc == 0 || nr == 0)
    stop("subjectData has no data")
  
  coltype <- sapply(subjectData, class)
  if (coltype[1] != "character")
    stop("first column of subject must be of type \"character\"")
  if (nc == 1)
    stop("subjectData has no phenotype data")
  if (coltype[2] == "character") {
    usesfid <- TRUE
    pc <- 3
    if (nc == 2)
      stop("subjectData has no phenotype data")
    if (nc == 3)
      stop("subjectData has no covariate data")
    
  } else {
    if (nc == 2)
      stop("subjectData has no covariate data")
  }

  colNumType <- coltype[pc:nc] == "numeric" | coltype[pc:nc] == "integer"
  if (all(colNumType) == FALSE)
    stop("Phenotype and covariates must be either integer or numeric values")
  
  completeSubjects <- subjectData[complete.cases(subjectData),]
  if (nrow(completeSubjects) == 0)
    stop("No subject in subjectData have complete data")
  
  if (usesfid) {
    fid <- completeSubjects[,1]
    iid <- completeSubjects[,2]
  } else {
    fid <- rep("", nrow(completeSubjects))
    iid <- completeSubjects[,1]
  }
  phenotype <- as.numeric(completeSubjects[,pc])
  covariates <- sapply(completeSubjects[,(pc + 1):nc], as.numeric)
  
  return (list(subjects = data.frame(fid, iid),
               usesFID = usesfid,
               phenotype = phenotype,
               covariates = as.matrix(covariates)))
}

# Tests and gets SNP indices
# If snpsToFind is an integer vector, it is the indices. It values
# are checked against the genetic file info to verify they are valid.
# If snpsToFind is a character vector, the values are found in the
# genetic file info SNPs data frame and the indices saved.
GetSNPIndices <- function(snpsToFind, snpsInGeneticFile, numSNPs) {
  if (is.vector(snpsToFind, 'character') == TRUE)
    snpIndices <- match(snpsToFind, snpsInGeneticFile)
  else if (is.vector(snpsToFind, 'integer') == FALSE)
    stop("snps is not a character or integer vector")
  else
    snpIndices <- snpsToFind

  if (anyNA(snpIndices))
    stop("Not all SNPs found")
  snpIndices <- sort(snpIndices)
  if (length(snpIndices) != length(unique(snpIndices)))
    stop("SNP values are not unique")
  if (min(snpIndices) < 1)
    stop("snp indices must be positive")
  if (max(snpIndices) > numSNPs)
    stop("Value of largest snp index greater than number of SNPs")
  
  return (snpIndices)  
}

# Matches subjects in subject data to subjects in genetic data
# returns a list that contains the phenotype, covariates, and indices
# for genetic data for subjects that have complete phenotype and
# covariate data that also have genetic data
MatchSubjectsAndGeneticData <- function(subjectData, geneticData) {
  if (subjectData$usesFID != geneticData$usesFID)
    stop("subjectData and geneticData must either both use family ID or neither use family ID")
  
  if (subjectData$usesFID == TRUE) {
    x <- paste(subjectData$subjects$fid, subjectData$subjects$iid, sep = "_gXe_")
    y <- paste(geneticData$Samples$FID, geneticData$Samples$SID, sep = "_gXe_")
    m1 <- match(x, y)
  } else {
    m1 <- match(subjectData$subjects$iid, geneticData$Samples$SID)
  }
  indices = m1[complete.cases(m1)]
  if (length(indices) == 0)
    stop("No subjects in subject data line up with subjects in genetic data")
  phenotype <- subjectData$phenotype[complete.cases(m1)]
  uniquePhenotypes <- sort(unique(phenotype))
  if (length(uniquePhenotypes) > 2)
    stop("There are more than 2 values for the phenotype")
  if (length(uniquePhenotypes) == 1)
    stop("All subjects have the same phenotype")
  if(all(uniquePhenotypes == c(0., 1.)) != TRUE)
    stop("Phenotype values must be coded as 0, 1")
  covariates <- as.matrix(subjectData$covariates[complete.cases(m1),])
  return (list(phenotype = phenotype,
               covariates = covariates,
               indices = indices))
}

# Standardize the covariates.
# This routines standardized the covariates. This makes fitting logistics
# models a little easier. The standard deviation is also returned because
# this is needed to convert the beta estimates back to the original scale.
StandardizeCovariates <- function(covariates) {
  sigmaE = sqrt(diag(var(covariates)))
  if (min(sigmaE) == 0.)
    stop("At least one of the covariates is constant")
  x <- scale(covariates)
  x = matrix(x, nrow(covariates), ncol(covariates))
  return (list(sigmaE = sigmaE,
               x = x))
}

# Get information needed to read in a binary dosage file
GetBDFileInfo <- function(bdInfo) {
  filename <- bdInfo$filename
  format <- integer(2)
  format[1] <- bdInfo$format
  format[2] <- bdInfo$version
  numSub <- integer(1)
  numSub[1] <- bdInfo$NumSamples
  numSNPs <- integer(1)
  numSNPs[1] <- bdInfo$NumSNPs
  bufferSize <- integer(1)
  bufferSize[1] <- 100000001L
  buffer <- integer((bufferSize / 4) + 1)
  sections <- integer(1)
  filesize <- file.info(filename)$size
  #  filesize <- 264829483501
  locations <- numeric(numSNPs + 1)
  sections <- GetLocations(bdInfo$Indices, locations, filesize, bufferSize)
  snpSection <- integer(numSNPs)
  fileLocation <- numeric(sections + 1)
  snpLocation <- integer(numSNPs)
  GetSections(locations, snpSection, fileLocation, snpLocation, bufferSize)
  currentSection <- integer(1)
  currentSection[1] <- -1L
  dosage <- numeric(numSub)
  p0 <- numeric(numSub)
  p1 <- numeric(numSub)
  p2 <- numeric(numSub)
  return (list(filename = filename,
               format = format,
               numSub = numSub,
               numSNPs = numSNPs,
               bufferSize = bufferSize,
               buffer = buffer,
               sections = sections,
               snpSection = snpSection,
               fileLocation = fileLocation,
               snpLocation = snpLocation,
               currentSection = currentSection,
               dosage = dosage,
               p0 = p0,
               p1 = p1,
               p2 = p2))
}

#' GxEScan - Perform a genomewide interaction scan
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
#' This is returned from BianryDosage::GetGeneticInfo
#' @param outputFile
#' Name of file to write the results to
#' @param skippedFilename
#' Name of file to write info about SNPs that were skipped. If this is blank
#' no file is written. If this is the same as outputFilename, the skipped SNPs
#' are written to the output file along with NA for all tests.
#' @param popminMaf
#' Population minimum allele frequency
#' @param sampleminMaf
#' Sample minimum allele frequency
#' @param snps
#' Vector of snps to be tested. Can either me a character vector of SNP names
#' or a vector of indices.
#' @return
#' 0 - success
#' 1 - failure
#' @export
#' @examples
#' subjectData <- simSubjectData
#' geneticData <- simGeneticData
#' geneticData$filename <- system.file("extdata", simGeneticData$filename, package = "GxEScanR")
#' outfile = "stdout"
#' skippedFile = ""
#' popminMaf = 0.05
#' sampleminMaf = 0.1
#' snps = 1L:110L
#' GxEScan(subjectData, geneticData, outfile, skippedFile, popminMaf, sampleminMaf, snps)
GxEScan <- function(subjectData, geneticData, outputFile, skippedFilename = "", popminMaf = 0.01, sampleminMaf = 0.05, snps) {
  if (missing(subjectData) == TRUE)
    stop("No subject data specified")
  subjectInfo <- ProcessSubjectData(subjectData = subjectData)

  # This needs to be update in the BinaryDosage package.
  # geneticData needs to be assigned a class and that
  # class needs to be checked for here
  if (missing(geneticData) == TRUE)
    stop("No genetic data specified")
  if ("genetic-file-info" %in% class(geneticData) == FALSE)
    stop("class of genetic data in not \"genetic-file-info\"")
  
  if (missing(outputFile) == TRUE)
    stop("No output file specified")
  if (is.character(outputFile) == FALSE)
    stop("outputFile must be of type character")

  if (is.character(skippedFilename) == FALSE)
    stop("skippedFilename must be of type character")

  if (is.numeric(popminMaf) == FALSE)
    stop("popminMaf must be a numeric value")
  if (length(popminMaf) != 1)
    stop("popminMaf must be a single value")
  if (popminMaf < 0. | popminMaf > 0.25)
    stop("popminMaf must be between 0 and 0.25 inclusive")
  
  if (is.numeric(sampleminMaf) == FALSE)
    stop("sampleminMaf must be a numeric value")
  if (length(sampleminMaf) != 1)
    stop("sampleminMaf must be a single value")
  if (sampleminMaf < 0. | sampleminMaf > 0.25)
    stop("sampleminMaf must be between 0 and 0.25 inclusive")

  if (missing(snps) == TRUE) {
    snpIndices <- seq(1L:geneticData$NumSNPs)
  } else {
    snpIndices <- GetSNPIndices(snps, geneticData$SNPs, geneticData$NumSNPs)
  }
  subjectInfo <- MatchSubjectsAndGeneticData(subjectInfo, geneticData)
  
  standardizedX <- StandardizeCovariates(subjectInfo$covariates)
  df <- data.frame(y = subjectInfo$phenotype, x = standardizedX$x)
  gge <- glm(y ~ ., data = df, family = "binomial")
  
  standardizedX$x <- cbind(1., standardizedX$x)
  p1d <- AllocateLRModConstantMemory(standardizedX$x, subjectInfo$phenotype)
  p2d <- AllocateLRModFixedMemory(p1d$n, p1d$p)
  p1d$beta <- gge$coefficients
  result <- InitializeLRMod(p1d$n, p1d$p, subjectInfo$phenotype, standardizedX$x,
                            p1d$beta, p1d$score, p1d$w, p1d$wInv,
                            p1d$yp, p1d$zt, p1d$k, p1d$ql, p1d$rtl,
                            p2d$abx, p2d$expabx, p2d$expabxp1, p2d$expitabx, p1d$logLikelihood)
  if (result != 0)
    stop("Error initialize D|E model")
  if (abs(p1d$logLikelihood - logLik(gge)) > 1e-7)
    stop("Error calculating log likelihood")
  rm(df)
  rm(gge)

  snpminmaf <- match(geneticData$SNPs$SNPID[abs(geneticData$SNPInfo$AAF - 0.5) < 0.5 - popminMaf], geneticData$SNPs$SNPID)
  snpminmaf <- snpminmaf[snpminmaf %in% snpIndices]
  bdInfo <- GetBDFileInfo(geneticData)
  
  snpSet <- matrix(data = 0., nrow = length(subjectInfo$indices), ncol = 400)
  p3dg <- AllocateLRModNotFixedMemory(p1d$n, p1d$p, 1)
  p3dgxe <- AllocateLRModNotFixedMemory(p1d$n, p1d$p, 2)
  # Memory needed for adding columns to previous LR model fit
  xr1 <- matrix(data = 0., nrow = p1d$n, ncol = 1)
  xr2 <- matrix(data = 0., nrow = p1d$n, ncol = 2)
  # Memory for results
  loglikelihoods <- matrix(data = 0., nrow = 100, ncol = 3)
  estimates <- matrix(data = 0., nrow = 100, ncol = 2)
  # Memory for SNP info
  snpID <- character(100)
  chromosome <- character(100)
  locations <- integer(100)
  refAllele <- character(100)
  altAllele <- character(100)
  # Number of subjects and number of cases
  numSub = length(subjectInfo$phenotype)
  numCases = sum(subjectInfo$phenotype)
  
  # Count the number of groups and loop over them
  numGroups = ceiling(length(snpminmaf) / 100)
  # Loop over the groups of SNPs
  OpenGxEOutFile(outputFile)
  
  for (i in 1:numGroups) {
    firstSNP <- (i - 1) * 100 + 1
    lastSNP <- min(firstSNP + 99, length(snpminmaf))

    readResult <- ReadSNP(snpminmaf[firstSNP:lastSNP], subjectInfo$indices,
                          bdInfo$filename, bdInfo$format,
                          bdInfo$numSub, bdInfo$numSNPs,
                          bdInfo$bufferSize, bdInfo$buffer,
                          bdInfo$sections, bdInfo$snpSection,
                          bdInfo$fileLocation, bdInfo$snpLocation,
                          bdInfo$currentSection, bdInfo$dosage,
                          bdInfo$p0, bdInfo$p1, bdInfo$p2, snpSet)
    if (readResult == 1)
      stop("Error reading genetic file")
    
    snpID[1:length(firstSNP:lastSNP)] <- geneticData$SNPs$SNPID[snpminmaf[firstSNP:lastSNP]]
    chromosome[1:length(firstSNP:lastSNP)] <- geneticData$SNPs$Chromosome[snpminmaf[firstSNP:lastSNP]]
    locations[1:length(firstSNP:lastSNP)] <- geneticData$SNPs$Location[snpminmaf[firstSNP:lastSNP]]
    refAllele[1:length(firstSNP:lastSNP)] <- geneticData$SNPs$Reference[snpminmaf[firstSNP:lastSNP]]
    altAllele[1:length(firstSNP:lastSNP)] <- geneticData$SNPs$Alternate[snpminmaf[firstSNP:lastSNP]]
    loglikelihoods[,] = NA
    scanResult <- ScanDisease(p1d$n, p1d$p, subjectInfo$phenotype, standardizedX$x,
                              snpSet, snpID, length(firstSNP:lastSNP), sampleminMaf,
                              p1d$beta, p1d$score, p1d$w, p1d$wInv, p1d$yp,
                              p1d$zt, p1d$k, p1d$ql, p1d$rtl,
                              p2d$abx, p2d$expabx, p2d$expabxp1, p2d$expitabx,
                              p2d$yp, p2d$zt, p2d$k, p2d$bt,
                              p3dg$xrw, p3dg$beta, p3dg$score, p3dg$zb, p3dg$bb,
                              p3dg$h, p3dg$rtr, p3dg$t, p3dg$qr, p3dg$rbr, p3dg$logLikelihood,
                              p3dgxe$xrw, p3dgxe$beta, p3dgxe$score, p3dgxe$zb, p3dgxe$bb,
                              p3dgxe$h, p3dgxe$rtr, p3dgxe$t, p3dgxe$qr, p3dgxe$rbr, p3dgxe$logLikelihood,
                              xr1, xr2, p1d$logLikelihood, loglikelihoods, estimates, skippedFilename)
    if (scanResult > 0) {
      stop("Error scanning genes")
    }
    AppendGxEResults(outputFile, snpID, chromosome, locations,
                     refAllele, altAllele, numSub, numCases,
                     loglikelihoods, estimates, length(firstSNP:lastSNP),
                     standardizedX$sigmaE[length(standardizedX$sigmaE)])
  }

  return(0)
}

