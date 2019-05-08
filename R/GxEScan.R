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
StandardizeCovariates <- function(subjectInfo, binCov) {
  sigmaE = sqrt(diag(var(subjectInfo$covariates)))
  if (min(sigmaE) == 0.)
    stop("At least one of the covariates is constant")
  if (binCov == TRUE) {
    uniqueCovs <- sort(unique(subjectInfo$covariates[,ncol(subjectInfo$covariates)]))
    if (length(uniqueCovs) > 2)
      stop("There are more than 2 values for the interacting binary covariate")
    if(all(uniqueCovs == c(0., 1.)) != TRUE)
      stop("Binary covariate values must be coded as 0, 1")
  }
  x <- scale(subjectInfo$covariates)
  x <- matrix(x, nrow(subjectInfo$covariates), ncol(subjectInfo$covariates))
  x <- cbind(1., x)
  eg <- as.matrix(x[,1:(ncol(x) - 1)])
  cases <- as.matrix(eg[subjectInfo$phenotype == 1,])
  controls <- as.matrix(eg[subjectInfo$phenotype == 0,])
  return (list(sigmaE = sigmaE,
               x = x,
               eg = eg,
               cases = cases,
               controls = controls))
}

# Get the outcome for the various models
GetOutcomes <- function(subjectInfo) {
  return (list(y = subjectInfo$phenotype,
               eg = subjectInfo$covariates[, ncol(subjectInfo$covariates)],
               cases = subjectInfo$covariates[subjectInfo$phenotype == 1, ncol(subjectInfo$covariates)],
               controls = subjectInfo$covariates[subjectInfo$phenotype == 0, ncol(subjectInfo$covariates)]))
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
#' @param outFile
#' Name of file to write the results to
#' @param skipFile
#' Name of file to write info about SNPs that were skipped. If this is blank
#' no file is written. If this is the same as outFile, the skipped SNPs
#' are written to the output file.
#' @param popminMaf
#' Population minimum allele frequency
#' @param sampleminMaf
#' Sample minimum allele frequency
#' @param binCov
#' Indicator if covariate interacting with the gene is a binary value. The
#' default value is TRUE.
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
#' outFile = "stdout"
#' skipFile = ""
#' popminMaf = 0.05
#' sampleminMaf = 0.1
#' binCov = TRUE
#' snps = 1L:110L
#' GxEScan(subjectData, geneticData, outFile, skipFile, popminMaf, sampleminMaf, binCov, snps)
GxEScan <- function(subjectData,
                    geneticData,
                    outFile,
                    skipFile = "",
                    popminMaf = 0.01,
                    sampleminMaf = 0.05,
                    binCov = TRUE,
                    snps) {
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
  
  if (missing(outFile) == TRUE)
    stop("No output file specified")
  if (is.character(outFile) == FALSE)
    stop("outFile must be of type character")

  if (is.character(skipFile) == FALSE)
    stop("skipFile must be of type character")

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

  if (is.logical(binCov) == FALSE)
    stop("binaryCovariate must be a logical value")
  
  if (missing(snps) == TRUE) {
    snpIndices <- seq(1L:geneticData$NumSNPs)
  } else {
    snpIndices <- GetSNPIndices(snps, geneticData$SNPs, geneticData$NumSNPs)
  }
  subjectInfo <- MatchSubjectsAndGeneticData(subjectInfo, geneticData)
  
  # Get the x's and y's for the regressions - standardizing the x's
  standardizedX <- StandardizeCovariates(subjectInfo, binCov)
  outcomes <- GetOutcomes(subjectInfo)
  # indices for cases and controls
  indicesCases <- which(outcomes$y == 1)
  indicesControls <- which(outcomes$y == 0)
  # Number of subjects and number of cases
  numSub <- length(outcomes$y)
  numCases <- sum(outcomes$y)
  numControls <- numSub - numCases
  
  # Select SNPs that meet minimum maf requirements  
  snpminmaf <- match(geneticData$SNPs$SNPID[abs(geneticData$SNPInfo$AAF - 0.5) < 0.5 - popminMaf], geneticData$SNPs$SNPID)
  snpminmaf <- snpminmaf[snpminmaf %in% snpIndices]
  bdInfo <- GetBDFileInfo(geneticData)
  
  # Memory allocation for models
  gxeMem <- AllocateLargeScaleLogRegMemory(outcomes$y, standardizedX$x, TRUE)
  if (binCov == TRUE) {
    egMem <- AllocateLargeScaleLogRegMemory(outcomes$eg, standardizedX$eg, FALSE)
    casesMem <- AllocateLargeScaleLogRegMemory(outcomes$cases, standardizedX$cases, FALSE)
    controlsMem <- AllocateLargeScaleLogRegMemory(outcomes$controls, standardizedX$controls, FALSE)
  } else {
    egMem <- AllocateLargeScaleLinRegMemory(outcomes$eg, standardizedX$eg, "E|G")
    casesMem <- AllocateLargeScaleLinRegMemory(outcomes$cases, standardizedX$cases, "case-only")
    controlsMem <- AllocateLargeScaleLinRegMemory(outcomes$controls, standardizedX$controls, "control-only")
  }

  # Memory needed for adding columns to previous LR model fit
  xr1 <- matrix(data = 0., nrow = gxeMem$p1$n, ncol = 1)
  xr2 <- matrix(data = 0., nrow = gxeMem$p1$n, ncol = 2)
  xrCases <- matrix(data = 0., nrow = numCases, ncol = 1)
  xrControls <- matrix(data = 0., nrow = numControls, ncol = 1)

  # Memory for results
  loglikelihoods <- matrix(data = 0., nrow = 100, ncol = 6)
  estimates <- matrix(data = 0., nrow = 100, ncol = 6)

  # Memory for SNPs
  snpSet <- matrix(data = 0., nrow = length(subjectInfo$indices), ncol = 400)
  snpCases <- matrix(data = 0., nrow = length(indicesCases), ncol = 400)
  snpControls <- matrix(data = 0., nrow = length(indicesControls), ncol = 400)
  # Memory for SNP info
  snpID <- character(100)
  chromosome <- character(100)
  locations <- integer(100)
  refAllele <- character(100)
  altAllele <- character(100)

  # Count the number of groups and loop over them
  numGroups = ceiling(length(snpminmaf) / 100)
  # Loop over the groups of SNPs
  OpenGxEOutFile(outFile)
  
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
    snpCases[1:nrow(snpCases), 1:400] = snpSet[indicesCases, 1:400]
    snpControls[1:nrow(snpControls), 1:400] = snpSet[indicesControls, 1:400]
    
    snpID[1:length(firstSNP:lastSNP)] <- geneticData$SNPs$SNPID[snpminmaf[firstSNP:lastSNP]]
    chromosome[1:length(firstSNP:lastSNP)] <- geneticData$SNPs$Chromosome[snpminmaf[firstSNP:lastSNP]]
    locations[1:length(firstSNP:lastSNP)] <- geneticData$SNPs$Location[snpminmaf[firstSNP:lastSNP]]
    refAllele[1:length(firstSNP:lastSNP)] <- geneticData$SNPs$Reference[snpminmaf[firstSNP:lastSNP]]
    altAllele[1:length(firstSNP:lastSNP)] <- geneticData$SNPs$Alternate[snpminmaf[firstSNP:lastSNP]]
    loglikelihoods[,] = NA
    estimates[,] = NA
    scanResult <- ScanDisease(gxeMem$p1$n, gxeMem$p1$p,
                              outcomes$y, standardizedX$x,
                              snpSet, snpID, length(firstSNP:lastSNP),
                              sampleminMaf,
                              gxeMem$p1$beta, gxeMem$p1$score, gxeMem$p1$w,
                              gxeMem$p1$wInv, gxeMem$p1$yp,
                              gxeMem$p1$zt, gxeMem$p1$k, gxeMem$p1$ql, gxeMem$p1$rtl,
                              gxeMem$p2$abx, gxeMem$p2$expabx, gxeMem$p2$expabxp1, gxeMem$p2$expitabx,
                              gxeMem$p2$yp, gxeMem$p2$zt, gxeMem$p2$k, gxeMem$p2$bt,
                              gxeMem$p3$xrw, gxeMem$p3$beta, gxeMem$p3$score, gxeMem$p3$zb, gxeMem$p3$bb,
                              gxeMem$p3$h, gxeMem$p3$rtr, gxeMem$p3$t, gxeMem$p3$qr, gxeMem$p3$rbr, gxeMem$p3$logLikelihood,
                              gxeMem$p3gxe$xrw, gxeMem$p3gxe$beta, gxeMem$p3gxe$score, gxeMem$p3gxe$zb, gxeMem$p3gxe$bb,
                              gxeMem$p3gxe$h, gxeMem$p3gxe$rtr, gxeMem$p3gxe$t, gxeMem$p3gxe$qr,
                              gxeMem$p3gxe$rbr, gxeMem$p3gxe$logLikelihood,
                              xr1, xr2, gxeMem$p1$logLikelihood, loglikelihoods,
                              estimates, skipFile)
    if (scanResult > 0)
      stop("Error scanning genes")
    estimates[,2] <- estimates[,2] / standardizedX$sigmaE[length(standardizedX$sigmaE)]
    if (binCov == TRUE) {
      scanResult <- ScanBinaryE(egMem$p1$n, egMem$p1$p,
                                outcomes$eg, standardizedX$eg,
                                snpSet, snpID, length(firstSNP:lastSNP),
                                sampleminMaf,
                                egMem$p1$beta, egMem$p1$score, egMem$p1$w,
                                egMem$p1$wInv, egMem$p1$yp,
                                egMem$p1$zt, egMem$p1$k, egMem$p1$ql, egMem$p1$rtl,
                                egMem$p2$abx, egMem$p2$expabx, egMem$p2$expabxp1, egMem$p2$expitabx,
                                egMem$p2$yp, egMem$p2$zt, egMem$p2$k, egMem$p2$bt,
                                egMem$p3$xrw, egMem$p3$beta, egMem$p3$score, egMem$p3$zb, egMem$p3$bb,
                                egMem$p3$h, egMem$p3$rtr, egMem$p3$t, egMem$p3$qr, egMem$p3$rbr, egMem$p3$logLikelihood,
                                xr1, egMem$p1$logLikelihood, loglikelihoods,
                                estimates, 3L, skipFile, "G|E")
      if (scanResult > 0)
        stop("Error scanning genes")
      scanResult <- ScanBinaryE(casesMem$p1$n, casesMem$p1$p,
                                outcomes$cases, standardizedX$cases,
                                snpCases, snpID, length(firstSNP:lastSNP),
                                sampleminMaf,
                                casesMem$p1$beta, casesMem$p1$score, casesMem$p1$w,
                                casesMem$p1$wInv, casesMem$p1$yp,
                                casesMem$p1$zt, casesMem$p1$k, casesMem$p1$ql, casesMem$p1$rtl,
                                casesMem$p2$abx, casesMem$p2$expabx, casesMem$p2$expabxp1, casesMem$p2$expitabx,
                                casesMem$p2$yp, casesMem$p2$zt, casesMem$p2$k, casesMem$p2$bt,
                                casesMem$p3$xrw, casesMem$p3$beta, casesMem$p3$score, casesMem$p3$zb, casesMem$p3$bb,
                                casesMem$p3$h, casesMem$p3$rtr, casesMem$p3$t, casesMem$p3$qr, casesMem$p3$rbr, casesMem$p3$logLikelihood,
                                xrCases, casesMem$p1$logLikelihood, loglikelihoods,
                                estimates, 4L, skipFile, "case-only")
      if (scanResult > 0)
        stop("Error scanning genes")
      scanResult <- ScanBinaryE(controlsMem$p1$n, controlsMem$p1$p,
                                outcomes$controls, standardizedX$controls,
                                snpControls, snpID, length(firstSNP:lastSNP),
                                sampleminMaf,
                                controlsMem$p1$beta, controlsMem$p1$score, controlsMem$p1$w,
                                controlsMem$p1$wInv, controlsMem$p1$yp,
                                controlsMem$p1$zt, controlsMem$p1$k, controlsMem$p1$ql, controlsMem$p1$rtl,
                                controlsMem$p2$abx, controlsMem$p2$expabx, controlsMem$p2$expabxp1, controlsMem$p2$expitabx,
                                controlsMem$p2$yp, controlsMem$p2$zt, controlsMem$p2$k, controlsMem$p2$bt,
                                controlsMem$p3$xrw, controlsMem$p3$beta, controlsMem$p3$score, controlsMem$p3$zb, controlsMem$p3$bb,
                                controlsMem$p3$h, controlsMem$p3$rtr, controlsMem$p3$t, controlsMem$p3$qr, controlsMem$p3$rbr, controlsMem$p3$logLikelihood,
                                xrControls, controlsMem$p1$logLikelihood, loglikelihoods,
                                estimates, 5L, skipFile, "control-only")
      if (scanResult > 0)
        stop("Error scanning genes")
    } else {
      scanResult <- ScanContinuousE(egMem$n,
                                    egMem$p,
                                    outcomes$eg,
                                    standardizedX$eg,
                                    snpSet,
                                    snpID,
                                    length(firstSNP:lastSNP),
                                    sampleminMaf,
                                    egMem$ql,
                                    egMem$rtl,
                                    egMem$k,
                                    egMem$bt,
                                    egMem$zb,
                                    egMem$bb,
                                    egMem$h,
                                    egMem$rtr,
                                    egMem$t,
                                    egMem$qr,
                                    egMem$rbr,
                                    egMem$logLikelihood,
                                    xr1,
                                    loglikelihoods,
                                    estimates,
                                    3L, skipFile, "G|E")
      if (scanResult > 0)
        stop("Error scanning genes")
      scanResult <- ScanContinuousE(casesMem$n,
                                    casesMem$p,
                                    outcomes$cases,
                                    standardizedX$cases,
                                    snpCases,
                                    snpID,
                                    length(firstSNP:lastSNP),
                                    sampleminMaf,
                                    casesMem$ql,
                                    casesMem$rtl,
                                    casesMem$k,
                                    casesMem$bt,
                                    casesMem$zb,
                                    casesMem$bb,
                                    casesMem$h,
                                    casesMem$rtr,
                                    casesMem$t,
                                    casesMem$qr,
                                    casesMem$rbr,
                                    casesMem$logLikelihood,
                                    xrCases,
                                    loglikelihoods,
                                    estimates,
                                    4L, skipFile, "case-only")
      if (scanResult > 0)
        stop("Error scanning genes")
      scanResult <- ScanContinuousE(controlsMem$n,
                                    controlsMem$p,
                                    outcomes$controls,
                                    standardizedX$controls,
                                    snpControls,
                                    snpID,
                                    length(firstSNP:lastSNP),
                                    sampleminMaf,
                                    controlsMem$ql,
                                    controlsMem$rtl,
                                    controlsMem$k,
                                    controlsMem$bt,
                                    controlsMem$zb,
                                    controlsMem$bb,
                                    controlsMem$h,
                                    controlsMem$rtr,
                                    controlsMem$t,
                                    controlsMem$qr,
                                    controlsMem$rbr,
                                    controlsMem$logLikelihood,
                                    xrControls,
                                    loglikelihoods,
                                    estimates,
                                    5L, skipFile, "control-only")
      if (scanResult > 0)
        stop("Error scanning genes")
    }
    AppendGxEScanResults(outFile, snpID, chromosome, locations,
                         refAllele, altAllele, numSub, numCases,
                         loglikelihoods, estimates, length(firstSNP:lastSNP),
                         standardizedX$sigmaE[length(standardizedX$sigmaE)])
  }

  return(0)
}

