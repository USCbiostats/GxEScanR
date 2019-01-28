#' Function to fit the D|G,E and D|G,E,GxE models and return the results
#' for the G, GxE, and the two degree of freedom test for both
#' 
#' Function to fit the D|G,E and D|G,E,GxE models and return the results
#' for the G, GxE, and the two degree of freedom test for both
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
#' @param minMAF1
#' Minimum minor allele frequency in complete data set. Has to be a value between 0.0001 and 0.25.
#' Default value is 0.05.
#' @param minMAF2
#' Minimum minor allele frequency in subjects with complete data. Has to be a value between 0.0001 and 0.25.
#' Default value is 0.05.
#' @param snps
#' Subset of SNPs to use. This can be a character vector of SNP names or an
#' integer vector of locations in the snps data frame in the genetic data list
#' @return
#' 0 - success
#' 1 - failure
#' @export
GxEOnly <- function(subjectData, geneticData, outputFile, minMAF1 = 0.05, minMAF2 = 0.05, snps) {
  # This part of the routine removes subjects without genetic data or missing covariate data
  # It also maps the SNPs in the genetic data file with the subjects in data file
  ###########################################################
  subjectDataToUse <- subjectData[complete.cases(subjectData),]
  m1 <- match(geneticData$Samples$SID, subjectDataToUse[,1])
  subjectDataToUse <- subjectDataToUse[m1[complete.cases(m1)],]
  SubjectToGene <- match(subjectDataToUse[,1], geneticData$Samples$SID)
  rm(m1)
  
  # This part of the routine preps the covariate and outcome data and fits the model without the gene
  ###########################################################

  # Standard deviation of the covariate in the interaction.
  # Needed to adjust the beta estimate of the GxE term.
  sigmaE = sqrt(var(subjectDataToUse[,ncol(subjectData)]))
  # Standardize the covariates
  x <- scale(as.matrix(subjectDataToUse[,3:ncol(subjectData)]))
  # Get the outcome
  y <- subjectDataToUse[,2]
  
  # Create a data frame with the outcome and covariate data
  rlrdf <- data.frame(y, x)
  
  # Fit the model using R's glm routine and save the log likelihood
  rlr <- glm(y ~ ., data = rlrdf, family = "binomial")
  
  # Fit logistic regression using C++ logReg routine.
  # This will start at the values glm converged to.
  ###########################################################
  # Add the column for the intercept term
  x <- cbind(1., x)
  # Allocate the memory needed to fit the model with no genes
  p1 <- GxEScanR:::AllocateLRModConstantMemory(x, y)
  p2 <- GxEScanR:::AllocateLRModFixedMemory(p1$n, p1$p)
  # Assign the betas from the glm results and assign the values
  # needed to fit the genes efficiently
  p1$beta <- rlr$coefficients
  result <- GxEScanR:::InitializeLRMod(p1$n, p1$p, y, x, p1$beta, p1$score, p1$w, p1$wInv,
                                       p1$yp, p1$zt, p1$k, p1$ql, p1$rtl,
                                       p2$abx, p2$expabx, p2$expabxp1, p2$expitabx, p1$logLikelihood)
  if (result == 1)
    return (1)
  # Clear unused memory
  rm(rlrdf)
  rm(rlr)
  
  # Prepare to read in the genes
  ###########################################################
  # Find SNPs that have sufficient minor all frequency
  snpminmaf <- match(geneticData$SNPs$SNPID[abs(geneticData$SNPInfo$AAF - 0.5) < 0.5 - minMAF1], geneticData$SNPs$SNPID)
  if (missing(snps) == FALSE)
    snpminmaf <- snpminmaf[snpminmaf %in% snps]
  # Allocate memory for reading the binary dosage file
  readInfo <- GxEScanR:::ReadBDInfo(geneticData)
  
  # Allocate space for 100 SNPs
  snpSet <- matrix(data = 0., nrow = length(SubjectToGene), ncol = 400)
  
  # Allocate memory to run G and GxE models
  p3g <- GxEScanR:::AllocateLRModNotFixedMemory(p1$n, p1$p, 1)
  p3gxe <- GxEScanR:::AllocateLRModNotFixedMemory(p1$n, p1$p, 2)
  # Memory needed for adding columns to previous LR model fit
  xr1 <- matrix <- matrix(data = 0., nrow = p1$n, ncol = 1)
  xr2 <- matrix <- matrix(data = 0., nrow = p1$n, ncol = 2)
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
  numSub = length(y)
  numCases = sum(y)
  
  # Count the number of groups and loop over them
  numGroups = ceiling(length(snpminmaf) / 100)
  # Loop over the groups of SNPs
  GxEScanR:::OpenGxEOutFile(outputFile)
  for (i in 1:numGroups) {
    firstSNP <- (i - 1) * 100 + 1
    lastSNP <- min(firstSNP + 99, length(snpminmaf))
    # print(lastSNP)
    # Read in 100 SNPs
    readResult <- GxEScanR:::ReadSNP(snpminmaf[firstSNP:lastSNP], SubjectToGene, readInfo$filename, readInfo$format, readInfo$numSub, readInfo$numSNPs,
                                     readInfo$bufferSize, readInfo$buffer, readInfo$sections, readInfo$snpSection,
                                     readInfo$fileLocation, readInfo$snpLocation, readInfo$currentSection,
                                     readInfo$dosage, readInfo$p0, readInfo$p1, readInfo$p2, snpSet)
    if (readResult == 1) {
      print("Error reading SNP group")
      print(i)
      return (i)
    }
    snpID[1:length(firstSNP:lastSNP)] <- geneticData$SNPs$SNPID[snpminmaf[firstSNP:lastSNP]]
    chromosome[1:length(firstSNP:lastSNP)] <- geneticData$SNPs$Chromosome[snpminmaf[firstSNP:lastSNP]]
    locations[1:length(firstSNP:lastSNP)] <- geneticData$SNPs$Location[snpminmaf[firstSNP:lastSNP]]
    refAllele[1:length(firstSNP:lastSNP)] <- geneticData$SNPs$Reference[snpminmaf[firstSNP:lastSNP]]
    altAllele[1:length(firstSNP:lastSNP)] <- geneticData$SNPs$Alternate[snpminmaf[firstSNP:lastSNP]]
    # Process the 100 SNPs
    scanResult <- GxEScanR:::ScanGenes(p1$n, p1$p, y, x, snpSet, length(firstSNP:lastSNP), minMAF2,
                                       p1$beta, p1$score, p1$w, p1$wInv, p1$yp,
                                       p1$zt, p1$k, p1$ql, p1$rtl,
                                       p2$abx, p2$expabx, p2$expabxp1, p2$expitabx,
                                       p2$yp, p2$zt, p2$k, p2$bt,
                                       p3g$xrw, p3g$beta, p3g$score, p3g$zb, p3g$bb,
                                       p3g$h, p3g$rtr, p3g$t, p3g$qr, p3g$rbr, p3g$logLikelihood,
                                       p3gxe$xrw, p3gxe$beta, p3gxe$score, p3gxe$zb, p3gxe$bb,
                                       p3gxe$h, p3gxe$rtr, p3gxe$t, p3gxe$qr, p3gxe$rbr, p3gxe$logLikelihood,
                                       xr1, xr2, p1$logLikelihood, loglikelihoods, estimates)
    if (scanResult > 0) {
      print("Error maximizing for SNP")
      print(snpminmaf[firstSNP + scanResult - 1])
      return (scanResult)
    }
    GxEScanR:::AppendGxEResults(outputFile, snpID, chromosome, locations, refAllele, altAllele,
                                numSub, numCases, loglikelihoods, estimates, length(firstSNP:lastSNP), sigmaE)
  }
  return (0)
}
