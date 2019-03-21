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
GEOnly <- function(subjectData, geneticData, outputFile, minMAF1 = 0.05, minMAF2 = 0.05, snps) {
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
  sigmaEcase = sqrt(var(subjectDataToUse[subjectDataToUse[,2] == 1,ncol(subjectData)]))
  sigmaEcontrol = sqrt(var(subjectDataToUse[subjectDataToUse[,2] == 0,ncol(subjectData)]))
  # Standardize the covariates
  xl <- scale(as.matrix(subjectDataToUse[,3:(ncol(subjectData) - 1)]))
  # Add the column for the intercept term
  xl <- cbind(1., xl)
  # Save the column for the interacting covariate
  xr <- scale(matrix(subjectDataToUse[,ncol(subjectData)]))
  # Repeat for cases
  xlcase <- scale(as.matrix(subjectDataToUse[subjectDataToUse[,2] == 1,3:(ncol(subjectData) - 1)]))
  xlcase <- cbind(1., xlcase)
  xrcase <- scale(matrix(subjectDataToUse[subjectDataToUse[,2] == 1, ncol(subjectData)]))
  # Repeat for controls
  xlcontrol <- scale(as.matrix(subjectDataToUse[subjectDataToUse[,2] == 0,3:(ncol(subjectData) - 1)]))
  xlcontrol <- cbind(1., xlcontrol)
  xrcontrol <- scale(matrix(subjectDataToUse[subjectDataToUse[,2] == 0, ncol(subjectData)]))

  # Get the outcome
  y <- subjectDataToUse[,2]

  # Allocate the memory
  # Linear regression
  geLinReg <- GxEScanR:::GELinRegMemory(nrow(xl), ncol(xl), ncol(xr))
  caseLinReg <- GxEScanR:::GELinRegMemory(nrow(xlcase), ncol(xlcase), ncol(xrcase))
  controlLinReg <- GxEScanR:::GELinRegMemory(nrow(xlcontrol), ncol(xlcontrol), ncol(xrcontrol))
  # Logistic regression
  geLogReg <- GxEScanR:::GELogRegMemory(nrow(xl), ncol(xl), ncol(xr))
  caseLogReg <- GxEScanR:::GELogRegMemory(nrow(xlcase), ncol(xlcase), ncol(xrcase))
  controlLogReg <- GxEScanR:::GELogRegMemory(nrow(xlcontrol), ncol(xlcontrol), ncol(xrcontrol))
  
  # Initialize memory  
  # Linear regression
  GxEScanR:::IntializeGELinReg(xl, xr, geLinReg$qL, geLinReg$qR, geLinReg$rTL, geLinReg$rTR, geLinReg$rBR)
  GxEScanR:::IntializeGELinReg(xlcase, xrcase, caseLinReg$qL, caseLinReg$qR, caseLinReg$rTL, caseLinReg$rTR, caseLinReg$rBR)
  GxEScanR:::IntializeGELinReg(xlcontrol, xrcontrol, controlLinReg$qL, controlLinReg$qR, controlLinReg$rTL, controlLinReg$rTR, controlLinReg$rBR)
  # Logistic regression - not needed
  
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
  
  # Memory for results
  loglikelihoods <- matrix(data = 0., nrow = 100, ncol = 6)
  estimates <- matrix(data = 0., nrow = 100, ncol = 6)
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
  GxEScanR:::OpenGEOutFile(outputFile)
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
    for (j in 1:length(firstSNP:lastSNP)) {
      g <- snpSet[,4 * j - 3]
      if (abs(0.5 - mean(g) * 0.5) > 0.5 - minMAF2) {
        loglikelihoods[j, 1] <- NA
        estimates[j, 1] <- 0.;
      } else {
        GxEScanR:::GELinReg(g, xl, xr,
                            geLinReg$betaT0, geLinReg$betaT, geLinReg$betaB,
                            geLinReg$qL, geLinReg$qR,
                            geLinReg$rTL, geLinReg$rTR, geLinReg$rBR,
                            geLinReg$sigma2, geLinReg$logLikelihood)
        loglikelihoods[j, 1] <- 2 * (geLinReg$logLikelihood[2] - geLinReg$logLikelihood[1])
        estimates[j, 1] <- geLinReg$betaB[1] 
        gCase <- snpSet[subjectDataToUse[,2] == 1, 4 * j - 3]
        if (abs(0.5 - mean(gCase) * 0.5) > 0.5 - minMAF2) {
          loglikelihoods[j, 2] <- -1.;
          estimates[j, 2] <- 0.;
        } else {
          GxEScanR:::GELinReg(gCase, xlcase, xrcase,
                            caseLinReg$betaT0, caseLinReg$betaT, caseLinReg$betaB,
                            caseLinReg$qL, caseLinReg$qR,
                            caseLinReg$rTL, caseLinReg$rTR, caseLinReg$rBR,
                            caseLinReg$sigma2, caseLinReg$logLikelihood)
          loglikelihoods[j, 2] <- 2 * (caseLinReg$logLikelihood[2] - caseLinReg$logLikelihood[1])
          estimates[j, 2] <- caseLinReg$betaB[1] 
        }
        gControl <- snpSet[subjectDataToUse[,2] == 0, 4 * j - 3]
        if (abs(0.5 - mean(gControl) * 0.5) > 0.5 - minMAF2) {
          loglikelihoods[j, 3] <- -1.;
          estimates[j, 3] <- 0.;
        } else {
          GxEScanR:::GELinReg(gControl, xlcontrol, xrcontrol,
                            controlLinReg$betaT0, controlLinReg$betaT, controlLinReg$betaB,
                            controlLinReg$qL, controlLinReg$qR,
                            controlLinReg$rTL, controlLinReg$rTR, controlLinReg$rBR,
                            controlLinReg$sigma2, controlLinReg$logLikelihood)
          loglikelihoods[j, 3] <- 2 * (controlLinReg$logLikelihood[2] - controlLinReg$logLikelihood[1])
          estimates[j, 3] <- controlLinReg$betaB[1] 
        }
        # G|E
        y2 <- g * 0.5;
        geLogReg$beta0[1] = log(geLinReg$betaT0[1] / (2. - geLinReg$betaT0[1]))
        geLogReg$beta0[2:ncol(xl)] = geLinReg$betaT0[2:ncol(xl)] / var(g);
#        geLogReg$beta0[2:ncol(xl)] = 0.
#        print(j)
        if (GxEScanR:::GELogReg(y2, xl, xr,
                                geLogReg$yp, geLogReg$xlw, geLogReg$xrw,
                                geLogReg$zt, geLogReg$zb, geLogReg$k, geLogReg$t, geLogReg$h,
                                geLogReg$beta0, geLogReg$beta0diff,
                                geLogReg$beta, geLogReg$betaT, geLogReg$betaB,
                                geLogReg$qL, geLogReg$qR,
                                geLogReg$rTL, geLogReg$rTR, geLogReg$rBR,
                                geLogReg$abx, geLogReg$expabx, geLogReg$expabxp1, geLogReg$expitabx,
                                geLogReg$w, geLogReg$wInv,
                                geLogReg$score, geLogReg$loglikelihood) == 0) {
          loglikelihoods[j, 4] <- 2 * (geLogReg$loglikelihood[2] - geLogReg$loglikelihood[1])
          estimates[j, 4] <- geLogReg$beta[length(geLogReg$beta)]
          # case only
          y2case <- gCase * 0.5;
          caseLogReg$beta0[1] = log(caseLinReg$betaT0[1] / (2. - caseLinReg$betaT0[1]))
          caseLogReg$beta0[2:ncol(xl)] = caseLinReg$betaT0[2:ncol(xl)] / var(gCase)
          if (GxEScanR:::GELogReg(y2case, xlcase, xrcase,
                                  caseLogReg$yp, caseLogReg$xlw, caseLogReg$xrw,
                                  caseLogReg$zt, caseLogReg$zb, caseLogReg$k, caseLogReg$t, caseLogReg$h,
                                  caseLogReg$beta0, caseLogReg$beta0diff,
                                  caseLogReg$beta, caseLogReg$betaT, caseLogReg$betaB,
                                  caseLogReg$qL, caseLogReg$qR,
                                  caseLogReg$rTL, caseLogReg$rTR, caseLogReg$rBR,
                                  caseLogReg$abx, caseLogReg$expabx, caseLogReg$expabxp1, caseLogReg$expitabx,
                                  caseLogReg$w, caseLogReg$wInv,
                                  caseLogReg$score, caseLogReg$loglikelihood) == 1) {
            loglikelihoods[j, 5] <- NA
            estimates[j, 5] <- NA
            
          } else {
            loglikelihoods[j, 5] <- 2 * (caseLogReg$loglikelihood[2] - caseLogReg$loglikelihood[1])
            estimates[j, 5] <- caseLogReg$beta[length(caseLogReg$beta)]
          }
          # control only
          y2control<- gControl * 0.5;
          controlLogReg$beta0[1] = log(controlLinReg$betaT0[1] / (2. - controlLinReg$betaT0[1]))
          controlLogReg$beta0[2:ncol(xl)] = controlLinReg$betaT0[2:ncol(xl)] / var(gControl)
          if (GxEScanR:::GELogReg(y2control, xlcontrol, xrcontrol,
                                  controlLogReg$yp, controlLogReg$xlw, controlLogReg$xrw,
                                  controlLogReg$zt, controlLogReg$zb, controlLogReg$k, controlLogReg$t, controlLogReg$h,
                                  controlLogReg$beta0, controlLogReg$beta0diff,
                                  controlLogReg$beta, controlLogReg$betaT, controlLogReg$betaB,
                                  controlLogReg$qL, controlLogReg$qR,
                                  controlLogReg$rTL, controlLogReg$rTR, controlLogReg$rBR,
                                  controlLogReg$abx, controlLogReg$expabx, controlLogReg$expabxp1, controlLogReg$expitabx,
                                  controlLogReg$w, controlLogReg$wInv,
                                  controlLogReg$score, controlLogReg$loglikelihood) == 1) {
            loglikelihoods[j, 6] <- NA
            estimates[j, 6] <- NA
          } else {
            loglikelihoods[j, 6] <- 2 * (controlLogReg$loglikelihood[2] - controlLogReg$loglikelihood[1])
            estimates[j, 6] <- controlLogReg$beta[length(controlLogReg$beta)]
          }
        } else {
          loglikelihoods[j, 4] <- NA
          estimates[j, 4] <- NA
          loglikelihoods[j, 5] <- NA
          estimates[j, 5] <- NA
          loglikelihoods[j, 6] <- NA
          estimates[j, 6] <- NA
        }
      }
    }

    GxEScanR:::AppendGEResults(outputFile, snpID, chromosome, locations, refAllele, altAllele,
                              numSub, numCases, loglikelihoods, estimates, length(firstSNP:lastSNP),
                              sigmaE, sigmaEcase, sigmaEcontrol);
  }
  return (0)
}
