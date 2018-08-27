#' Function to subset SNPS
#' 
#' Function to subset SNPs by minor allele frequency and r squared
#' from a binary dosage file in format 4 or greater. If the binary
#' dosage file contains more than one group of imputation results,
#' the minimum r squared and minor allele frequency of the groups
#' are used for comparison with the values supplied.
#' 
#' @param bdInfo
#' A binary dosage info list returned from GetBinaryDosageInfo
#' This information has to be from a binary dosage file in
#' format 4 or greatef.
#' @param minMAF
#' Mimimum minor allele frequency of the SNP. A value of zero
#' indicates to not subset based on minor allele frequency.
#' Default value = 0
#' @param minRsq
#' Minimum r squared of the imputation for the SNP. A value of
#' zero indicates to not subset based on r squared value.
#' Default value = 0
#' @param useAltFreq
#' If TRUE, minor allele frequency is calculated from the alternate
#' allele frequency. If FALSE, the minor allele frequency is used.
#' Default value is TRUE.
#' @return
#' Vector of SNP indices to use in the imputation - Usually passed
#' to GxEScan as the snps value
#' @export
SubsetBDSNPs <- function(bdInfo, minMAF = 0, minRsq = 0, useAltFreq = TRUE) {
  if (missing(bdInfo) == TRUE)
    stop ("No binary dosage information provided")
  if (bdInfo$filetype != 'BinaryDosage')
    stop ("List provided is not a binary dosage list")
  if (bdInfo$format != 4)
    stop ("Binary dosage file is not in format 4 or greater")
  if (minMAF == 0 & minRsq == 0)
    stop ("A nonzero value has to be provided for either minMAF or minRsq")
  if (minMAF < 0 | minMAF > 0.5)
    stop ("Value for minMAF must be from 0 to 0.5")
  if (minRsq < 0 | minRsq > 1)
    stop ("Value for minRsq must be from 0 to 1")

  snpInfoColNames <- colnames(bdInfo$SNPinfo)
  if (bdInfo$Groups > 1) {
    if (useAltFreq == TRUE)
      mafColNames <- paste(rep("AAF"), 1:bdInfo$Groups, sep = '.')
    else
      mafColNames <- paste(rep("MAF"), 1:bdInfo$Groups, sep = '.')
    rsqColNames <- paste(rep("Rsq"), 1:bdInfo$Groups, sep = '.')
  } else {
    if (useAltFreq == TRUE)
      mafColNames <- c("AAF")
    else
      mafColNames <- c("MAF")
    rsqColNames <- c("Rsq")
  }
  
  if (minMAF != 0) {
    mafColNum <- match(mafColNames, snpInfoColNames)
    if (anyNA(mafColNum) == TRUE)
      stop ("Cannot find minor allele frequncy in bdInfo")
  }
  if (minRsq != 0) {
    rsqColNum <- match(rsqColNames, snpInfoColNames)
    if (anyNA(rsqColNum))
      stop ("Cannot find r squared value in bdInfo")
  }
  
  if (minRsq == 0) {
    if (useAltFreq)
      return (which(apply((0.5 - abs(bdInfo$SNPinfo[,mafColNum] - 0.5)), 1, 'min') > minMAF))
    return (which(apply(bdInfo$SNPinfo[,mafColNum], 1, 'min') > minMAF))
  }
  if (minMAF == 0)
    return (which(apply(bdInfo$SNPinfo[,rsqColNum], 1, 'min') > minRsq))
  if (useAltFreq)
    return (which(apply(bdInfo$SNPinfo[,rsqColNum], 1, 'min') > minRsq & apply((0.5 - abs(bdInfo$SNPinfo[,mafColNum] - 0.5)), 1, 'min') > minMAF))
  return (which(apply(bdInfo$SNPinfo[,rsqColNum], 1, 'min') > minRsq & apply(bdInfo$SNPinfo[,mafColNum], 1, 'min') > minMAF))
}
  