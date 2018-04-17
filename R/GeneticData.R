FamilyFileCheck <- function(familyFile) {
  df <- read.table(familyFile, stringsAsFactors = FALSE, colClasses = c("character", "character", "character", "character", "numeric", "numeric"), header = FALSE, sep = '\t')
  if (ncol(df) != 6)
    stop("Family file does not have 6 columns")
  colnames(df) <- c("FID", "IID")
  return (df[,c(1:2)])
}

MapFileCheck <- function(mapFile) {
  df <- read.table(mapFile, stringsAsFactors = FALSE) #, c("character", "character", "numeric", "numeric", "character", "character"), header = FALSE)
  if (ncol(df) != 6)
    stop("Map file does not have 6 columns")
  if (is.numeric(df[,3]) == FALSE | is.numeric(df[,4]) == FALSE)
    stop("Third or fourth column of map file is not numeric")
  df <- df[,c(2,1,4:6)]
  colnames(df) <- c("SNP", "CHR", "BP", "A1", "A2")
  return (df)
}

#' Function read information about genetic data needed for GxEScan
#' 
#' Function read information about genetic data needed for GxEScan
#' This information is then passed to GxEScan
#' 
#' @param geneticFile
#' Name of file with genetic data - currently only binary dosage
#' files are supported.
#' @param familyFile
#' Name of file with family data corresponding to the genetic data
#' file. This must be in plink format.
#' @param mapFile
#' Name of the map file associated with the genetic data file. This
#' must be in extended plink map file format.
#' @param usesFID
#' Indicator if family IDs are used. If this is set to FALSE, the value
#' for FID is set to "". The .fam file must still contian a column of
#' family IDs. These values are ignored if usesFID is set to TRUE.
#' @return
#' List with information about genetic data need by the GxEScan routine
#' @export
GetBinaryDosageInfo <- function(geneticFile, familyFile, mapFile, usesFID = TRUE) {
  format <- 0

  if (missing(geneticFile) == TRUE)
    stop("No binary dosage file specified")
  if (missing(familyFile) == TRUE) {
    if (missing(mapFile) == TRUE) {
      format <- 4
    } else {
      stop("If map file is specified, family file must be specified")
    }
  } else {
    if (missing(mapFile)) {
      stop("If family file is specified, map file must be specified")
    }
  }

  if (format == 4) {
    res <- GetBinaryDosageInformation(geneticFile, 0, 0)
    return (res)
  }

  famdf <- FamilyFileCheck(familyFile)
  if (is.data.frame(famdf) == FALSE)
    stop("Error reading family file")
  mapdf <- MapFileCheck(mapFile)
  if (is.data.frame(mapdf) == FALSE)
    stop("Error reading map file")
  res <- GetBinaryDosageInformation(geneticFile, nrow(famdf), nrow(mapdf))
  if (res$format == 4)
    return (res)
  if (usesFID == FALSE) {
    famdf$FID = ""
    res$usesFID = FALSE
  }
  res$numSubjects <- nrow(famdf)
  res$subjects <- famdf
  res$numSNPs <- nrow(mapdf)
  res$snps <- mapdf
#  res$subjects$useFID <- TRUE
#  res$subjects$Info <- famdf[,c(1:2)]
#  colnames(res$subjects$Info) <- c("FID", "IID")
#  Info <- mapdf[,c(2,1,4:6)]
#  colnames(Info) <- c("ID", "Chromosome", "bp", "refAllele", "altAllele")
#  snps <- list(Info)
#  res$SNPs <- snps;
  return (res)
}