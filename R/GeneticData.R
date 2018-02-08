FamilyFileCheck <- function(familyFile) {
  df <- read.table(familyFile, stringsAsFactors = FALSE)
  if (ncol(df) != 6) {
    print("Family file does not have 6 columns")
    return (0);
  }
  return (df[,c(1:2)])
}

MapFileCheck <- function(mapFile) {
  df <- read.table(mapFile, stringsAsFactors = FALSE)
  if (ncol(df) != 6) {
    print("Map file does not have 6 columns")
    return (0);
  }
  if (is.numeric(df[,3]) == FALSE | is.numeric(df[,4]) == FALSE) {
    print("Third or fourth column of map file is not numeric")
    return (0)
  }
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
#' @return
#' List with information about genetic data need by the GxEScan routine
#' @export
GetGeneticInfo <- function(geneticFile, familyFile, mapFile) {
  format <- 0

  if (missing(geneticFile))
    return (0)
  if (missing(familyFile)) {
    if (missing(mapFile)) {
      format <- 4
    } else {
      return (0)
    }
  } else {
    if (missing(mapFile)) {
      return (0)
    }
  }

  if (format == 4) {
    res <- GetBinaryDosageInformation(geneticFile, 0, 0)
    return (res)
  }

  famdf <- FamilyFileCheck(familyFile)
  if (is.data.frame(famdf) == FALSE)
    return (0)
  mapdf <- MapFileCheck(mapFile)
  if (is.data.frame(mapdf) == FALSE)
    return (0)
  res <- GetBinaryDosageInformation(geneticFile, nrow(famdf), nrow(mapdf))
  if (res$success == FALSE)
    return (0)
  if (res$format == 4)
    return (res)
  res$subjects$useFID <- TRUE
  res$subjects$FID <- famdf[,1]
  res$subjects$IID <- famdf[,2]
  info <- mapdf[,c(2,1,4:6)]
  colnames(info) <- c("ID", "Chromosome", "bp", "refAllele", "altAllele")
  snps <- list(info)
  res$SNPs <- snps;
  return (res)
}