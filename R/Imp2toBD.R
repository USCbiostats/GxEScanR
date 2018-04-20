#' Function to create a binary dosage file from an Impute 2 style format.
#'
#' Function to create a binary dosage file from an Impute 2 style format.
#'
#' @param i2info
#' Information about Impute 2 file. List returned from GetI2Info
#' @param filename
#' Base name of binary dosage file to be written.
#' The binary dosage file will have the extension .bdose.
#' If format written is less 3 a family and map file will be created with
#' the family file having the extension .fam and the map file having the
#' extenstion .map
#' @param format
#' Format of the binary dosage file to generate.
#' @param subformat
#' Subformat of the binary dosage file to generate.
#' @return
#' 0 - success, 1 - failure
#' @export
Imp2toBD <- function (i2info, filename, format = 1, subformat = 1) {
  if (format < 4) {
    famDF <- i2info$subjects
    if (i2info$usesFID == FALSE)
      famDF$FID <- seq(1,i2info$numSubjects)
    famDF$PID <- 0
    famDF$MID <- 0
    famDF$sex <- 9
    famDF$pheno <- 9
    famDF <- famDF[,c("FID", "IID", "PID", "MID", "sex", "pheno")]
    famFile <- paste(filename, ".fam", sep = "")
    write.table(famDF, famFile, quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)
    
    mapDF <- i2info$snps
    if (i2info$snpCol[2] == 0)
      mapDF$CHR <- "#"
    if (i2info$snpCol[3] == 0)
      mapDF$BP <- 0
    if (i2info$snpCol[4] == 0)
      mapDF$A1 <- '1'
    if (i2info$snpCol[5] == 0)
      mapDF$A2 <- '2'
    mapDF$cm <- 0
    mapDF <- mapDF[,c("CHR", "SNP", "BP", "cm", "A1", "A2")]
    mapFile <- paste(filename, ".map", sep = "")
    write.table(mapDF, mapFile, quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)
  }
  return (0)
}
  