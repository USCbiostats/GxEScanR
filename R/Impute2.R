#' Function to get information about a file in an Impute 2 format
#'
#' Function to get information about a file in an Impute 2 like format.
#' An Impute 2 like format has a group of columns with SNP information
#' followed by a group of dosage or genetic probabilities for the subjects.
#' There must be a column with the SNP name that uniquely identifies the SNP.
#' Optional column values can include chromosome, location in base pairs,
#' reference and alternate alleles. The file may contain a header. If there
#' is a header, the columns with SNP data can be found using the column names
#' or by specifying the column numbers. If there is a header, the subject and
#' family IDs are read from the header. If there is not header, the subject
#' and family IDs must be specified.
#'
#' @param i2file
#' Impute 2 file name. Must be provided.
#' @param header
#' Does the Impute 2 file have a header? Default value is TRUE.
#' @param snpCol
#' Must be provided if there is no header.
#' If this case this must be a integer vector of the column numbers for SNP name,
#' chromosome, location in BP, reference allele, and alternate allele respectively.
#' All must be provided a value of 0 indicates the value isn't in the file. SNP name
#' column is required.
#' If provided and there is a header, the vector may integer vector as described
#' above or it may be a character array with the column names of the appropriate
#' values. Unused columns can be indicated the an empty string or NA. If snpCol
#' is not provided and there is a header, snpCol takes on the default values of
#' 'SNP', 'CHR', 'BP', 'A1' and 'A2'.
#' @param subID
#' A data frame with two columns. Each column must be a character. The first column
#' is the family ID and the second column is the subject ID. If family IDs are not
#' used the family ID column must be filled with empty strings and not NAs.
#' @param startcol
#' Column number where dosages or genetic probabilities begin. Must be provided.
#' @param format
#' Number of genetic probabilities in file, must be 1, 2 or 3
#' 1 indicates that only the dosage is in the file Default value is 2.
#' @param usesFID
#' Does the impute file contain both the family and subject ID? Not used if
#' there is not header. Default value is TRUE.
#' @param sep
#' Separator used in file. Currently only space and tabs are allowed.
#' Default value is tab.
#' @return
#' List with information needed by GxEScan routine to read SNP data from file
#' @export
GetI2Info <- function(i2file, header = TRUE, snpCol = c("SNP", "CHR", "BP", "A1", "A2"), subID, startCol, format = 2L, usesFID = TRUE, sep = '\t') {
  fileinfo <- list(filetype = "Impute2",
                   filename = "",
                   header = header,
                   snpCol = rep(0L,5),
                   startCol = 0,
                   numSubjects = 0,
                   usesFID = usesFID,
                   subjects = data.frame(FID = character(), IID = character(), stringsAsFactors = FALSE),
                   format = format,
                   numSNPs = 0,
                   snps = data.frame(SNP = character(), CHR = character(), BP = integer(), A1 = character(), A2 = character(), stringsAsFactors = FALSE),
                   sep = sep)
  # Is a file name provided?
  if (missing(i2file) == TRUE)
    stop("No file name specified")
  # Is file name a string
  if (is.character(i2file) == FALSE)
    stop("File name is not a string")
  fileinfo$filename <- i2file
  
  # Is the header a logical value?
  if (is.logical(header) == FALSE)
    stop("Value for header is not a logical value")
  # If header is TRUE then subID must not be present
  if (header == TRUE & missing(subID) == FALSE)
    stop("Cannot specify subject IDs if header is TRUE")
  # If header is FALSE then subID must be specified
  if (header == FALSE & missing(subID) == TRUE)
    stop("If header is FALSE then subject IDs must be specified")
  
  # Is usesFID a logical value?
  if (is.logical(usesFID) == FALSE)
    stop("Value for usesFID is not a logical value")
  # Verify subID is good - has 2 columns and they are character values
  if (missing(subID) == FALSE) {
    if (is.data.frame(subID) == FALSE)
      stop("subID is not a data frame")
    if (ncol(subID) != 2)
      stop("subID must have two columns")
    if(all(sapply(subID, typeof) == "character") == FALSE)
      stop("columns in subID must be character values")
    numSub = nrow(subID)
    fileinfo$numSubjects <- numSub
    fileinfo$subjects <- subID
    colnames(fileinfo$subjects) <- c("FID", "IID")
  }
  
  # Is a genetic data starting column number provided?
  if (missing(startCol) == TRUE)
    stop("No data start column specified")
  fileinfo$startCol <- startCol
  # The starting column value must be an integer
  if ((startCol == as.integer(startCol)) == FALSE)
    stop("Starting data column must be an integer")
  startCol <- as.integer(startCol)
  # Starting column number must be greater than 1 - 1 column is needed to identify the SNP
  if (startCol < 2)
    stop("Starting data column must be greater than 1")
  
  # Is snpCol a vector of length 5?
  if (is.vector(snpCol) == FALSE)
    stop("snpCol is not a vector")
  if (length(snpCol) != 5)
    stop("Length of snpCol vector must be 5")
  # If there is a header the snpCol vector may be character string
  # in which case no checks are need. Otherwise checks that it is
  # an integer string are required
  if (all(sapply(snpCol, class) == "character") == TRUE) {
    if (header == FALSE)
      stop("If snpCol is a character string, there must be a header")
    # SNP column names must be unique
    snpColNonBlank <- snpCol[snpCol != '']
    if (length(snpColNonBlank) != length(unique(snpColNonBlank)))
      stop("Non blank SNP column names must be unique")
    if (length(snpColNonBlank) >= startCol)
      stop("Number of column names provided for SNP data is equal to or larger than the data start column number")
    # A column name for the SNP name must be provided
    if (snpCol[1] == '')
      stop("A column name must be provided the SNP name")
    snpChrString <- TRUE
  } else {
    # Column numbers must be integers
    if (all(snpCol == as.integer(snpCol)) == FALSE)
      stop("SNP column numbers are not integers or characters")
    snpCol <- as.integer(snpCol)
    fileinfo$snpCol <- snpCol
    # Column numbers must be nonnegative
    if (min(snpCol) < 0)
      stop("SNP column numbers cannot be negative")
    # SNP name column must exist, id est, > 0
    if (min(snpCol[1]) < 1)
      stop("Column for SNP name must be positive")
    # SNP column numbers must be unique
    snpColNonZeros <- snpCol[snpCol > 0]
    if (length(snpColNonZeros) != length(unique(snpColNonZeros)))
      stop("Nonzero SNP column numbers must be unique")
    # Is the starting data column number greater than the largest
    # column number in  the snp column array?
    if (startCol <= max(snpCol))
      stop("Data start column must be greater than largest SNP column value")
    snpChrString <- FALSE
  }
  
  # Is format an integer from 1 to 3?
  if ((as.integer(format) == format) == FALSE)
    stop("Value for format must be an integer from 1 to 3")
  format <- as.integer(format)
  if (format < 1 | format > 3)
    stop("Value for format must be an integer from 1 to 3")
  # Is sep a character and an allowed value
  if (is.character(sep) == FALSE)
    stop("sep is not a character")
  if ((sep %in% c('\t', ' ', ',')) == FALSE)
    stop ("Invalid value for sep. Must be ' ', '\\t' or ','")
  
  if (header == TRUE) {
    # Read the header
    con <- file(description = i2file, open = "r")
    headerLine <- scan(file = con,  nlines = 1, quiet = TRUE, what = 'character', sep = sep)
    close(con)
    # Are there enough columns?
    if (length(headerLine) < startCol)
      stop("Number of columns in file is less than the genetic data starting column")
    # Get the values that can be SNP data columns
    snpDataHeader <- headerLine[1:(startCol - 1)]
    # Are the column names unique?
    if (length(snpDataHeader) != length(unique(snpDataHeader)))
      stop("Column names in file for SNP data are not unique")
    # Find the columns with the appropriate names
    if (snpChrString == TRUE) {
      snpColUsed <- match(snpColNonBlank, snpDataHeader)
      if (length(snpColUsed[is.na(snpColUsed)]) > 0)
        stop("Not all column names found")
      isnpCol <- rep(0L, 5)
      isnpCol[snpCol != ''] <- snpColUsed
      fileinfo$snpCol <- isnpCol
      snpCol <- isnpCol
    }
    # Get the subject IDs
    if (usesFID == FALSE) {
      df <- data.frame(FID = rep("", length(headerLine) - startCol + 1), IID = headerLine[startCol:length(headerLine)], stringsAsFactors = FALSE)
    } else {
      if (((length(headerLine) - startCol) %% 2) != 1)
        stop("Odd number of entries for family and subject IDs")
      df <- data.frame(FID = headerLine[seq(startCol, length(headerLine) - 1, 2)], IID = headerLine[seq(startCol + 1, length(headerLine), 2)], stringsAsFactors = FALSE)
    }
    fileinfo$subjects <- df
    numSub <- nrow(df)
    fileinfo$numSubjects <- numSub
    # Read the first line of data
    con <- file(description = i2file, open = "r")
    firstDataLine <- scan(file = con,  skip = 1, nlines = 1, quiet = TRUE, what = character(), sep = sep)
    close(con)
    if (length(firstDataLine) != startCol + format * numSub - 1)
      stop("Number of entries on first data file does not agree with number of subjects and format")
  } else {
    # Read in the first line of data
    con <- file(description = i2file, open = "r")
    firstDataLine <- scan(file = con,  nlines = 1, quiet = TRUE, what = character(), sep = sep)
    close(con)
    if (length(firstDataLine) != startCol + format * numSub - 1)
      stop("Number of entries on first data file does not agree with number of subjects and format")
  }
  
  # Read in the information about the SNPs
  numCol <- fileinfo$numSubjects * format + startCol - 1
  colClasses <- c(rep("NULL", numCol))
  snpColNonZeros <- snpCol[snpCol > 0]
  colClasses[snpColNonZeros] <- "character"
  if (snpCol[3] > 0)
    colClasses[snpCol[3]] <- "integer"
  snpColNames <- c("SNP", "CHR", "BP", "A1", "A2")
  if (header == TRUE)
    snps <- read.table(i2file, sep = sep, skip = 1, colClasses = colClasses, quote = "")
  else
    snps <- read.table(i2file, sep = sep, header = FALSE, colClasses = colClasses, quote = "")
  colnames(snps) <- snpColNames[snpColNonZeros]
  fileinfo$snps <- snps
  fileinfo$numSNPs <- nrow(snps)
  
  return (fileinfo)
}
