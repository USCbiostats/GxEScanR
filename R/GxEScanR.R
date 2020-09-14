#' @useDynLib GxEScanR, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom prodlim row.match
#' @importFrom stats complete.cases glm var
NULL

subsetsnps <- function(snps, snplist) {
  if (length(snps) == 0) {
    stop("No SNPs selected")
  }
  if (is.character(snps) == TRUE) {
    if (length(snps) == 1 & snps[1] == "all")
      return(rep(TRUE, length(snplist)))
    snps2 <- match(snps, snplist)
    snps2 <- snps2[!is.na(snps2)]
    if (length(snps2) == 0)
      stop("No SNPs found")
    if (length(snps) != length(snps2))
      print(paste(length(snps) - length(snps2), " SNP(s) were not found"))
    snps <- snps2
  }
  if (is.numeric(snps) == FALSE)
    stop("snps must be character or integer array")
  if (is.integer(snps) == FALSE) {
    if (all(floor(snps) == snps) == FALSE)
      stop("snps must be character or integer array")
  }
  if (min(snps) < 1)
    stop("numeric snp values must be positive")
  if (max(snps) > length(snplist))
    stop("at least one snp value is invalid")
  snpstouse <- rep(FALSE, length(snplist))
  snpstouse[snps] <- TRUE
  
  return (snpstouse)
}
# Subsets the covariate data to complete cases
# Also gets subject indices in the genetic data
subsetdata <- function(subdata, bdinfo, binary, mincov) {
  # There must be at least two columns in the subject data
  if(ncol(subdata) < 2)
    stop("There must me at least two columns in the subject data")
  
  # Check if the first column is a character value
  if (is.character(subdata[,1]) == FALSE)
    stop("First column of subject data must be a character value")
  
  # Remove subjects without complete data
  covdata <- subdata[complete.cases(subdata),]
  colnames(covdata) <- colnames(subdata)
  if (nrow(covdata) == 0)
    stop("No subjects have complete phenotype/covariate data")
  
  # Determine if family IDs are used
  # and get the indices of the subjects in the genetic data that are
  # in the covariate data
  if (bdinfo$usesfid == FALSE) {
    covdata$genindex <- match(covdata[,1], bdinfo$samples$sid)
    phenocol <- 2
  } else {
    # When family ID is used, there must be at least 3 columns
    # in the subject data
    if (ncol(covdata) < 3)
      stop("When using family ID, subject data must have at least 3 columns")
    # When family ID is used, the second columns must be a character value
    if (is.character(covdata[,2]) == FALSE)
      stop("When using family ID, the first two columns must be character values")
    covdata$genindex <- row.match(covdata[,1:2], bdinfo$samples)
    phenocol <- 3
  }
  if (ncol(covdata) < phenocol + mincov)
    stop("Subject data has no covariates")
  for (i in phenocol:ncol(covdata)) {
    if (is.numeric(covdata[,i]) == FALSE)
      stop("Phenotype and covariate values must be numeric")
  }

  # Drop subjects that don't have genetic values
  covdata <- covdata[complete.cases(covdata),]
  # Check if any subjects have complete data
  if (nrow(covdata) == 0)
    stop("No subjects have complete data")

  # If the outcome is binary sort by outcome and
  if (binary == TRUE) {
    covdata <- covdata[order(covdata[, phenocol]),]
    # Check there are two outcome values
    if (length(unique(covdata[,phenocol])) != 2)
      stop("When using a binary outcome there must only be two values")
    # check if outcome values are 0, 1
    if (all(unique(covdata[, phenocol]) == c(0,1)) == FALSE)
      stop("When using a binary outcome must be coded 0,1")
  }
  
  phenocov = as.matrix(covdata[,phenocol:(ncol(covdata)-1)])
  dimnames(phenocov) <- list(covdata[,ncol(covdata)],
                             colnames(covdata)[phenocol:(ncol(covdata)-1)])
  return (phenocov)
}

# Create blocks used for faster reading of binary dosage files
assignblocks <- function(nsub, nsnps, snploc, snpbytes, reqblksize) {
  if (nsub < 10000) {
    blksnps <- 5000
  } else if (nsub < 25000) {
    blksnps <- 2000
  } else if (nsub < 50000) {
    blksnps <- 1000
  } else if (nsub < 100000) {
    blksnps <- 500
  } else if (nsub < 250000) {
    blksnps <- 200
  } else if (nsub < 500000) {
    blksnps <- 100
  } else {
    blksnps <- 50
  }
  
  if (reqblksize > 2 * blksnps & nsnps > 2 * blksnps)
    stop("Requested block size greater than twice the recommended block size")
  if (reqblksize != 0)
    blksnps <- reqblksize
  nblks <- ceiling(nsnps / blksnps)
  if (nblks == 1) {
    blksnps <- nsnps
    blkloc <- snploc[1]
    blkbytes <- snploc[nsnps] - snploc[1] + snpbytes[length(snpbytes)]
  } else {
    fsnp <- seq(1, (nblks - 1) * blksnps + 1, blksnps)
    blkloc <- snploc[fsnp]
    blkbytes <- numeric(nblks)
    blkbytes[1:(nblks - 1)] <- snploc[fsnp[2:nblks]] - snploc[fsnp[1:(nblks - 1)]]
    blkbytes[nblks] <- snploc[nsnps] - snploc[fsnp[nblks]] + snpbytes[length(snpbytes)]
  }
  return(list(snpsperblk = blksnps,
              blkloc = blkloc,
              blkbytes = blkbytes))
}

#####################################################
###          Check for valid input values
#####################################################
validateinput <- function(data, bdinfo, outfile, skipfile,
                          minmaf, blksize, binary) {
  # Check if input values are of correct type
  if (is.data.frame(data) == FALSE)
    stop("data must be a data frame")
  if (class(bdinfo) != "genetic-info")
    stop("bdinfo not a genetic-info class")
  if (class(bdinfo$additionalinfo) != "bdose-info")
    stop("bdinfo does not have information about a binary dosage file")
  if (is.character(outfile) == FALSE)
    stop("outfile must be a character value")
  if (length(outfile) != 1)
    stop("outfile must be a character vector of length 1")
  if (is.character(skipfile) == FALSE)
    stop("outfile must be a character value")
  if (length(skipfile) != 1)
    stop("outfile must be a character vector of length 1")
  if (is.numeric(minmaf) == FALSE)
    stop("minmaf must be a numeric value")
  if (length(minmaf) != 1)
    stop("minmaf must be a numeric vector of length 1")
  if (minmaf < 0 | minmaf > 0.25)
    stop("minmaf must be a value from 0 to 0.25, inclusive")
  if (is.numeric(blksize) == FALSE)
    stop("blksize must be an integer")
  if (length(blksize) != 1)
    stop("blksize must be an integer vector of length 1")
  if (blksize != floor(blksize))
    stop("blksize must be an integer")
  blksize <- as.integer(blksize)
  if (blksize < 0)
    stop("blksize must be greater than or equal to 0")
  if (is.logical(binary) == FALSE)
    stop("binary must be a logical value")
  if (length(binary) != 1)
    stop("binary must be a logical vector of length 1")
}

#' gwas
#'
#' Run a gwas using genetic data from a binary dosage file
#' @param data Data frame containing the subject ID, phenotype
#' and covariates
#' @param bdinfo Information about the binary dosage file returned
#' from the BinaryDosage::getbdinfo routine
#' @param outfile The file name for the results output.
#' @param minmaf Minimum minor allele frequency of SNPs to include
#' in analysis. SNPS that have less than 20 minor alleles observed
#' will be excluded from the analysis regardless of the value of
#' minmaf. A value of 0 indicates to use all the SNPs that have 20
#' minor alleles observed. Default value is 0.
#' @param blksize Size of blocks of SNPs to read in at one time.
#' Larger blocks can improve overall speed but require larger
#' amounts of computer memory. A value of 0 indicates to use the
#' recommended block size. Default value is 0.
#' @param binary Logical value indicating if the phenotype
#' is a binary value. Default value is false.
#' @return
#' 0
#' @export
#'
#' @examples
#' bdinfo <- readRDS(system.file("extdata/pdata_4_1.bdinfo", package = "GxEScanR"))
#' covdata <- readRDS(system.file("extdata/covdata.rds", package = "GxEScanR"))
#'
#' results <- gwas(data = covdata, bdinfo = bdinfo, binary = FALSE)
gwas <- function(data, bdinfo, snps, outfile, skipfile,
                 minmaf, blksize, binary) {
  # Set missing values to default values
  if (missing(snps) == TRUE)
    snps = "all"
  if (missing(outfile) == TRUE)
    outfile <- ""
  if (missing(skipfile) == TRUE)
    skipfile <- ""
  if (missing(minmaf) == TRUE)
    minmaf <- 0.
  if (missing(blksize))
    blksize <- 0L
  if (missing(binary) == TRUE)
    binary = TRUE
  validateinput(data, bdinfo, outfile, skipfile, minmaf, blksize, binary)
  snps <- subsetsnps(snps = snps,
                     snplist = bdinfo$snps$snpid)

  #####################################################
  ###       Subset data to subjects with complete
  ###       data and run baseline regression
  #####################################################
  data <- subsetdata(subdata = data,
                     bdinfo = bdinfo,
                     binary = binary,
                     mincov = 0)
  stddata <- matrix(data = 0,
                    nrow = nrow(data),
                    ncol = ncol(data))
  means <- numeric(ncol(data))
  stddevs <- numeric(ncol(data))
  stdmat(data, stddata, means, stddevs)
  print(paste(nrow(data), "subjects have complete data"))
  if (binary == TRUE)
    modfamily = "binomial"
  else
    modfamily = "gaussian"
  modformula = paste(colnames(data)[1], '~', '.')
  basemodel = glm(formula = modformula,
                  family = modfamily,
                  data = as.data.frame(data))

  #####################################################
  ###       Get information about the binary
  ###       dosage file and create the blocks
  ###       needed for reading the data
  #####################################################
  blkinfo <- assignblocks(nsub = nrow(bdinfo$samples),
                          nsnps = length(bdinfo$snps$snpid),
                          snploc = bdinfo$indices,
                          snpbytes = bdinfo$datasize,
                          reqblksize = blksize)
  if (bdinfo$additionalinfo$format == 1) {
    if (bdinfo$additionalinfo$subformat == 1)
      base = 0L
    else
      base = 1L
  } else {
    base <- 2L
  }
  subindex <- as.integer(rownames(data))
  
  if (binary == TRUE) {
    ncov <- length(basemodel$coefficients)
    beta0 <- numeric(ncov)
    beta0[2:ncov] <- basemodel$coefficients[2:ncov] * stddevs[2:ncov]
    beta0[1] <- basemodel$coefficients[1] +
                  sum(means[2:ncov] * basemodel$coefficients[2:ncov])
    return (logreggwas(bdinfo = bdinfo,
                       blkinfo = blkinfo,
                       snps = snps,
                       stddata = stddata,
                       subindex = subindex,
                       outfile = outfile,
                       skipfile = skipfile,
                       beta0 = beta0,
                       minmaf = minmaf,
                       base = base))
  }
  
  return (linreggwas(bdinfo = bdinfo,
                     blkinfo = blkinfo,
                     snps = snps,
                     stddata = stddata,
                     subindex = subindex,
                     outfile = outfile,
                     skipfile = skipfile,
                     minmaf = minmaf,
                     base = base))
}

#####################################################
###                      GWIS
#####################################################

#' gwis
#'
#' Run a gwas using genetic data from a binary dosage file
#' @param data Data frame containing the subject ID, phenotype
#' and covariates
#' @param bdinfo Information about the binary dosage file returned
#' from the BinaryDosage::getbdinfo routine
#' @param outfile The file name for the results output.
#' @param minmaf Minimum minor allele frequency of SNPs to include
#' in analysis. SNPS that have less than 20 minor alleles observed
#' will be excluded from the analysis regardless of the value of
#' minmaf. A value of 0 indicates to use all the SNPs that have 20
#' minor alleles observed. Default value is 0.
#' @param blksize Size of blocks of SNPs to read in at one time.
#' Larger blocks can improve overall speed but require larger
#' amounts of computer memory. A value of 0 indicates to use the
#' recommended block size. Default value is 0.
#' @param binary Logical value indicating if the phenotype
#' is a binary value. Default value is false.
#' @return
#' 0
#' @export
#'
#' @examples
#' bdinfo <- readRDS(system.file("extdata/pdata_4_1.bdinfo", package = "GxEScanR"))
#' covdata <- readRDS(system.file("extdata/covdata.rds", package = "GxEScanR"))
#'
#' results <- gwis(data = covdata, bdinfo = bdinfo)
gwis <- function(data, bdinfo, snps, outfile, skipfile,
                 minmaf, blksize, binary) {
  # Set missing values to default values
  if (missing(snps) == TRUE)
    snps = "all"
  if (missing(outfile) == TRUE)
    outfile <- ""
  if (missing(skipfile) == TRUE)
    skipfile <- ""
  if (missing(minmaf) == TRUE)
    minmaf <- 0.
  if (missing(blksize))
    blksize <- 0L
  if (missing(binary) == TRUE)
    binary = TRUE
  validateinput(data, bdinfo, outfile, skipfile, minmaf, blksize, binary)
  snps <- subsetsnps(snps = snps,
                     snplist = bdinfo$snps$snpid)
  
  #####################################################
  ###       Subset data to subjects with complete
  ###       data and run baseline regression
  #####################################################
  data <- subsetdata(subdata = data,
                     bdinfo = bdinfo,
                     binary = binary,
                     mincov = 1)
  nsub <- nrow(data)
  ncov <- ncol(data)
  stddata <- matrix(data = 0,
                    nrow = nsub,
                    ncol = ncov)
  means <- numeric(ncov)
  stddevs <- numeric(ncov)
  stdmat(data, stddata, means, stddevs)
  print(paste(nsub, "subjects have complete data"))

  #####################################################
  ###       Get information about the binary
  ###       dosage file and create the blocks
  ###       needed for reading the data
  #####################################################
  blkinfo <- assignblocks(nsub = nrow(bdinfo$samples),
                          nsnps = length(bdinfo$snps$snpid),
                          snploc = bdinfo$indices,
                          snpbytes = bdinfo$datasize,
                          reqblksize = blksize)
  blkbuffer <- integer((max(blkinfo$blkbytes) + 3) %/% 4)
  numblks <- length(blkinfo$blkloc)
  if (bdinfo$additionalinfo$format == 1) {
    if (bdinfo$additionalinfo$subformat == 1)
      base = 0L
    else
      base = 1L
  } else {
    base <- 2L
  }
  
  subindex <- as.integer(rownames(data))
  if (binary == TRUE)
    return (logreggwis(bdinfo = bdinfo,
                       blkinfo = blkinfo,
                       snps = snps,
                       stddata = stddata,
                       subindex = subindex,
                       outfile = outfile,
                       skipfile = skipfile,
                       minmaf = minmaf,
                       base = base,
                       e = data[,ncol(data)]))
  return (linreggwis(bdinfo = bdinfo,
                     blkinfo = blkinfo,
                     snps = snps,
                     stddata = stddata,
                     subindex = subindex,
                     minmaf = minmaf,
                     outfile = outfile,
                     skipfile = skipfile,
                     base = base,
                     estddev = stddevs[length(stddevs)]))
}