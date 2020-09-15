test_that("gwas input", {
  ##############################################################
  #                 Data sets with good values
  ##############################################################
  
  # Data for phenotype/covariate data
  substouse <- c(1,2,4,8,16,32,64,128,256)
  sid <- paste(rep("I", length(substouse)), substouse, sep = "")
  fid <- paste(rep("F", length(substouse)), substouse, sep = "")
  y = c(0, 0, 0, 0, 1, 1, 1, 1, 1)
  x = c(1, 2, 4, 3, 2, 5, 3, 4, 5)
  
  # Subject and genetic data sets using only subject IDs
  data <- data.frame(sid = sid,
                     y = y,
                     x = x,
                     stringsAsFactors = FALSE)
  bdinfofile <- system.file("extdata", "pdata_4_1.bdinfo", package = "GxEScanR")
  bdinfo <- readRDS(bdinfofile)
  bdinfo$filename <- system.file("extdata", "pdata_4_1.bdose", package = "GxEScanR")
  # Large bdinfo data for testing blksize
  bdinfofile <- system.file("extdata", "largebdinfo.rds", package = "GxEScanR")
  bdinfobig <- readRDS(bdinfofile)
  
  # Subject and genetic data sets using subject IDs and family IDs
  dataf <- data.frame(fid = fid,
                      sid = sid,
                      y = y,
                      x = x,
                      stringsAsFactors = FALSE)
  bdinfof <- bdinfo
  bdinfof$usesfid <- TRUE
  bdinfof$samples$fid <- paste(rep("F", nrow(bdinfo$samples)),
                               1:nrow(bdinfo$samples),
                               sep = "")

  # Other data values for testing that are valid  
  outfile <- ""
  minmaf <- 0.05
  blksize <- 0
  binary <- FALSE
 
  ##############################################################
  #                 Testing the genetic data
  ##############################################################
  
  # Testing if bdinfo has information about genetic data
  expect_error(gwas(data = data,
                    bdinfo = 1,
                    outfile = outfile),
               "bdinfo not a genetic-info class")

  # Testing bdinfo is information about a binary dosage file
  addinfoclass <- class(bdinfo$additionalinfo)
  class(bdinfo$additionalinfo) <- "xyz"
  expect_error(gwas(data = data,
                    bdinfo = bdinfo,
                    outfile = outfile),
               "bdinfo does not have information about a binary dosage file")
  class(bdinfo$additionalinfo) <- addinfoclass

  ##############################################################
  #                 Testing outfile value
  ##############################################################
  
  # Testing the outfile value is a string value
  expect_error(gwas(data = data,
                    bdinfo = bdinfo,
                    outfile = 1),
               "outfile must be a character value")
  
  # Testing the outfile value is a single string value
  expect_error(gwas(data = data,
                    bdinfo = bdinfo,
                    outfile = c("a", "b")),
               "outfile must be a character vector of length 1")
  
  # Testing the skipfile value is a string value
  expect_error(gwas(data = data,
                    bdinfo = bdinfo,
                    skipfile = 1),
               "skipfile must be a character value")
  
  # Testing the skipfile value is a single string value
  expect_error(gwas(data = data,
                    bdinfo = bdinfo,
                    skipfile = c("a", "b")),
               "skipfile must be a character vector of length 1")

  ##############################################################
  #                 Testing minmaf value
  ##############################################################
  
  # Testing minmaf is a numeric value
  expect_error(gwas(data = data,
                    bdinfo = bdinfo,
                    outfile = outfile,
                    minmaf = "a"),
               "minmaf must be a numeric value")
  
  # Testing minmaf is a single numeric value
  expect_error(gwas(data = data,
                    bdinfo = bdinfo,
                    outfile = outfile,
                    minmaf = 1:2),
               "minmaf must be a numeric vector of length 1")
  
  # Testing minmaf is a number value from 0 to 0.25 inclusive
  expect_error(gwas(data = data,
                    bdinfo = bdinfo,
                    outfile = outfile,
                    minmaf = 1.2),
               "minmaf must be a value from 0 to 0.25, inclusive")
  
  ##############################################################
  #                 Testing blocksize value
  ##############################################################
  
  # Testing if blocksize is an numeric value
  expect_error(gwas(data = data,
                    bdinfo = bdinfo,
                    outfile = outfile,
                    blksize = "a"),
               "blksize must be an integer")
  
  # Testing if blocksize is a single integer value
  expect_error(gwas(data = data,
                    bdinfo = bdinfo,
                    outfile = outfile,
                    blksize = 1:2),
               "blksize must be an integer vector of length 1")
  
  # Testing if blksize is an integer value
  expect_error(gwas(data = data,
                    bdinfo = bdinfo,
                    outfile = outfile,
                    blksize = 1.2),
               "blksize must be an integer")
  
  # Testing if blksize is an positive integer value
  expect_error(gwas(data = data,
                    bdinfo = bdinfo,
                    outfile = outfile,
                    blksize = -2),
               "blksize must be greater than or equal to 0")
  
  # Testing if the blksize is too large
  expect_error(gwas(data = data,
                    bdinfo = bdinfobig,
                    outfile = outfile,
                    minmaf = minmaf,
                    blksize = 10001),
               "Requested block size greater than twice the recommended block size")
  
  ##############################################################
  #                 Testing binary value
  ##############################################################
  
  # Checking if binary is a logical value
  expect_error(gwas(data = data,
                    bdinfo = bdinfo,
                    outfile = outfile,
                    binary = 1),
               "binary must be a logical value")

  # Checking if binary is a single logical value
  expect_error(gwas(data = data,
                    bdinfo = bdinfo,
                    outfile = outfile,
                    binary = c(FALSE, TRUE)),
               "binary must be a logical vector of length 1")

  ##############################################################
  #                 Testing the subject data
  ##############################################################
  
  # Testing data is a data frame 
  expect_error(gwas(data = 1,
                    bdinfo = bdinfo,
                    outfile = outfile),
               "data must be a data frame")
  
  # Check if subject data has at least two columns
  dataerror <- data.frame(sid = data$sid, stringsAsFactors = FALSE)  
  expect_error(gwas(data = dataerror,
                    bdinfo = bdinfo,
                    outfile = outfile,
                    binary = binary),
               "There must me at least two columns in the subject data")

  # Check if subject data first column is a character value
  dataerror <- data
  dataerror$sid <- rep(1, nrow(data))
  expect_error(gwas(data = dataerror,
                    bdinfo = bdinfo,
                    outfile = outfile,
                    binary = binary),
               "First column of subject data must be a character value")
  
  # Check if any subjects have complete data
  dataerror <- data
  dataerror$y[1] <- NA
  dataerror$x[2:nrow(data)] <- NA
  expect_error(gwas(data = dataerror,
                    bdinfo = bdinfo,
                    outfile = outfile,
                    binary = binary),
               "No subjects have complete phenotype/covariate data")
  
  # Check if subject data has at least three columns if family id is used
  dataerror <- dataf[,1:2]
  expect_error(gwas(data = dataerror,
                    bdinfo = bdinfof,
                    outfile = outfile,
                    binary = binary),
               "When using family ID, subject data must have at least 3 columns")
  
  # Check if family id is a character value
  dataerror <- dataf[,c(1,3:4)]
  expect_error(gwas(data = dataerror,
                    bdinfo = bdinfof,
                    outfile = outfile,
                    binary = binary),
               "When using family ID, the first two columns must be character values")
  
  # Check if phenotype and covariate values are numeric
  dataerror <- data
  dataerror$x[1:(nrow(dataerror) - 1)] <- NA
  expect_error(gwas(data = dataerror,
                    bdinfo = bdinfo,
                    outfile = outfile,
                    binary = binary),
               "No subjects have complete data")
  
  # Check if there are two phenotype values
  dataerror <- data
  dataerror$y[1] <- 9
  expect_error(gwas(data = dataerror,
                    bdinfo = bdinfo,
                    outfile = outfile,
                    binary = TRUE),
               "When using a binary outcome there must only be two values")
  
  # Check if the two phenotype values are 0,1
  dataerror <- data
  dataerror$y <- dataerror$y + 1
  expect_error(gwas(data = dataerror,
                    bdinfo = bdinfo,
                    outfile = outfile,
                    binary = TRUE),
               "When using a binary outcome must be coded 0,1")

  # Check if snps is a character or positive integer array
  expect_error(subsetsnps(snps = TRUE,
                          snplist = bdinfo$snps$snpid),
               "snps must be a character or integer array")
  expect_error(subsetsnps(snps = 1.2,
                          snplist = bdinfo$snps$snpid),
               "snps must be a character or integer array")
  expect_error(subsetsnps(snps = 0,
                          snplist = bdinfo$snps$snpid),
               "snp indices must be positive")
  expect_error(subsetsnps(snps = 10,
                          snplist = bdinfo$snps$snpid),
               "at least one snp index is greater than the number of SNPs available")
  expect_error(subsetsnps(snps = "1:10010",
                          snplist = bdinfo$snps$snpid),
               "No matching SNPs found")
  expect_error(subsetsnps(snps = character(),
                          snplist = bdinfo$snps$snpid),
               "No SNPs selected")
  dataerror <- data[,1:2]
  expect_error(gwis(data = dataerror,
                    bdinfo = bdinfo),
               "Subject data has no covariates")
  
  
  x <- subsetsnps(snps = c("1:10001", "1:10004", "1:10010"),
                           snplist = bdinfo$snps$snpid)
  expect_equal(x[1], TRUE)
  expect_equal(x[2], FALSE)
  expect_equal(x[3], FALSE)
  expect_equal(x[4], TRUE)
  expect_equal(x[5], FALSE)
  
  # Checking subsetting of subject data
  expect_error(gwas(data = data,
                    bdinfo = bdinfof),
               "When using family ID, the first two columns must be character values")
  
  
  bdinfofile <- system.file("extdata", "pdata_4_1.bdinfo", package = "GxEScanR")
  bdinfo <- readRDS(bdinfofile)
  bdinfo$filename <- system.file("extdata", "pdata_4_1.bdose", package = "GxEScanR")
  covdatafile <- system.file("extdata", "covdata.rds", package = "GxEScanR")
  covdata <- readRDS(covdatafile)
  results <- gwas(data = covdata, bdinfo = bdinfo)
  results <- gwas(data = covdata, bdinfo = bdinfo, binary = FALSE)
  results <- gwis(data = covdata, bdinfo = bdinfo)
  results <- gwis(data = covdata, bdinfo = bdinfo, binary = FALSE)
  
})
