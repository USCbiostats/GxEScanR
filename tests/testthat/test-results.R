test_that("results", {
  bdinfofile <- system.file("extdata", "pdata_4_1.bdinfo", package = "GxEScanR")
  bdinfo <- readRDS(bdinfofile)
  bdinfo$filename <- system.file("extdata", "pdata_4_1.bdose", package = "GxEScanR")
  covdatafile <- system.file("extdata", "covdata.rds", package = "GxEScanR")
  covdata <- readRDS(covdatafile)
  covdata2 <- covdata
  covdata2$e <- covdata2$e + 1
  
  outfile <- tempfile()
  skipfile <- tempfile()
  
  # Linear regression GWAS
  lingwas1 <- gwas(data = covdata,
                   bdinfo = bdinfo,
                   binary = FALSE)
  lingwas2 <- gwas(data = covdata,
                   bdinfo = bdinfo,
                   outfile = outfile,
                   minmaf = 0.05,
                   blksize = 2,
                   binary = FALSE)
  expect_equal(lingwas2, 0)
  lingwas2 <- read.table(outfile, header = TRUE, sep = '\t')
  lingwas3 <- gwas(data = covdata,
                   bdinfo = bdinfo,
                   snps = 1:2,
                   skipfile = skipfile,
                   blksize = 2,
                   binary = FALSE)
  lingwas3skip <- read.table(skipfile, header = TRUE, sep = '\t')
  expect_true(all(lingwas3skip$reason == 1))
  expect_true(all(abs(lingwas1$betag - lingwas2$betag) < 1e-6))
  expect_true(all(abs(lingwas1$betag[1:2] - lingwas3$betag) < 1e-6))

  # Linear regression GWIS
  lingwis1 <- gwis(data = covdata,
                   bdinfo = bdinfo,
                   binary = FALSE)
  lingwis2 <- gwis(data = covdata,
                   bdinfo = bdinfo,
                   outfile = outfile,
                   minmaf = 0.05,
                   blksize = 2,
                   binary = FALSE)
  expect_equal(lingwis2, 0)
  lingwis2 <- read.table(outfile, header = TRUE, sep = '\t')
  lingwis3 <- gwis(data = covdata,
                   bdinfo = bdinfo,
                   snps = 4:5,
                   skipfile = skipfile,
                   blksize = 2,
                   binary = FALSE)
  lingwis3skip <- read.table(skipfile, header = TRUE, sep = '\t')
  expect_true(all(lingwis3skip$reason == 1))
  expect_true(all(abs(lingwis1$betadg - lingwis2$betadg) < 1e-6))
  expect_true(all(abs(lingwis1$betadg[4:5] - lingwis3$betadg) < 1e-6))
  
  # Logistic regression GWAS
  loggwas1 <- gwas(data = covdata,
                   bdinfo = bdinfo)
  loggwas2 <- gwas(data = covdata,
                   bdinfo = bdinfo,
                   outfile = outfile,
                   minmaf = 0.001,
                   blksize = 2)
  expect_equal(loggwas2, 0)
  loggwas2 <- read.table(outfile, header = TRUE, sep = '\t')
  loggwas3 <- gwas(data = covdata,
                   bdinfo = bdinfo,
                   snps = 1:2,
                   skipfile = skipfile,
                   blksize = 2)
  loggwas3skip <- read.table(skipfile, header = TRUE, sep = '\t')
  expect_true(all(loggwas3skip$reason == 1))
  expect_true(all(abs(loggwas1$betag - loggwas2$betag) < 1e-6))
  expect_true(all(abs(loggwas1$betag[1:2] - loggwas3$betag) < 1e-6))

  # Logistic regression GWIS - binary covariate
  loggwis1 <- gwis(data = covdata,
                   bdinfo = bdinfo)
  loggwis2 <- gwis(data = covdata,
                   bdinfo = bdinfo,
                   outfile = outfile,
                   minmaf = 0.05,
                   blksize = 2)
  expect_equal(loggwis2, 0)
  loggwis2 <- read.table(outfile, header = TRUE, sep = '\t')
  loggwis3 <- gwis(data = covdata,
                   bdinfo = bdinfo,
                   snps = 4:5,
                   skipfile = skipfile,
                   blksize = 2)
  loggwis3skip <- read.table(skipfile, header = TRUE, sep = '\t')
  expect_true(all(loggwis3skip$reason == 1))
  expect_true(all(abs(loggwis1$betadg - loggwis2$betadg) < 1e-6))
  expect_true(all(abs(loggwis1$betadg[4:5] - loggwis3$betadg) < 1e-6))
  
  # Logistic regression GWIS - continuous covariate
  loggwis1c <- gwis(data = covdata2,
                    bdinfo = bdinfo)
  loggwis2c <- gwis(data = covdata2,
                    bdinfo = bdinfo,
                    outfile = outfile,
                    minmaf = 0.001,
                    blksize = 2)
  expect_equal(loggwis2c, 0)
  loggwis2c <- read.table(outfile, header = TRUE, sep = '\t')
  loggwis3c <- gwis(data = covdata,
                    bdinfo = bdinfo,
                    snps = 4:5,
                    skipfile = skipfile,
                    blksize = 2)
  loggwis3cskip <- read.table(skipfile, header = TRUE, sep = '\t')
  expect_true(all(loggwis3cskip$reason == 1))
  expect_true(all(abs(loggwis1c$betadg - loggwis2c$betadg) < 1e-6))
  expect_true(all(abs(loggwis1c$betadg[4:5] - loggwis3c$betadg) < 1e-6))
})
