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

  # Linear regression GWEIS
  lingweis1 <- gweis(data = covdata,
                     bdinfo = bdinfo,
                     binary = FALSE)
  lingweis2 <- gweis(data = covdata,
                     bdinfo = bdinfo,
                     outfile = outfile,
                     minmaf = 0.05,
                     blksize = 2,
                     binary = FALSE)
  expect_equal(lingweis2, 0)
  lingweis2 <- read.table(outfile, header = TRUE, sep = '\t')
  lingweis3 <- gweis(data = covdata,
                     bdinfo = bdinfo,
                     snps = 4:5,
                     skipfile = skipfile,
                     blksize = 2,
                     binary = FALSE)
  lingweis3skip <- read.table(skipfile, header = TRUE, sep = '\t')
  expect_true(all(lingweis3skip$reason == 1))
  expect_true(all(abs(lingweis1$betadg - lingweis2$betadg) < 1e-6))
  expect_true(all(abs(lingweis1$betadg[4:5] - lingweis3$betadg) < 1e-6))
  
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

  # Logistic regression GWEIS - binary covariate
  loggweis1 <- gweis(data = covdata,
                     bdinfo = bdinfo)
  loggweis2 <- gweis(data = covdata,
                     bdinfo = bdinfo,
                     outfile = outfile,
                     minmaf = 0.05,
                     blksize = 2)
  expect_equal(loggweis2, 0)
  loggweis2 <- read.table(outfile, header = TRUE, sep = '\t')
  loggweis3 <- gweis(data = covdata,
                     bdinfo = bdinfo,
                     snps = 4:5,
                     skipfile = skipfile,
                     blksize = 2)
  loggweis3skip <- read.table(skipfile, header = TRUE, sep = '\t')
  expect_true(all(loggweis3skip$reason == 1))
  expect_true(all(abs(loggweis1$betadg - loggweis2$betadg) < 1e-6))
  expect_true(all(abs(loggweis1$betadg[4:5] - loggweis3$betadg) < 1e-6))
  
  # Logistic regression GWEIS - continuous covariate
  loggweis1c <- gweis(data = covdata2,
                      bdinfo = bdinfo)
  loggweis2c <- gweis(data = covdata2,
                      bdinfo = bdinfo,
                      outfile = outfile,
                      minmaf = 0.001,
                      blksize = 2)
  expect_equal(loggweis2c, 0)
  loggweis2c <- read.table(outfile, header = TRUE, sep = '\t')
  loggweis3c <- gweis(data = covdata,
                      bdinfo = bdinfo,
                      snps = 4:5,
                      skipfile = skipfile,
                      blksize = 2)
  loggweis3cskip <- read.table(skipfile, header = TRUE, sep = '\t')
  expect_true(all(loggweis3cskip$reason == 1))
  expect_true(all(abs(loggweis1c$betadg - loggweis2c$betadg) < 1e-6))
  expect_true(all(abs(loggweis1c$betadg[4:5] - loggweis3c$betadg) < 1e-6))
})
