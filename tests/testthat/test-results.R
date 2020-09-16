test_that("results", {
  bdinfofile <- system.file("extdata", "pdata_4_1.bdinfo", package = "GxEScanR")
  bdinfo <- readRDS(bdinfofile)
  bdinfo$filename <- system.file("extdata", "pdata_4_1.bdose", package = "GxEScanR")
  covdatafile <- system.file("extdata", "covdata.rds", package = "GxEScanR")
  covdata <- readRDS(covdatafile)
  
  outfile <- tempfile()
  skipfile <- tempfile()
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
  
  results <- gwas(data = covdata, bdinfo = bdinfo)
  results <- gwas(data = covdata, bdinfo = bdinfo, binary = FALSE)
  results <- gwis(data = covdata, bdinfo = bdinfo)
  results <- gwis(data = covdata, bdinfo = bdinfo, binary = FALSE)
  expect_equal(ncol(results), 6L)
})
