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
  explingwas <- readRDS(system.file("extdata", "lingwas.rds", package = "GxEScanR"))
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
  expect_true(all(abs(lingwas1$betag - explingwas$betag) < 1e-6))
  expect_true(all(abs(lingwas2$betag - explingwas$betag) < 1e-6))
  expect_true(all(abs(lingwas3$betag - explingwas$betag[1:2]) < 1e-6))
  expect_true(all(abs(lingwas1$lrtg - explingwas$lrtg) < 1e-6))
  expect_true(all(abs(lingwas2$lrtg - explingwas$lrtg) < 1e-6))
  expect_true(all(abs(lingwas3$lrtg - explingwas$lrtg[1:2]) < 1e-6))
  
  # Linear regression GWEIS
  explingweis <- readRDS(system.file("extdata", "lingweis.rds", package = "GxEScanR"))
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
  expect_true(all(abs(explingweis$betadg - lingweis1$betadg) < 1e-6))
  expect_true(all(abs(explingweis$betadg - lingweis2$betadg) < 1e-6))
  expect_true(all(abs(explingweis$betadg[4:5] - lingweis3$betadg) < 1e-6))
  expect_true(all(abs(explingweis$lrtdg - lingweis1$lrtdg) < 1e-6))
  expect_true(all(abs(explingweis$lrtdg - lingweis2$lrtdg) < 1e-6))
  expect_true(all(abs(explingweis$lrtdg[4:5] - lingweis3$lrtdg) < 1e-6))
  expect_true(all(abs(explingweis$betagxe - lingweis1$betagxe) < 1e-6))
  expect_true(all(abs(explingweis$betagxe - lingweis2$betagxe) < 1e-6))
  expect_true(all(abs(explingweis$betagxe[4:5] - lingweis3$betagxe) < 1e-6))
  expect_true(all(abs(explingweis$lrtgxe - lingweis1$lrtgxe) < 1e-6))
  expect_true(all(abs(explingweis$lrtgxe - lingweis2$lrtgxe) < 1e-6))
  expect_true(all(abs(explingweis$lrtgxe[4:5] - lingweis3$lrtgxe) < 1e-6))
  expect_true(all(abs(explingweis$lrt2df - lingweis1$lrt2df) < 1e-6))
  expect_true(all(abs(explingweis$lrt2df - lingweis2$lrt2df) < 1e-6))
  expect_true(all(abs(explingweis$lrt2df[4:5] - lingweis3$lrt2df) < 1e-6))
  
  # Logistic regression GWAS
  exploggwas <- readRDS(system.file("extdata", "loggwas.rds", package = "GxEScanR"))
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
  expect_true(all(abs(exploggwas$betag - loggwas1$betag) < 1e-6))
  expect_true(all(abs(exploggwas$betag - loggwas2$betag) < 1e-6))
  expect_true(all(abs(exploggwas$betag[1:2] - loggwas3$betag) < 1e-6))
  expect_true(all(abs(exploggwas$lrtg - loggwas1$lrtg) < 1e-6))
  expect_true(all(abs(exploggwas$lrtg - loggwas2$lrtg) < 1e-6))
  expect_true(all(abs(exploggwas$lrtg[1:2] - loggwas3$lrtg) < 1e-6))
  
  # Logistic regression GWEIS - binary covariate
  exploggweis <- readRDS(system.file("extdata", "loggweis.rds", package = "GxEScanR"))
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
  expect_true(all(abs(exploggweis$betadg - loggweis1$betadg) < 1e-6))
  expect_true(all(abs(exploggweis$betadg - loggweis2$betadg) < 1e-6))
  expect_true(all(abs(exploggweis$betadg[4:5] - loggweis3$betadg) < 1e-6))
  expect_true(all(abs(exploggweis$lrtdg - loggweis1$lrtdg) < 1e-6))
  expect_true(all(abs(exploggweis$lrtdg - loggweis2$lrtdg) < 1e-6))
  expect_true(all(abs(exploggweis$lrtdg[4:5] - loggweis3$lrtdg) < 1e-6))
  expect_true(all(abs(exploggweis$betagxe - loggweis1$betagxe) < 1e-6))
  expect_true(all(abs(exploggweis$betagxe - loggweis2$betagxe) < 1e-6))
  expect_true(all(abs(exploggweis$betagxe[4:5] - loggweis3$betagxe) < 1e-6))
  expect_true(all(abs(exploggweis$lrtgxe - loggweis1$lrtgxe) < 1e-6))
  expect_true(all(abs(exploggweis$lrtgxe - loggweis2$lrtgxe) < 1e-6))
  expect_true(all(abs(exploggweis$lrtgxe[4:5] - loggweis3$lrtgxe) < 1e-6))
  expect_true(all(abs(exploggweis$lrt2df - loggweis1$lrt2df) < 1e-6))
  expect_true(all(abs(exploggweis$lrt2df - loggweis2$lrt2df) < 1e-6))
  expect_true(all(abs(exploggweis$lrt2df[4:5] - loggweis3$lrt2df) < 1e-6))
  expect_true(all(abs(exploggweis$betaeg - loggweis1$betaeg) < 1e-6))
  expect_true(all(abs(exploggweis$betaeg - loggweis2$betaeg) < 1e-6))
  expect_true(all(abs(exploggweis$betaeg[4:5] - loggweis3$betaeg) < 1e-6))
  expect_true(all(abs(exploggweis$lrteg - loggweis1$lrteg) < 1e-6))
  expect_true(all(abs(exploggweis$lrteg - loggweis2$lrteg) < 1e-6))
  expect_true(all(abs(exploggweis$lrteg[4:5] - loggweis3$lrteg) < 1e-6))
  expect_true(all(abs(exploggweis$lrt3df - loggweis1$lrt3df) < 1e-6))
  expect_true(all(abs(exploggweis$lrt3df - loggweis2$lrt3df) < 1e-6))
  expect_true(all(abs(exploggweis$lrt3df[4:5] - loggweis3$lrt3df) < 1e-6))
  expect_true(all(abs(exploggweis$betacase - loggweis1$betacase) < 1e-6))
  expect_true(all(abs(exploggweis$betacase - loggweis2$betacase) < 1e-6))
  expect_true(all(abs(exploggweis$betacase[4:5] - loggweis3$betacase) < 1e-6))
  expect_true(all(abs(exploggweis$lrtcase - loggweis1$lrtcase) < 1e-6))
  expect_true(all(abs(exploggweis$lrtcase - loggweis2$lrtcase) < 1e-6))
  expect_true(all(abs(exploggweis$lrtcase[4:5] - loggweis3$lrtcase) < 1e-6))
  expect_true(all(abs(exploggweis$betactrl - loggweis1$betactrl) < 1e-6))
  expect_true(all(abs(exploggweis$betactrl - loggweis2$betactrl) < 1e-6))
  expect_true(all(abs(exploggweis$betactrl[4:5] - loggweis3$betactrl) < 1e-6))
  expect_true(all(abs(exploggweis$lrtctrl - loggweis1$lrtctrl) < 1e-6))
  expect_true(all(abs(exploggweis$lrtctrl - loggweis2$lrtctrl) < 1e-6))
  expect_true(all(abs(exploggweis$lrtctrl[4:5] - loggweis3$lrtctrl) < 1e-6))
  
  # Logistic regression GWEIS - continuous covariate
  exploggweisc <- readRDS(system.file("extdata", "loggweisc.rds", package = "GxEScanR"))
  loggweis1c <- gweis(data = covdata2,
                      bdinfo = bdinfo)
  loggweis2c <- gweis(data = covdata2,
                      bdinfo = bdinfo,
                      outfile = outfile,
                      minmaf = 0.001,
                      blksize = 2)
  expect_equal(loggweis2c, 0)
  loggweis2c <- read.table(outfile, header = TRUE, sep = '\t')
  loggweis3c <- gweis(data = covdata2,
                      bdinfo = bdinfo,
                      snps = 4:5,
                      skipfile = skipfile,
                      blksize = 2)
  loggweis3cskip <- read.table(skipfile, header = TRUE, sep = '\t')
  expect_true(all(loggweis3cskip$reason == 1))
  expect_true(all(abs(exploggweisc$betadg - loggweis1c$betadg) < 1e-6))
  expect_true(all(abs(exploggweisc$betadg - loggweis2c$betadg) < 1e-6))
  expect_true(all(abs(exploggweisc$betadg[4:5] - loggweis3c$betadg) < 1e-6))
  expect_true(all(abs(exploggweisc$lrtdg - loggweis1c$lrtdg) < 1e-6))
  expect_true(all(abs(exploggweisc$lrtdg - loggweis2c$lrtdg) < 1e-6))
  expect_true(all(abs(exploggweisc$lrtdg[4:5] - loggweis3c$lrtdg) < 1e-6))
  expect_true(all(abs(exploggweisc$betagxe - loggweis1c$betagxe) < 1e-6))
  expect_true(all(abs(exploggweisc$betagxe - loggweis2c$betagxe) < 1e-6))
  expect_true(all(abs(exploggweisc$betagxe[4:5] - loggweis3c$betagxe) < 1e-6))
  expect_true(all(abs(exploggweisc$lrtgxe - loggweis1c$lrtgxe) < 1e-6))
  expect_true(all(abs(exploggweisc$lrtgxe - loggweis2c$lrtgxe) < 1e-6))
  expect_true(all(abs(exploggweisc$lrtgxe[4:5] - loggweis3c$lrtgxe) < 1e-6))
  expect_true(all(abs(exploggweisc$lrt2df - loggweis1c$lrt2df) < 1e-6))
  expect_true(all(abs(exploggweisc$lrt2df - loggweis2c$lrt2df) < 1e-6))
  expect_true(all(abs(exploggweisc$lrt2df[4:5] - loggweis3c$lrt2df) < 1e-6))
  expect_true(all(abs(exploggweisc$betaeg - loggweis1c$betaeg) < 1e-6))
  expect_true(all(abs(exploggweisc$betaeg - loggweis2c$betaeg) < 1e-6))
  expect_true(all(abs(exploggweisc$betaeg[4:5] - loggweis3c$betaeg) < 1e-6))
  expect_true(all(abs(exploggweisc$lrteg - loggweis1c$lrteg) < 1e-6))
  expect_true(all(abs(exploggweisc$lrteg - loggweis2c$lrteg) < 1e-6))
  expect_true(all(abs(exploggweisc$lrteg[4:5] - loggweis3c$lrteg) < 1e-6))
  expect_true(all(abs(exploggweisc$lrt3df - loggweis1c$lrt3df) < 1e-6))
  expect_true(all(abs(exploggweisc$lrt3df - loggweis2c$lrt3df) < 1e-6))
  expect_true(all(abs(exploggweisc$lrt3df[4:5] - loggweis3c$lrt3df) < 1e-6))
  expect_true(all(abs(exploggweisc$betacase - loggweis1c$betacase) < 1e-6))
  expect_true(all(abs(exploggweisc$betacase - loggweis2c$betacase) < 1e-6))
  expect_true(all(abs(exploggweisc$betacase[4:5] - loggweis3c$betacase) < 1e-6))
  expect_true(all(abs(exploggweisc$lrtcase - loggweis1c$lrtcase) < 1e-6))
  expect_true(all(abs(exploggweisc$lrtcase - loggweis2c$lrtcase) < 1e-6))
  expect_true(all(abs(exploggweisc$lrtcase[4:5] - loggweis3c$lrtcase) < 1e-6))
  expect_true(all(abs(exploggweisc$betactrl - loggweis1c$betactrl) < 1e-6))
  expect_true(all(abs(exploggweisc$betactrl - loggweis2c$betactrl) < 1e-6))
  expect_true(all(abs(exploggweisc$betactrl[4:5] - loggweis3c$betactrl) < 1e-6))
  expect_true(all(abs(exploggweisc$lrtctrl - loggweis1c$lrtctrl) < 1e-6))
  expect_true(all(abs(exploggweisc$lrtctrl - loggweis2c$lrtctrl) < 1e-6))
  expect_true(all(abs(exploggweisc$lrtctrl[4:5] - loggweis3c$lrtctrl) < 1e-6))
})
