context("test-snpindices")

test_that("Checking SNP indices", {
  snpsToFindI <- c(seq(1L, 3L), 5L)
  snpsToFindC <- paste("SNP", snpsToFindI, sep = "")
  numSNPs <- 6L
  snpsInGeneticFile <- paste("SNP", seq(1L,numSNPs), sep = "")
  
  expect_error(GetSNPIndices(snpsToFind = 3.,
                             snpsInGeneticFile = snpsInGeneticFile,
                             numSNPs = numSNPs
                             ), "snps is not a character or integer vector")
  expect_error(GetSNPIndices(snpsToFind = c(snpsToFindC, "X"),
                             snpsInGeneticFile = snpsInGeneticFile,
                             numSNPs = numSNPs
                             ), "Not all SNPs found")
  expect_error(GetSNPIndices(snpsToFind = c(snpsToFindC, "SNP1"),
                             snpsInGeneticFile = snpsInGeneticFile,
                             numSNPs = numSNPs
                             ), "SNP values are not unique")
  expect_error(GetSNPIndices(snpsToFind = c(snpsToFindI, 5L),
                             snpsInGeneticFile = snpsInGeneticFile,
                             numSNPs = numSNPs
                             ), "SNP values are not unique")
  expect_error(GetSNPIndices(snpsToFind = c(snpsToFindI, 0L),
                             snpsInGeneticFile = snpsInGeneticFile,
                             numSNPs = numSNPs
                             ), "snp indices must be positive")
  expect_error(GetSNPIndices(snpsToFind = c(snpsToFindI, 8L),
                             snpsInGeneticFile = snpsInGeneticFile,
                             numSNPs = numSNPs
                             ), "Value of largest snp index greater than number of SNPs")
  expect_equal(GetSNPIndices(snpsToFind = snpsToFindC,
                             snpsInGeneticFile = snpsInGeneticFile,
                             numSNPs = numSNPs
                             ), snpsToFindI)
})
