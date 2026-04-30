library(GxEScanR)

bdinfo <- BinaryDosage::getbdinfo(system.file("extdata", "gendata.bdose", package = "GxEScanR"))
subdata <- readRDS(system.file("extdata", "subdata.rds", package = "GxEScanR"))
subdata <- subdata[complete.cases(subdata), ]
subdata <- subdata[subdata$subid %in% bdinfo$samples$sid, ]

linearmodel <- glm(y_linear ~ x2 + x1, data = subdata)
logisticmodel <- glm(y_logistic ~ x2 + x1, family = binomial, data = subdata)

test_that("gweis.mem returns 1 for invalid tests argument", {
  result <- gweis.mem(gemdl = linearmodel, subids = subdata$subid, tests = "invalid")
  expect_equal(result, 1)
})

test_that("gweis.mem errors when bg_go requested without gomdl", {
  expect_error(
    gweis.mem(gemdl = linearmodel, subids = subdata$subid, tests = "bg_go"),
    "gomdl must be provided"
  )
})

test_that("gweis.mem returns a list with expected structure", {
  mem <- gweis.mem(gemdl = linearmodel, subids = subdata$subid,
                   tests = c("bg_ge", "bg_gxe", "bgxe", "joint"))
  expect_type(mem, "list")
  expect_named(mem, c("tests", "subids", "gomdlmem", "gemdlmem", "gxemdlmem",
                      "gxe0mdlmem", "egmdlmem", "caseids", "casemdlmem",
                      "ctrlids", "ctrlmdlmem"))
  expect_equal(length(mem$subids), nrow(subdata))
  expect_length(mem$tests, 8)
})

test_that("gweis.mem allocates correct model memory for requested tests", {
  mem <- gweis.mem(gemdl = linearmodel, subids = subdata$subid, tests = "bg_ge")
  expect_false(is.null(mem$gemdlmem))
  expect_null(mem$gomdlmem)
  expect_null(mem$gxemdlmem)
})

test_that("rungweis writes output file with correct columns (linear)", {
  mem <- gweis.mem(gemdl = linearmodel, subids = subdata$subid,
                   tests = c("bg_ge", "bgxe", "joint"))
  outfile <- tempfile(fileext = ".txt")
  on.exit(unlink(outfile))

  rungweis(gweismem = mem, bdinfo = bdinfo,
           snps = seq_len(nrow(bdinfo$snps)), outfilename = outfile)

  expect_true(file.exists(outfile))
  results <- read.table(outfile, header = TRUE, sep = "\t")
  expected_cols <- c("snpid", "chr", "loc", "ref", "alt", "aaf",
                     "bg_ge", "bg_ge_lrt", "bg_gxe", "bgxe", "bgxe_lrt", "joint_lrt")
  expect_true(all(expected_cols %in% names(results)))
  expect_gt(nrow(results), 0)
})

test_that("rungweis writes output file with correct columns (logistic)", {
  mem <- gweis.mem(gemdl = logisticmodel, subids = subdata$subid,
                   tests = c("bg_ge", "bg_eg", "bg_case", "bg_ctrl"))
  outfile <- tempfile(fileext = ".txt")
  on.exit(unlink(outfile))

  rungweis(gweismem = mem, bdinfo = bdinfo,
           snps = seq_len(nrow(bdinfo$snps)), outfilename = outfile)

  expect_true(file.exists(outfile))
  results <- read.table(outfile, header = TRUE, sep = "\t")
  expect_true(all(c("aaf", "aaf_case", "aaf_ctrl",
                    "bg_ge", "bg_ge_lrt",
                    "bg_eg", "bg_eg_lrt",
                    "bg_case", "bg_case_lrt",
                    "bg_ctrl", "bg_ctrl_lrt") %in% names(results)))
})

test_that("rungweis respects maf filter", {
  mem <- gweis.mem(gemdl = linearmodel, subids = subdata$subid, tests = "bg_ge")
  outfile <- tempfile(fileext = ".txt")
  on.exit(unlink(outfile))

  rungweis(gweismem = mem, bdinfo = bdinfo,
           snps = seq_len(nrow(bdinfo$snps)), outfilename = outfile, maf = 0.49)

  results <- read.table(outfile, header = TRUE, sep = "\t")
  expect_true(all(results$aaf >= 0.49 & results$aaf <= 0.51))
})
