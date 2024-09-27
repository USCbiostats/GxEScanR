#' @importFrom lsReg lsReg
#' @importFrom lsReg runtest
#' @importFrom BinaryDosage getsnp
NULL

#' Routine to allocate memory needed to perform a GWEIS.
#'
#' @param gomdl The results from glm for the gene-only model. This model
#' contains the outcome and all the covariates except the covariate that
#' the gene interaction is being tested for. This can be NULL.
#' @param gemdl The results from glm for the gene-environment model. This model
#' contains the outcome and all covariates of interest with the last covariate
#' listed in the model being the covaraiate that the gene interaction is being
#' tested for.
#' @param subids A character vector of subject IDs that line up with the
#' data that was used in the models that are passed to this routine
#' @param tests The list of tests to perform. These can be any combination of
#' the following values "bg_go", "bg_ge", "bg_gxe", "bgxe", "joint", "bg_eg",
#' "bg_case", "bg_ctrl"
#'
#' @return List containing allocated memory to perform the specified GWEIS.
#' This value is passed to the rungweis routine.
#' @export
gweis.mem <- function(gomdl, gemdl, subids, tests) {
  x <- match(tests, c("bg_go", "bg_ge", "bg_gxe", "bgxe", "joint", "bg_eg", "bg_case", "bg_ctrl"))
  if (all(is.na(x) == FALSE) == FALSE)
    return (1)
  mdls <- logical(8)
  mdls[x] <- TRUE

  if (mdls[1] == TRUE)
    gomdlmem <- lsReg::lsReg(gomdl, 1, "lrt")
  else
    gomdlmem <- NULL
  if (mdls[2] == TRUE || mdls[4] == TRUE)
    gemdlmem <- lsReg::lsReg(gemdl, 1, "lrt")
  else
    gemdlmem <- NULL
  if (mdls[3] == TRUE || mdls[4] == TRUE|| mdls[5] == TRUE)
    gxemdlmem <- lsReg::lsReg(gemdl, 2, "lrt")
  else
    gxemdlmem <- NULL
  if (mdls[3] == TRUE)
    gxe0mdlmem <- lsReg::lsReg(gemdl, 1, "lrt")
  else
    gxe0mdlmem <- NULL

  f <- deparse(gemdl$formula)
  f <- gsub("\"", "", f)
  f <- unlist(strsplit(f, "~"))
  f <- unlist(strsplit(f, "\\+"))
  f2 <- paste(f[2:(length(f) -1)], collapse = "+")
  formula <- paste(f[length(f)],
                   f2,
                   sep = "~")

  evalues <- sort(unique(gemdl$model[,ncol(gemdl$model)]))
  if (length(evalues) == 2) {
    if (all(evalues == c(0, 1)) == TRUE)
      family = binomial
    else
      family = gaussian
  } else {
    family = gaussian
  }

  if (mdls[6] == TRUE) {
    emdl <- glm(formula = formula,
                data = gemdl$data,
                family = family)
    egmdlmem <- lsReg::lsReg(emdl, 1, "lrt")
  } else {
    egmdlmem <- NULL
  }

  if (gemdl$family$family == "binomial") {
    if (mdls[7] == TRUE) {
      emdl <- glm(formula = formula,
                  data = gemdl$data[gemdl$model[,1] == 1,],
                  family = family)
      casemdlmem <- lsReg::lsReg(emdl, 1, "lrt")
      caseids <- subids[gemdl$model[,1] == 1]
    } else {
      casemdlmem <- NULL
      caseids <- NULL
    }
    if (mdls[8] == TRUE) {
      emdl <- glm(formula = formula,
                  data = gemdl$data[gemdl$model[,1] == 0,],
                  family = family)
      ctrlmdlmem <- lsReg::lsReg(emdl, 1, "lrt")
      ctrlids <- subids[gemdl$model[,1] == 0]
    } else {
      ctrlmdlmem <- NULL
      ctrlids <- NULL
    }
  } else {
    mdls[7] <- FALSE
    mdls[8] <- FALSE
    casemdlmem <- NULL
    ctrlmdlmem <- NULL
    caseids <- NULL
    ctrlids <- NULL
  }
  return (list(tests = mdls,
               subids = subids,
               gomdlmem = gomdlmem,
               gemdlmem = gemdlmem,
               gxemdlmem = gxemdlmem,
               gxe0mdlmem = gxe0mdlmem,
               egmdlmem = egmdlmem,
               caseids = caseids,
               casemdlmem = casemdlmem,
               ctrlids = ctrlids,
               ctrlmdlmem = ctrlmdlmem))
}

#' Routine to run a GWEIS
#'
#' @param gweismem Models and memory allocated by gweis-mem to run the GWEIS.
#' @param bdinfo Information about a Binary Dosage file that contains the
#' genetic data to run the GWEIS
#' @param snps List of SNPs in the Binary Dosage file to perform the GWEIS on.
#' @param outfilename Name of the file to contain the output
#' @param maf Minimum minor allele frequency of SNPs needed to run test on.
#'
#' @return None
#' @export
rungweis <- function(gweismem, bdinfo, snps, outfilename, maf) {
  if (missing(maf) == TRUE)
    maf <- 0.01
  minaaf <- maf
  maxaaf <- 1. - maf

  xr <- matrix(0, length(gweismem$subids), 1)
  xr2 <- matrix(0, nrow(xr), 2)
  xrgxe <- matrix(0, nrow(xr), 1)

  submatch <- match(gweismem$subids, bdinfo$samples$sid)
  if (length(gweismem$caseids) > 0) {
    casematch <- match(gweismem$caseids, bdinfo$samples$sid)
    xrcase <- matrix(0, length(gweismem$caseids))
  }
  if (length(gweismem$ctrlids) > 0) {
    ctrlmatch <- match(gweismem$ctrlids, bdinfo$samples$sid)
    xrctrl <- matrix(0, length(gweismem$ctrlids))
  }

  statnames <-  c("aaf", "aaf_e0", "aaf_e1",
                  "bg_go", "bg_go_lrt",
                  "bg_ge", "bg_ge_lrt",
                  "bg_gxe", "bg_gxe_lrt",
                  "bgxe", "bgxe_lrt",
                  "joint_lrt",
                  "bg_eg", "bg_eg_lrt",
                  "bg_case", "bg_case_lrt",
                  "bg_ctrl", "bg_ctrl_lrt")
  statsout <- c(TRUE, FALSE, FALSE, gweismem$tests[c(1,1,2,2,3,3,4,4,5,6,6,7,7,8,8)])
  if (gweismem$tests[5] == TRUE) {
    statsout[8] <- TRUE
    statsout[10] <- TRUE
  }
  snpinfo <- paste("snpid", "chr", "loc", "ref", "alt", sep = '\t')
  if (gweismem$test[8] == TRUE)
    statsout[2] <- TRUE
  if (gweismem$test[7] == TRUE)
    statsout[3] <- TRUE

  outline <- paste(statnames[statsout], collapse = "\t")
  outfile <- file(outfilename, "w")
  writeLines(paste(snpinfo, outline, sep = "\t"), outfile)

  numsnps <- length(snps)
  #  if (numsnps > 10)
  #    numsnps <- 10
  outvalues <- numeric(18)
  for (i in 1:numsnps) {
    g <- BinaryDosage::getsnp(bdinfo = bdinfo,
                snp = snps[i])$dosage
    xr[,1] <- g[submatch]
    outvalues[1] <- mean(xr) / 2
    if (outvalues[1] < minaaf | outvalues[1] > maxaaf) {
      next
    }

    if (length(gweismem$gomdlmem) > 0) {
      lsReg::runtest(gweismem$gomdlmem, xr)
    }
    if (length(gweismem$gemdlmem) > 0) {
      lsReg::runtest(gweismem$gemdlmem, xr)
    }
    if (length(gweismem$gxemdlmem) > 0) {
      xr2[,1] <- xr
      xr2[,2] <- xr * gweismem$gxemdlmem$fitdata$xl[,ncol(gweismem$gxemdlmem$fitdata$xl)]
      lsReg::runtest(gweismem$gxemdlmem, xr2)
    }
    if (length(gweismem$gxe0mdlmem) > 0) {
      xrgxe[,1] <- xr2[,2]
      lsReg::runtest(gweismem$gxe0mdlmem, xrgxe)
    }
    if (length(gweismem$egmdlmem) > 0) {
      lsReg::runtest(gweismem$egmdlmem, xr)
    }
    if (length(gweismem$casemdlmem) > 0) {
      xrcase[,1] <- g[casematch]
      outvalues[2] <- mean(xrcase) / 2
      lsReg::runtest(gweismem$casemdlmem, xrcase)
    }
    if (length(gweismem$ctrlmdlmem) > 0) {
      xrctrl[,1] <- g[ctrlmatch]
      outvalues[3] <- mean(xrctrl) / 2
      lsReg::runtest(gweismem$ctrlmdlmem, xrctrl)
    }
    if (gweismem$tests[1] == TRUE) {
      outvalues[4] <- gweismem$gomdlmem$fitdata$betab[1]
      outvalues[5] <- gweismem$gomdlmem$testvalue
    }
    if (gweismem$tests[2] == TRUE) {
      outvalues[6] <- gweismem$gemdlmem$fitdata$betab[1]
      outvalues[7] <- gweismem$gemdlmem$testvalue
    }
    if (gweismem$tests[3] == TRUE) {
      outvalues[8] <- gweismem$gxemdlmem$fitdata$betab[1]
      outvalues[9] <- 2*(gweismem$gxemdlmem$loglike[2] - gweismem$gxe0mdlmem$loglike[2])
    }
    if (gweismem$tests[4] == TRUE) {
      outvalues[10] <- gweismem$gxemdlmem$fitdata$betab[2]
      outvalues[11] <- 2*(gweismem$gxemdlmem$loglike[2] - gweismem$gemdlmem$loglike[2])
    }
    if (gweismem$tests[5] == TRUE) {
      outvalues[8] <- gweismem$gxemdlmem$fitdata$betab[1]
      outvalues[10] <- gweismem$gxemdlmem$fitdata$betab[2]
      outvalues[12] <- gweismem$gxemdlmem$testvalue
    }
    if (gweismem$tests[6] == TRUE) {
      if (gweismem$egmdlmem$fitdata$family == "gaussian")
        outvalues[13] <- gweismem$egmdlmem$fitdata$bb[1]
      else
        outvalues[13] <- gweismem$egmdlmem$fitdata$betab[1]
      outvalues[14] <- gweismem$egmdlmem$testvalue
    }
    if (gweismem$tests[7] == TRUE) {
      if (gweismem$casemdlmem$fitdata$family == "gaussian")
        outvalues[15] <- gweismem$casemdlmem$fitdata$bb[1]
      else
        outvalues[15] <- gweismem$casemdlmem$fitdata$betab[1]
      outvalues[16] <- gweismem$casemdlmem$testvalue
    }
    if (gweismem$tests[8] == TRUE) {
      if (gweismem$egmdlmem$fitdata$family == "gaussian")
        outvalues[17] <- gweismem$ctrlmdlmem$fitdata$bb[1]
      else
        outvalues[17] <- gweismem$ctrlmdlmem$fitdata$betab[1]
      outvalues[18] <- gweismem$ctrlmdlmem$testvalue
    }
    snpinfo <- paste(bdinfo$snps[snps[i], c(3,1,2,4,5)], collapse = "\t")
    writeLines(paste(snpinfo, paste(outvalues[statsout],collapse = "\t"), sep = "\t"), outfile)
  }
  close(outfile)
}
