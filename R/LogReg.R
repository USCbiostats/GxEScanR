#####################################################
###       Perform Logistic Regression GWAS
#####################################################
#
# Routine to do a GWAS with a binary outcome
#
logreggwas <- function(bdinfo, blkinfo, snps, stddata, subindex,
                       outfile, skipfile, beta0, minmaf, base) {
  #####################################################
  ###       Initialize memory
  #####################################################
  dosages <- matrix(0,
                    nrow = nrow(bdinfo$samples),
                    ncol = blkinfo$snpsperblk)
  gxr <- matrix(0,
                nrow = nrow(stddata),
                ncol = blkinfo$snpsperblk)
  tmploglh <- numeric(blkinfo$snpsperblk)
  tmpbeta <- matrix(0, blkinfo$snpsperblk, 1)
  if (outfile == "") {
    lrt <- numeric(nrow(bdinfo$snps))
    lrt[1:nrow(bdinfo$snps)] <- NA
    beta <- matrix(0, nrow(bdinfo$snps), 1)
  }
  
  #####################################################
  ###       Calculate the minimum number of
  ###       observed genes
  #####################################################
  if (minmaf == 0)
    minsum <- 2 * nrow(stddata) * 0.05
  else
    minsum <- 2 * nrow(stddata) * minmaf

  if (minsum < 3 * ncol(stddata))
    minsum = 3 * ncol(stddata)

  maxsum <- 2*nrow(stddata) - minsum
  
  #####################################################
  ###       Do initial regression
  #####################################################
  
  y <- stddata[,1]
  stddata[,1] <- 1.
  logreg0 <- initlslogreg(y = y,
                          xl = stddata,
                          beta = beta0)
  
  #####################################################
  ###       Write file headers
  #####################################################
  
  if (outfile != "")
    write("snp\tbetag\tlrtg", outfile);
  if (skipfile != "")
    write("snp\treason", skipfile)
  
  #####################################################
  ###       Loop over blocks
  #####################################################
  firstsnp <- 1
  numblks <- length(blkinfo$blkloc)
#  numblks <- 1
  blkbuffer <- integer((max(blkinfo$blkbytes) + 3) %/% 4)
  subindexm1 <- subindex - 1
  maxn <- blkinfo$snpsperblk
  for (i in 1:numblks) {
    # Set up first and last SNP and adjust memory
    # for last group
    lastsnp <- firstsnp + blkinfo$snpsperblk - 1
    if (lastsnp > nrow(bdinfo$snps)) {
      lastsnp <- nrow(bdinfo$snps)
      maxn <- lastsnp - firstsnp + 1
      tmploglh <- numeric(maxn)
      tmpbeta <- matrix(0, maxn, 1)
    }
    # Do the regression if there are SNPs in the block
    # that are not skipped
    if (any(snps[firstsnp:lastsnp] == TRUE)) {
      # Read in the dosages
      readblock(filename = bdinfo$filename,
                blkloc = blkinfo$blkloc[i],
                blkbytes = blkinfo$blkbytes[i],
                blkbuffer = blkbuffer)
      # Calculate the dosages
      getdosages(dosages = dosages,
                 blkbuffer = blkbuffer,
                 fileloc = blkinfo$blkloc[i],
                 indices = bdinfo$indices,
                 firstsnp = firstsnp,
                 lastsnp = lastsnp,
                 base = base)
      # Make the xr matrix by subsetting out the SNPs
      # for the selected subjects
      makegxr(dest = gxr,
              src = dosages,
              idx = subindexm1)
      # Fit D|G
      logregres <- lslogreg(y = y,
                            xl = stddata,
                            xr = gxr,
                            beta0 = beta0,
                            yp0 = logreg0$yp,
                            ql = logreg0$ql,
                            rtl = logreg0$rtl,
                            k0 = logreg0$k,
                            w = logreg0$w,
                            winv = logreg0$winv,
                            q = 1L,
                            skipped = snps,
                            skipoffset = firstsnp - 1,
                            maxn = maxn,
                            minsum = minsum,
                            maxsum = maxsum,
                            loglike = tmploglh,
                            beta = tmpbeta)
    } else {
      # Set values indicating that the SNPs were not processed
      # because they were not in the list the user provided
      tmploglh[1:length(tmploglh)] <- NA
      tmpbeta[1:length(tmpbeta)] <- 1
    }

    if (outfile != "") {
      # Write the results to the given output file
      beta <- tmpbeta[!is.na(tmploglh)]
      if (length(beta) > 0) {
        lrt <- 2*(tmploglh[!is.na(tmploglh)] - logreg0$loglike)
        snpids <- bdinfo$snps$snpid[firstsnp:lastsnp]
        snpids <- snpids[!is.na(tmploglh)]
        outstring <- paste(snpids, beta, lrt, sep = '\t')
        write(outstring, outfile, append = TRUE)
      }
    } else {
      # Save the results in R
      calculatelrt(lrt = lrt,
                   idx1 = firstsnp,
                   idx2 = lastsnp,
                   loglike = tmploglh,
                   loglike0 = logreg0$loglike)
      copybeta(dest = beta,
               src = tmpbeta,
               colnum = 1,
               startrow = firstsnp,
               numrows = maxn)
    }
    
    if (skipfile != "") {
      # Output the reason SNPs were skipped
      reason <- as.integer(tmpbeta[is.na(tmploglh)])
      snpids <- bdinfo$snps$snpid[firstsnp:lastsnp]
      snpids <- snpids[is.na(tmploglh)]
      if (length(snpids) > 0) {
        outstring <- paste(snpids, reason, sep = '\t')
        write(outstring, skipfile, append = TRUE)
      }
    }
    
    firstsnp <- lastsnp + 1
  }
  # Return 0 if results were written to file
  if (outfile != "")
    return(0)
  
  # Return results as data frame
  results <- data.frame(snp = bdinfo$snps$snpid[!is.na(lrt)],
                        betag = beta[!is.na(lrt),1],
                        lrtg = lrt[!is.na(lrt)],
                        stringsAsFactors = FALSE)
  return (results)
}

b0logreg <- function(x, ebinary) {
  beta0 <- vector("list", 5)
  modfamily = "binomial"
  for (i in 1:5) {
    df <- data.frame(x[[i]])
    colnames(df)[1] <- "X1"
    if (i == 3) {
      if (ebinary == FALSE)
        modfamily <- "gaussian"
    }
    modformula <- "X1 ~ ."
    basemodel <- glm(formula = modformula,
                     family = modfamily,
                     data = df)
    beta0[[i]] <- basemodel$coefficients
  }

  return (beta0)
}

initlogreggweis <- function(y, x, beta, ebinary) {
  yp0 <- vector("list", 5)
  ql <- vector("list", 5)
  rtl <- vector("list", 5)
  k0 <- vector("list", 5)
  w <- vector("list", 5)
  winv <- vector("list", 5)
  loglike0 <- numeric(5)
  for (i in 1:2) {
    initvalues <- initlslogreg(y = y[[i]],
                               xl = x[[i]],
                               beta = beta[[i]])
    yp0[[i]] <- initvalues$yp
    ql[[i]] <- initvalues$ql
    rtl[[i]] <- initvalues$rtl
    k0[[i]] <- initvalues$k
    w[[i]] <- initvalues$w
    winv[[i]] <- initvalues$winv
    loglike0[i] <- initvalues$loglike
  }
  if (ebinary == TRUE) {
    for (i in 3:5) {
      initvalues <- initlslogreg(y = y[[i]],
                                 xl = x[[i]],
                                 beta = beta[[i]])
      yp0[[i]] <- initvalues$yp
      ql[[i]] <- initvalues$ql
      rtl[[i]] <- initvalues$rtl
      k0[[i]] <- initvalues$k
      w[[i]] <- initvalues$w
      winv[[i]] <- initvalues$winv
      loglike0[i] <- initvalues$loglike
    }
  } else {
    for (i in 3:5) {
      initvalues <- initlslinreg(y = y[[i]],
                                 x = x[[i]])
      ql[[i]] <- initvalues$ql
      rtl[[i]] <- initvalues$rtl
      k0[[i]] <- initvalues$k
      loglike0[i] <- initvalues$loglike
    }
  }
  return (list(yp0 = yp0,
               ql = ql,
               rtl = rtl,
               k0 = k0,
               w = w,
               winv = winv,
               loglike0 = loglike0))
}

minmafgweis <- function(datalist, minmaf) {
  minsum <- numeric(5)
  maxsum <- numeric(5)
  for (i in 1:5) {
    if (minmaf == 0)
      minsum[i] <- 2 * nrow(datalist[[i]]) * 0.05
    else
      minsum[i] <- 2 * nrow(datalist[[i]]) * minmaf
    
    if (minsum[i] < 3 * ncol(datalist[[i]]))
      minsum[i] <- 3 * ncol(datalist[[i]])
    
    maxsum[i] <- 2 * nrow(datalist[[i]]) - minsum[i]
  }
  return (list(minsum = minsum,
               maxsum = maxsum))
}
#####################################################
###       Perform Logistic Regression GWEIS
#####################################################
#
# Routine to do a GWEIS with a binary outcome
#
logreggweis <- function(bdinfo, blkinfo, snps, stddata, subindex,
                        outfile, skipfile, minmaf, base, e) {
  #####################################################
  ###       Set up needed arrays
  #####################################################
  # Determine if e is binary
  evalues <- unique(e)
  ebinary <- FALSE
  if (length(evalues) == 2) {
    evalues <- sort(evalues)
    if (evalues[1] == 0 & evalues[2] == 1) {
      ebinary <- TRUE
    }
  }
  
  # Set up the x values that don't change, xl
  x <- vector("list", length = 5)
  x[[1]] <- stddata
  x[[2]] <- stddata
  x[[3]] <- as.matrix(stddata[,1:(ncol(stddata)-1)])
  x[[3]][,1] <- e
  x[[4]] <- as.matrix(stddata[stddata[,1] == 1, 1:(ncol(stddata)-1)])
  x[[4]][,1] <- e[stddata[,1] == 1]
  x[[5]] <- as.matrix(stddata[stddata[,1] == 0, 1:(ncol(stddata)-1)])
  x[[5]][,1] <- e[stddata[,1] == 0]
  # The first column was set to y to use the glm routine to
  # get initial values for beta. This is done because it is
  # already written and produces good warnings
  beta0 <- b0logreg(x, ebinary)
  # The first column is changed from the outcome to 1
  # The first column is now the intercept
  for (i in 1:5)
    x[[i]][,1] <- 1;

  # Set up the vector of outcomes, y  
  y <- vector("list", length = 5)
  y[[1]] <- stddata[,1] 
  y[[2]] <- stddata[,1] 
  y[[3]] <- e
  y[[4]] <- e[stddata[,1] == 1]
  y[[5]] <- e[stddata[,1] == 0]

  # Save the needed values from the initial regressions
  # that are needed to do the large scale regression part
  initvalues <- initlogreggweis(y = y,
                                x = x,
                                beta = beta0,
                                ebinary = ebinary)
  
  # The number of columns added
  q <- c(1L, 2L, 1L, 1L, 1L)
  # Calculate the minimum and maximum number of SNPs that
  # can be observed
  minmaxg <- minmafgweis(x, minmaf)
  
  # Indices for cases and controls in binary dosage file
  caseindex <- subindex[stddata[,1] == 1]
  controlindex <- subindex[stddata[,1] == 0]
  
  #####################################################
  ###       Initialize memory
  #####################################################
  # Buffer to hold all dosage data
  dosages <- matrix(0,
                    nrow = nrow(bdinfo$samples),
                    ncol = blkinfo$snpsperblk)
  xr <- vector("list", length = 5)
  xr[[1]] <- matrix(0,
                    nrow = nrow(x[[1]]),
                    ncol = blkinfo$snpsperblk)
  xr[[2]] <- matrix(0,
                    nrow = nrow(x[[2]]),
                    ncol = 2 * blkinfo$snpsperblk)
  xr[[3]] <- matrix(0,
                    nrow = nrow(x[[3]]),
                    ncol = blkinfo$snpsperblk)
  xr[[4]] <- matrix(0,
                    nrow = nrow(x[[4]]),
                    ncol = blkinfo$snpsperblk)
  xr[[5]] <- matrix(0,
                    nrow = nrow(x[[5]]),
                    ncol = blkinfo$snpsperblk)

  tmploglh <- vector("list", length = 5)  
  tmpbeta <- vector("list", length = 5)
  for (i in 1:5) {
    tmploglh[[i]] <- numeric(blkinfo$snpsperblk)
    tmpbeta[[i]] <- matrix(0, blkinfo$snpsperblk, q[i])
  }
  
  if (outfile == "") {
    betag <- numeric(nrow(bdinfo$snps))
    lrtg <- numeric(nrow(bdinfo$snps))
    lrtg[1:nrow(bdinfo$snps)] <- NA
    betagxe <- numeric(nrow(bdinfo$snps))
    lrtgxe <- numeric(nrow(bdinfo$snps))
    lrtgxe[1:nrow(bdinfo$snps)] <- NA
    lrt2df <- numeric(nrow(bdinfo$snps))
    lrt2df[1:nrow(bdinfo$snps)] <- NA
    betaeg <- numeric(nrow(bdinfo$snps))
    lrteg <- numeric(nrow(bdinfo$snps))
    lrteg[1:nrow(bdinfo$snps)] <- NA
    lrt3df <- numeric(nrow(bdinfo$snps))
    lrt3df[1:nrow(bdinfo$snps)] <- NA
    betacase <- numeric(nrow(bdinfo$snps))
    lrtcase <- numeric(nrow(bdinfo$snps))
    lrtcase[1:nrow(bdinfo$snps)] <- NA
    betactrl <- numeric(nrow(bdinfo$snps))
    lrtctrl <- numeric(nrow(bdinfo$snps))
    lrtctrl[1:nrow(bdinfo$snps)] <- NA
  }
  
  #####################################################
  ###       Write file headers
  #####################################################
  
  if (outfile != "")
    write("snp\tbetadg\tlrtdg\tbetagxe\tlrtgxe\tlrt2df\tbetaeg\tlrteg\tlrt3df\tbetacase\tlrtcase\tbetactrl\tlrtctrl", outfile);
  if (skipfile != "")
    write("snp\treasondg\treasongxe\treasoneg\treasoncase\treasonctrl", skipfile)
  
  #####################################################
  ###       Loop over blocks
  #####################################################
  firstsnp <- 1
  numblks <- length(blkinfo$blkloc)
  # numblks <- 1
  blkbuffer <- integer((max(blkinfo$blkbytes) + 3) %/% 4)
  subindexm1 <- subindex - 1
  estd <- stddata[,ncol(stddata)]
  estddev <- sqrt(var(e))
  maxn <- blkinfo$snpsperblk
  for (i in 1:numblks) {
    # Set up first and last SNP and adjust memory
    # for last group
    lastsnp <- firstsnp + blkinfo$snpsperblk - 1
    if (lastsnp > nrow(bdinfo$snps)) {
      lastsnp <- nrow(bdinfo$snps)
      maxn <- lastsnp - firstsnp + 1
      rm(list=(c("tmploglh","tmpbeta")))
      gc()
      tmploglh <- vector("list", length = 5)  
      tmpbeta <- vector("list", length = 5)
      for (j in 1:5) {
        tmploglh[[j]] <- numeric(maxn)
        tmpbeta[[j]] <- matrix(0, maxn, q[j])
      }
    }
    # Do the regression if there are SNPs in the block
    # that are not skipped
    if (any(snps[firstsnp:lastsnp] == TRUE)) {
      # Read in the dosages
      readblock(filename = bdinfo$filename,
                blkloc = blkinfo$blkloc[i],
                blkbytes = blkinfo$blkbytes[i],
                blkbuffer = blkbuffer)
      # Calculate the dosages
      getdosages(dosages = dosages,
                 blkbuffer = blkbuffer,
                 fileloc = blkinfo$blkloc[i],
                 indices = bdinfo$indices,
                 firstsnp = firstsnp,
                 lastsnp = lastsnp,
                 base = base)
      xrgweis2(xr1 = xr[[1]],
               xr2 = xr[[2]],
               xr3 = xr[[3]],
               xr4 = xr[[4]],
               xr5 = xr[[5]],
               idx = subindexm1,
               src1 = dosages,
               src2 = estd)
      for (j in 1:5) {
        if (j > 2 & ebinary == FALSE) {
          linregres <- lslinreg(y = y[[j]],
                                x = x[[j]],
                                xr = xr[[j]],
                                ql = initvalues$ql[[j]],
                                rtl = initvalues$rtl[[j]],
                                k = initvalues$k[[j]],
                                q = q[j],
                                skipped = snps,
                                skipoffset = firstsnp - 1,
                                maxn = maxn,
                                minsum = minmaxg$minsum[j],
                                maxsum = minmaxg$maxsum[j],
                                loglike = tmploglh[[j]],
                                beta = tmpbeta[[j]])
        } else {
          logregres <- lslogreg(y = y[[j]],
                                xl = x[[j]],
                                xr = xr[[j]],
                                beta0 = beta0[[j]],
                                yp0 = initvalues$yp0[[j]],
                                ql = initvalues$ql[[j]],
                                rtl = initvalues$rtl[[j]],
                                k0 = initvalues$k0[[j]],
                                w = initvalues$w[[j]],
                                winv = initvalues$winv[[j]],
                                q = q[j],
                                skipped = snps,
                                skipoffset = firstsnp - 1,
                                maxn = maxn,
                                minsum = minmaxg$minsum[j],
                                maxsum = minmaxg$maxsum[j],
                                loglike = tmploglh[[j]],
                                beta = tmpbeta[[j]])
        }
      }
    } else {
      for (j in 1:5) {
        tmploglh[[j]][1:length(tmploglh[[j]])] <- NA
        tmpbeta[[j]][1:nrow(tmpbeta[[j]]), 1] <- 1
      }
    }
    if (outfile == "") {
      lrtgweis2(lrtg = lrtg,
                lrtgxe = lrtgxe,
                lrt2df = lrt2df,
                lrteg = lrteg,
                lrt3df = lrt3df,
                lrtcase = lrtcase,
                lrtctrl = lrtctrl,
                loglike0 = initvalues$loglike0,
                loglhg = tmploglh[[1]],
                loglhgxe = tmploglh[[2]],
                loglheg = tmploglh[[3]],
                loglhcase = tmploglh[[4]],
                loglhctrl = tmploglh[[5]],
                offset = firstsnp)
      betagweis2(betag = betag,
                 betagxe = betagxe,
                 betaeg = betaeg,
                 betacase = betacase,
                 betactrl = betactrl,
                 tmpbetag = tmpbeta[[1]],
                 tmpbetagxe = tmpbeta[[2]],
                 tmpbetaeg = tmpbeta[[3]],
                 tmpbetacase = tmpbeta[[4]],
                 tmpbetactrl = tmpbeta[[5]],
                 estddev = estddev,
                 offset = firstsnp)
    } else {
      tokeep <- !is.na(tmploglh[[1]])
      if (length(tokeep) > 0) {
        snpids <- bdinfo$snps$snpid[firstsnp:lastsnp]
        snpids <- snpids[tokeep]
        betag <- tmpbeta[[1]][tokeep,1]
        lrtg <- 2*(tmploglh[[1]][tokeep] - initvalues$loglike0[1])
        betagxe <- tmpbeta[[2]][tokeep,2] / estddev
        lrt2df <- 2*(tmploglh[[2]][tokeep] - initvalues$loglike0[2])
        lrtgxe <- lrt2df - lrtg
        betaeg <- tmpbeta[[3]][tokeep,1]
        lrteg <- 2*(tmploglh[[3]][tokeep] - initvalues$loglike0[3])
        lrt3df <- lrt2df + lrteg
        betacase <- tmpbeta[[4]][tokeep,1]
        lrtcase <- 2*(tmploglh[[4]][tokeep] - initvalues$loglike0[4])
        betactrl <- tmpbeta[[5]][tokeep,1]
        lrtctrl <- 2*(tmploglh[[5]][tokeep] - initvalues$loglike0[5])
        outstring <- paste(snpids, betag, lrtg, betagxe, lrtgxe, lrt2df,
                           betaeg, lrteg, lrt3df, betacase, lrtcase,
                           betactrl, lrtctrl, sep = '\t')
        write(outstring, outfile, append = TRUE)
      }
    }
    if (skipfile != "") {
      # Output the reason SNPs were skipped
      tokeep <- is.na(tmploglh[[1]]) | is.na(tmploglh[[2]]) | 
        is.na(tmploglh[[3]]) | is.na(tmploglh[[4]]) | is.na(tmploglh[[5]])
      reason1 <- as.integer(tmpbeta[[1]][tokeep])
      reason2 <- as.integer(tmpbeta[[2]][tokeep, 1])
      reason3 <- as.integer(tmpbeta[[3]][tokeep])
      reason4 <- as.integer(tmpbeta[[4]][tokeep])
      reason5 <- as.integer(tmpbeta[[5]][tokeep])
      snpids <- bdinfo$snps$snpid[firstsnp:lastsnp]
      snpids <- snpids[tokeep]
      if (length(snpids) > 0) {
        outstring <- paste(snpids, reason1, reason2, reason3,
                           reason4, reason5, sep = '\t')
        write(outstring, skipfile, append = TRUE)
      }
    }
    firstsnp <- lastsnp + 1
  }
  
  if (outfile == "") {
    tokeep <- !is.na(lrtg)
    results <- data.frame(snp = bdinfo$snps$snpid[tokeep],
                          betadg = betag[tokeep],
                          lrtdg = lrtg[tokeep],
                          betagxe = betagxe[tokeep],
                          lrtgxe = lrtgxe[tokeep],
                          lrt2df = lrt2df[tokeep],
                          betaeg = betaeg[tokeep],
                          lrteg = lrteg[tokeep],
                          lrt3df = lrt3df[tokeep],
                          betacase = betacase[tokeep],
                          lrtcase = lrtcase[tokeep],
                          betactrl = betactrl[tokeep],
                          lrtctrl = lrtctrl[tokeep])
  } else {
    results <- 0
  }

  return (results)
}