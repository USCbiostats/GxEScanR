#####################################################
###       Perform Linear Regression GWAS
#####################################################
#
# Routine to do a GWAS with a continuous outcome
#
linreggwas <- function(bdinfo, blkinfo, snps, stddata, subindex,
                       outfile, skipfile, minmaf, base) {
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
  ###       Calcualte the minimum number
  ###       of observed genes
  #####################################################
  if (minmaf == 0)
    minsum <- 2 * nrow(stddata) * 0.05
  else
    minsum <- 2 * nrow(stddata) * minmaf
  maxsum <- 2*nrow(stddata) - minsum
  
  #####################################################
  ###       Do initial regression
  #####################################################
  
  y <- stddata[,1]
  stddata[,1] <- 1.
  
  linreg0 <- initlslinreg(y = y, x = stddata)

  #####################################################
  ###       Write file headers
  #####################################################
  
  if (outfile != "")
    write("SNPID\tbetag\tlrtg", outfile);
  if (skipfile != "")
    write("SNPID\treason", skipfile)
  
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
      # Fit D|G using large scale linear regression
      linregres <- lslinreg(y = y,
                            x = stddata,
                            xr = gxr,
                            ql = linreg0$ql,
                            rtl = linreg0$rtl,
                            k = linreg0$k,
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
        lrt <- 2*(tmploglh[!is.na(tmploglh)] - linreg0$loglike)
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
                   loglike0 = linreg0$loglike)
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
  
  # Return 0 if the results were written to a file
  if (outfile != "")
    return (0)
  
  # Return the results if not output to a file
  results <- data.frame(snp = bdinfo$snps$snpid[!is.na(lrt)],
                        betag = beta[!is.na(lrt),1],
                        lrtg = lrt[!is.na(lrt)],
                        stringsAsFactors = FALSE)
  return (results)
}

#####################################################
###       Perform Linear Regression GWIS
#####################################################
#
# Routine to do a GWIS with a continuous outcome
#
linreggwis <- function(bdinfo, blkinfo, snps, stddata, subindex,
                       outfile, skipfile, minmaf, base, estddev) {
  #####################################################
  ###       Useful values to have
  #####################################################
  nsub <- nrow(stddata)
  nbdsub <- nrow(bdinfo$samples)
  nsnps <- nrow(bdinfo$snps)
  maxdoses <- blkinfo$snpsperblk
  #####################################################
  ###       Initialize memory
  #####################################################
  dosages <- matrix(0, nrow = nbdsub, ncol = maxdoses)
  
  # Memory for D|G
  gxr <- matrix(0, nrow = nsub, ncol = maxdoses)
  tmploglhg <- numeric(maxdoses)
  tmpbetag <- matrix(0, nrow = maxdoses, ncol = 1)
  if (outfile == "") {
    lrtg <- numeric(nsnps)
    betag <- matrix(0, nrow = nsnps, ncol = 1)
  }
  
  # Memory for D|G,E,GxE
  gxexr <- matrix(0, nrow = nsub, ncol = 2*maxdoses)
  tmploglhgxe <- numeric(maxdoses)
  tmpbetagxe <- matrix(0, nrow = maxdoses, ncol = 2)
  if (outfile == "") {
    lrtgxe <- numeric(nsnps)
    lrt2df <- numeric(nsnps)
    betagxe <- matrix(0, nrow = nsnps, ncol = 2)
  }
  
  #####################################################
  ###       Calcualte the minimum number of
  ###       observed genes
  #####################################################
  if (minmaf == 0)
    minsum <- 2 * nsub * 0.05
  else
    minsum <- 2 * nsub * minmaf
  maxsum <- 2 * nsub - minsum
  
  #####################################################
  ###       Do initial regression
  #####################################################

  # Save y and set the intercept column to 1
  y <- stddata[,1]
  stddata[,1] <- 1.
  
  # Initialize large scale linear regression
  linreg0 <- initlslinreg(y = y, x = stddata)
  # Fill cubes with g and gxe
  seqg <- seq(1, 2*maxdoses - 1, 2)
  seqgxe <- seq(2, 2*maxdoses, 2)
  
  #####################################################
  ###       Write file headers
  #####################################################
  
  if (outfile != "")
    write("SNPID\tbetag\tlrtg\tbetagxe\tlrtgxe\tlrt2df", outfile);
  if (skipfile != "")
    write("SNPID\treason", skipfile)
  
  #####################################################
  ###       Loop over blocks
  #####################################################
  numblks <- length(blkinfo$blkloc)
#  numblks <- 1
  blkbuffer <- integer((max(blkinfo$blkbytes) + 3) %/% 4)
  firstsnp <- 1
  subindexm1 <- subindex - 1
  maxn <- blkinfo$snpsperblk
  for (i in 1:numblks) {
    # Set up first and last SNP and adjust memory
    # for last group
    lastsnp <- firstsnp + blkinfo$snpsperblk - 1
    if (lastsnp > nrow(bdinfo$snps)) {
      lastsnp <- nrow(bdinfo$snps)
      maxn <- lastsnp - firstsnp + 1
      tmploglhg <- numeric(maxn)
      tmpbetag <- matrix(0, maxn, 1)
      tmploglhgxe <- numeric(maxn)
      tmpbetagxe <- matrix(0, nrow = maxn, ncol = 2)
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
      # Fit D|G using large scale linear regression
      linregres <- lslinreg(y = y,
                            x = stddata,
                            xr = gxr,
                            ql = linreg0$ql,
                            rtl = linreg0$rtl,
                            k = linreg0$k,
                            q = 1L,
                            skipped = snps,
                            skipoffset = firstsnp - 1,
                            maxn = maxn,
                            minsum = minsum,
                            maxsum = maxsum,
                            loglike = tmploglhg,
                            beta = tmpbetag)
      # Make the xr matrix using the previously calculated one
      # and the x value from the data
      makegxexr(dest = gxexr,
                src1 = gxr,
                src2 = stddata)
      linregres <- lslinreg(y = y,
                            x = stddata,
                            xr = gxexr,
                            ql = linreg0$ql,
                            rtl = linreg0$rtl,
                            k = linreg0$k,
                            q = 2L,
                            skipped = snps,
                            skipoffset = firstsnp - 1,
                            maxn = lastsnp - firstsnp + 1,
                            minsum = minsum,
                            maxsum = maxsum,
                            loglike = tmploglhgxe,
                            beta = tmpbetagxe)
    } else {
      # Set values indicating that the SNPs were not processed
      # because they were not in the list the user provided
      tmploglhg[1:length(tmploglhg)] <- NA
      tmpbetag[1:length(tmpbetag)] <- 1
      tmploglhgxe[1:length(tmploglhgxe)] <- NA
      tmpbetagxe[1:length(tmpbetagxe)] <- 1
    }
    
    if (outfile != "") {
      # Write the results to the given output file
      tokeep <- !is.na(tmploglhg)
      if (length(tokeep) > 0) {
        betag <- tmpbetag[tokeep]
        lrtg <- 2*(tmploglhg[tokeep] - linreg0$loglike)
        lrt2df <- 2*(tmploglhgxe[tokeep] - linreg0$loglike)
        lrtgxe <- lrt2df - lrtg
        betagxe <- tmpbetagxe[tokeep, 2] / estddev
        betagxe[is.na(lrtgxe)] <- NA
        snpids <- bdinfo$snps$snpid[firstsnp:lastsnp]
        snpids <- snpids[tokeep]
        outstring <- paste(snpids, betag, lrtg, betagxe,
                           lrtgxe, lrt2df, sep = '\t')
        write(outstring, outfile, append = TRUE)
      }
    } else {
      calculatelrtgxe(lrtg = lrtg,
                      lrtgxe = lrtgxe,
                      lrt2df = lrt2df,
                      idx1 = firstsnp,
                      idx2 = lastsnp,
                      loglikeg = tmploglhg,
                      loglikegxe = tmploglhgxe,
                      loglike0 = linreg0$loglike)
      copybeta(dest = betag,
               src = tmpbetag,
               colnum = 1,
               startrow = firstsnp,
               numrows = maxn)
      copybeta(dest = betagxe,
               src = tmpbetagxe,
               colnum = 2,
               startrow = firstsnp,
               numrows = maxn)
    }
 
    firstsnp <- lastsnp + 1
  }
  
  if (outfile == "") {
    tokeep <- !is.na(lrtg)
    results <- data.frame(snp = bdinfo$snps$snpid[tokeep],
                          betadg = betag[tokeep],
                          lrtdg = lrtg[tokeep],
                          betagxe = betagxe[tokeep] / estddev,
                          lrtgxe = lrtgxe[tokeep],
                          lrt2df = lrt2df[tokeep],
                          stringsAsFactors = FALSE)
    
    return (results)
  }
  
  return (0)
}
