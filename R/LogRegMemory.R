# Allocates memory for large scale logistic regression
# does not change from iteration to iteration
AllocateLogRegConstantMemory <- function(x, y) {
  n <- nrow(x)
  p <- ncol(x)
  beta <- numeric(p)
  score <- numeric(p)
  w <- numeric(n)
  wInv <- numeric(n)
  yp <- numeric(n)
  zt <- numeric(p)
  k <- numeric(p)
  ql <- matrix(0., nrow = n, ncol = p)
  rtl <- matrix(0., nrow = p, ncol = p)
  logLikelihood <- numeric(1)
  return (list(n = n,
               p = p,
               beta = beta,
               score = score,
               w = w,
               wInv = wInv,
               yp = yp,
               zt = zt,
               k = k,
               ql = ql,
               rtl = rtl,
               logLikelihood = logLikelihood))
}

# Allocates memory for large scale logistic regression
# that is fixed in size but can vary iteration to iteration
AllocateLogRegFixedMemory <- function(n, p) {
  abx <- numeric(n)
  expabx <- numeric(n)
  expabxp1 <- numeric(n)
  expitabx <- numeric(n)
  yp <- numeric(n)
  zt <- numeric(p)
  k <- numeric(p)
  bt <- numeric(p)
  return (list(abx = abx,
               expabx = expabx,
               expabxp1 = expabxp1,
               expitabx = expitabx,
               yp = yp,
               zt = zt,
               k = k,
               bt = bt))
}

# Allocates memory for large scale logistic regression
# that is not fixed in size and can vary iteration to iteration
AllocateLogRegNotFixedMemory <- function(n, p, q) {
  xrw <- matrix(0., nrow = n, ncol = q)
  beta <- numeric(p + q)
  score <- numeric(p + q)
  zb <- vector(mode = "numeric", length = q)
  bb <- vector(mode = "numeric", length = q)
  h <- matrix(0, nrow = p, ncol = q)
  rtr <- matrix(0, nrow = p, ncol = q)
  t <- matrix(0., nrow = n, ncol = q)
  qr <- matrix(0., nrow = n, ncol = q)
  rbr <- matrix(0., nrow = q, ncol = q)
  logLikelihood <- numeric(1)
  return (list(xrw = xrw,
               beta = beta,
               score = score,
               zb = zb,
               bb = bb,
               h = h,
               rtr = rtr,
               t = t,
               qr = qr,
               rbr = rbr,
               logLikelihood = logLikelihood))
}

# Allocates memory need for large scale logistic regression
AllocateLargeScaleLogRegMemory <- function(y, x, gxe, errorInfo) {
  if (ncol(x) == 1) {
    df <- data.frame(y = y)
  } else {
    df <- data.frame(y = y, x = x[,2:ncol(x)])
  }
  rlogreg <- glm(y ~ ., data = df, family = "binomial")
  
  p1 <- AllocateLogRegConstantMemory(x, y)
  p2 <- AllocateLogRegFixedMemory(p1$n, p1$p)
  p1$beta <- rlogreg$coefficients
  result <- InitializeLargeScaleLogReg(p1$n, p1$p, y, x,
                                       p1$beta, p1$score, p1$w, p1$wInv,
                                       p1$yp, p1$zt, p1$k, p1$ql, p1$rtl,
                                       p2$abx, p2$expabx, p2$expabxp1, p2$expitabx,
                                       p1$logLikelihood)
  if (result != 0)
    stop("Error initializing D|E model")
  if (is.na(p1$logLikelihood) == TRUE) {
    if (errorInfo) {
      GxEErrorData <- list(p1 = p1,
                           p2 = p2,
                           rlogreg = rlogreg,
                           message = "Error initializing D|E model")
      saveRDS(GxEErrorData, "GxEScanErrorData")
    }
    stop("Error initializing D|E model")
  }
  if (abs(p1$logLikelihood - logLik(rlogreg)) > 1e-7)
    stop("Error calculating log likelihood for D|E")
  p3 <- AllocateLogRegNotFixedMemory(p1$n, p1$p, 1)
  if (gxe == TRUE) {
    p3gxe <- AllocateLogRegNotFixedMemory(p1$n, p1$p, 2)
    return (list(p1 = p1,
                 p2 = p2,
                 p3 = p3,
                 p3gxe = p3gxe))
  }
  return (list(p1 = p1,
               p2 = p2,
               p3 = p3))
}