# Allocates memory for large scale logistic regression
# does not change from iteration to iteration
AllocateLRModConstantMemory <- function(x, y) {
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
  logLikelihood <- 0.
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
AllocateLRModFixedMemory <- function(n, p) {
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
AllocateLRModNotFixedMemory <- function(n, p, q) {
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