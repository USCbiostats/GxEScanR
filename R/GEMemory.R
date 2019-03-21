GELinRegMemory <- function(n, p, q) {
  betaT0 <- numeric(p)
  betaT <- numeric(p)
  betaB <- numeric(q)
  qL <- matrix(0., n, p)
  qR <- matrix(0., n, q)
  rTL <- matrix(0., p, p)
  rTR <- matrix(0., p, q)
  rBR <- matrix(0., q, q)
  sigma2 <- numeric(2)
  logLikelihood <- numeric(2)
  return (list(betaT0 = betaT0,
               betaT = betaT,
               betaB = betaB,
               qL = qL,
               qR = qR,
               rTL = rTL,
               rTR = rTR,
               rBR = rBR,
               sigma2 = sigma2,
               logLikelihood = logLikelihood))
}

GELogRegMemory <- function(n, p, q) {
  yp <- numeric(n)
  xlw <- matrix(0., n, p)
  xrw <- matrix(0., n, q)
  zt <- numeric(p)
  zb <- numeric(q)
  k <- numeric(p)
  t <- matrix(0., n, q)
  h <- matrix(0., p, q)
  beta <- numeric(p + q)
  beta0 <- numeric(p)
  beta0diff <- numeric(p)
  betaT <- numeric(p)
  betaB <- numeric(q)
  qL <- matrix(0., n, p)
  qR <- matrix(0., n, q)
  rTL <- matrix(0., p, p)
  rTR <- matrix(0., p, q)
  rBR <- matrix(0., q, q)
  abx <- numeric(n)
  expabx <- numeric(n)
  expabxp1 <- numeric(n)
  expitabx <- numeric(n)
  w <- numeric(n)
  wInv <- numeric(n)
  score <- numeric(p + q)
  loglikelihood <- numeric(2)
  return (list(yp = yp,
               xlw = xlw,
               xrw = xrw,
               zt = zt,
               zb = zb,
               k = k,
               t = t,
               h = h,
               beta0 = beta0,
               beta0diff = beta0diff,
               beta = beta,
               betaT = betaT,
               betaB = betaB,
               qL = qL,
               qR = qR,
               rTL = rTL,
               rTR = rTR,
               rBR = rBR,
               abx = abx,
               expabx = expabx,
               expabxp1 = expabxp1,
               expitabx = expitabx,
               w = w,
               wInv = wInv,
               score = score,
               loglikelihood = loglikelihood))
}

  