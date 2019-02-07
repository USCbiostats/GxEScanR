GELinRegMemory <- function(n, p, q) {
  betaT <- numeric(p)
  betaB <- numeric(q)
  qL <- matrix(0., n, p)
  qR <- matrix(0., n, q)
  rTL <- matrix(0., p, p)
  rTR <- matrix(0., p, q)
  rBR <- matrix(0., q, q)
  sigma2 <- numeric(2)
  logLikelihood <- numeric(2)
  return (list(betaT = betaT,
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
  z <- numeric(p)
  beta <- numeric(p + q)
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
               z = z,
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

  