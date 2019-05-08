# Allocates memory need for large scale logistic regression
AllocateLinRegMemory <- function(y, x) {
  n <- numeric(1)
  n[1] <- length(y)
  p <- numeric(1)
  p[1] <- ncol(x)
  q <- numeric(1)
  q[1] <- 1
  ql <- matrix(0., n, p)
  rtl <- matrix(0., p, p)
  zt <- numeric(p)
  k <- numeric(p)
  rtr <- matrix(0., p, q)
  t <- matrix(0., n, q)
  qr <- matrix(0., n, q)
  rbr <- matrix(0., q, q)
  h <- matrix(0., p, q)
  zb <- numeric(q)
  bb <- numeric(q)
  bt <- numeric(p)
  logLikelihood <- numeric(2)
  return(list(n = n,
              p = p,
              q = q,
              ql = ql,
              qr = qr,
              rtl = rtl,
              rtr = rtr,
              rbr = rbr,
              bt = bt,
              bb = bb,
              zt = zt,
              zb = zb,
              k = k,
              t = t,
              h = h,
              logLikelihood = logLikelihood))
}

# Allocates memory need for large scale logistic regression
AllocateLargeScaleLinRegMemory <- function(y, x, modelName) {
  if (ncol(x) == 1) {
    df <- data.frame(y = y)
  } else {
    df <- data.frame(y = y, x = x[,2:ncol(x)])
  }
  rlinreg <- glm(y ~ ., data = df)
  
  p1 <- AllocateLinRegMemory(y, x)
  p1$bt <- rlinreg$coefficients
  result <- InitializeLargeScaleLinReg(y, x, p1$bt,
                                       p1$ql, p1$rtl,
                                       p1$zt, p1$k, 
                                       p1$logLikelihood)
  if (result != 0)
    stop(paste("Error initializing", modelName, "model"))
  if (abs(p1$logLikelihood[1] - logLik(rlinreg)) > 1e-7)
    stop("Error calculating log likelihood for D|E")
  return (p1)
}