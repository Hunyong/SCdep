dBvNB2.check <- function(x, y, a0, a1, a2, b) {
  h <- function(k, m) {result <- lgamma(a0 + x + y - k - m) - lgamma(a0 + y - m) - lgamma(x - k + 1) +
    lgamma(a0 + y - m) - lgamma(a0) - lgamma(y - m + 1) +
    lgamma(a1 + k) - lgamma(a1) - lgamma(k + 1) +
    lgamma(a2 + m) - lgamma(a2) - lgamma(m + 1) +
    (x + y) * log(b) + (k + m - x - y - a0) * log(1 + 2*b) - (k + a1 + m + a2) * log(1 + b)
  return(exp(result))
  }
  h0 <- function(k, m) {h(k,m) * digamma(x+y-k-m+a0)}
  h1 <- function(k, m) {h(k,m) * digamma(k + a1)}
  h2 <- function(k, m) {h(k,m) * digamma(m + a2)}
  
  result <- sum(sapply(0:y, function(m) sapply(0:x, h, m=m)))
  result0 <- sum(sapply(0:y, function(m) sapply(0:x, h0, m=m)))
  result1 <- sum(sapply(0:y, function(m) sapply(0:x, h1, m=m)))
  result2 <- sum(sapply(0:y, function(m) sapply(0:x, h2, m=m)))
  
  return(list(den=result, h0=result0, h1=result1, h2=result2))
}

BvNB2.Expt <- function (x, y, a0, a1, a2, b) {
  den <- dBvNB2.check(x, y, a0, a1, a2, b)
  h0 <- den$h0
  h1 <- den$h1
  h2 <- den$h2
  den <- den$den
  den0 <- dBvNB2(x, y, a0 + 1, a1, a2, b, log=FALSE)
  den1 <- dBvNB2(x, y, a0, a1 + 1, a2, b, log=FALSE)
  den2 <- dBvNB2(x, y, a0, a1, a2 + 1, b, log=FALSE)
  ER0 <- a0 * b * den0 / den
  ER1 <- a1 * b * den1 / den
  ER2 <- a2 * b * den2 / den
  ELR0 <- h0/den + log(b) - log(2*b+1)
  ELR1 <- h1/den + log(b) - log(b+1)
  ELR2 <- h2/den + log(b) - log(b+1)
  return(c(den=log(den), ER0 = ER0, ER1 = ER1, ER2=ER2, ELR0 = ELR0,  ELR1 = ELR1,  ELR2 = ELR2))
}
BvNB2.Expt.vec <- Vectorize(BvNB2.Expt)
BvNB2.Expt.mean <- function (xvec, yvec, a0, a1, a2, b) {
  xy.reduced <- as.data.frame(table(xvec,yvec))
  names(xy.reduced) <- c("x", "y","freq")
  xy.reduced <- xy.reduced[xy.reduced$freq != 0,]
  xy.reduced$x <- as.numeric(as.character(xy.reduced$x))
  xy.reduced$y <- as.numeric(as.character(xy.reduced$y))
  xy.reduced$freq <- as.numeric(as.character(xy.reduced$freq))
  n <- sum(xy.reduced$freq)
  result <- as.vector(BvNB2.Expt.vec(xy.reduced$x, xy.reduced$y, a0, a1, a2, b) %*% xy.reduced$freq)/n
  result[1] <- result[1] * n
  names(result) <- c("den", "ER0", "ER1", "ER2", "ELR0", "ELR1", "ELR2")
  return(result)
}

ML.BvNB2.EM <- function (xvec, yvec, tol=1e-8, showFlag=FALSE, initial=NULL) {
  xy.reduced <- as.data.frame(table(xvec,yvec))
  names(xy.reduced) <- c("x", "y","freq")
  xy.reduced <- xy.reduced[xy.reduced$freq != 0,]
  xy.reduced$x <- as.numeric(as.character(xy.reduced$x))
  xy.reduced$y <- as.numeric(as.character(xy.reduced$y))
  xy.reduced$freq <- as.numeric(as.character(xy.reduced$freq))
  if (max(xvec)==0 & max(yvec)==0) {return(rep(1e-10,4))}
  
  # initial guess
  if (is.null(initial)) {
    xbar <- mean(xvec); ybar <- mean(yvec); xybar <- mean(c(xbar, ybar))
    s2.x <- var(xvec); s2.y <- var(yvec); if(is.na(s2.x)) {s2.x <- s2.y <- 1}
    s2.mean <- mean(s2.x, s2.y)
    cor.xy <- cor(xvec,yvec); if (is.na(cor.xy)) {cor.xy <- 0}
    initial <- rep(NA,4)
    initial[4] <- s2.mean /xybar
    initial[2:3] <- c(xbar,ybar)/initial[4]
    initial[1] <- min(initial[2],initial[3]) * abs(cor.xy)
    initial[2:3] <-  initial[2:3] - initial[1]
    initial <- pmax(initial, 1e-5)
    #print(initial) ###    
  }
  
  iter = 0
  param = initial
  if (showFlag) {print(c("iter", "a0", "a1", "a2", "b"))}
  repeat {
    iter = iter + 1
    param.old <- param
    a0 = param[1]; a1 = param[2]; a2 = param[3]; b = param[4]
    Eval <- BvNB2.Expt.mean(xvec, yvec, a0, a1, a2, b)
    # opt.b <- function(bb) {sum(Eval[2:4]) - sum(igamma(Eval[5:7] - log(bb))) * bb}
    opt.vec <- function(param) {
      a0 <- param[1]; a1 <- param[2]; a2 <- param[3]; b <- param[4]
      r1 <- sum(Eval[2:4]) - (a0 + a1 + a2) *b
      r2 <- Eval[5:7] - digamma(c(a0,a1,a2)) - log(b)
      return(c(r1,r2))
    }
    param <- multiroot(opt.vec, start=c(a0,a1,a2,b))$root %>%print
    # param[4] <- multiroot(opt.b, start=b)$root[1] %>% print
    # param[1:3] <- igamma(Eval[5:7] - log(param[4]))
  # print(c("this!!!",param)) ##
    if (showFlag) {print(c(iter, round(param,5), Eval[1]))} #lik: lik of previous iteration
    if (max(abs(param - param.old)) <= tol) {
      names(param) <- c("a0", "a1", "a2", "b")
      return(param)
      break
    }
  }
  
}
ML.BvNB2.EM(c(10,1,1),c(10,1,2), showFlag=TRUE)
# 1.19297    0.00272    0.00393    3.47934   -2.61458 
# 1.18967     0.00008     0.00025     3.50169   -12.92425 

ML.BvNB2(c(10,1,1),c(10,1,2), showFlag=TRUE)
lik.BvNB2(c(10,1,1),c(10,1,2), c(1.19297,  0.00272,    0.00393,    3.47934))
lik.BvNB2(c(10,1,1),c(10,1,2), c(1.18971, 0.0001450095, 0.0004408919, 3.501135))
ML.BvZINB3(extractor(1), extractor(3), showFlag=TRUE)
