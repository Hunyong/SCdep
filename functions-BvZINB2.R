### This is a zero-inflated version of BvNB2
# Needs dBvNB2 function
library(distr)

### 1. Density, likelihood
dBvZINB2 <- function(x, y, a0, a1, a2, b, pp, log=FALSE) {
  result <- dBvNB2(x=x, y=y, a0=a0, a1=a1, a2=a2, b=b, log=FALSE) * (1-pp) + ifelse(x + y == 0, pp, 0)
  return(ifelse(log, log(result), result))
}
# for cross-check
dBvZINB2.check <- function(x, y, a0, a1, a2, b, pp, log=FALSE) {
  h <- function(k, m) {result <- lgamma(a0 + x + y - k - m) - lgamma(a0 + y - m) - lgamma(x - k + 1) +
      lgamma(a0 + y - m) - lgamma(a0) - lgamma(y - m + 1) +
      lgamma(a1 + k) - lgamma(a1) - lgamma(k + 1) +
      lgamma(a2 + m) - lgamma(a2) - lgamma(m + 1) +
      (x + y) * log(b) + (k + m - x - y - a0) * log(1 + 2*b) - (k + a1 + m + a2) * log(1 + b)
      return(exp(result))
  }
  result <- sum(sapply(0:y, function(m) sapply(0:x, h, m = m)))
  return(ifelse(log, log(result), result))
}
if (FALSE) {
  dBvNB2(1,1,1,1,1,1)
  dBvZINB2(0,1,2,1,3,1,.5) ;dBvNB2(0,1,2,1,3,1) 
  tmp <- sapply(0:50, function(r) sapply(0:50, function(s) dBvZINB2(s,r,1,1,1,1,.5)))
  sum(tmp) #1
}
dBvZINB2.vec <- Vectorize(dBvZINB2)
lik.BvZINB2 <- function(x, y, param) {
  sum(log(dBvZINB2.vec(x, y, param[1], param[2], param[3], param[4], param[5])))
}


### 2.EM

dBvZINB2.Expt <- function(x, y, a0, a1, a2, b, pp) {
  p1 = (2 * b + 1) /(b + 1)
  f00 = (b + 1)^-(a1 + a2) * (2* b + 1)^-a0
  if (x == 0 & y == 0) {
    den <- f00 *(1-pp) + pp
    rate1 <- f00 *(1-pp) / den
    rate2 <- pp / den
    xi0.E <- (a0 * b) * (rate1 /(2*b + 1) + rate2)
    xi1.E <- (a1 * b) * (rate1 /(b + 1) + rate2)
    xi2.E <- (a2 * b) * (rate1 /(b + 1) + rate2)
    z.E   <- rate2
    v.E   <- (a0 + a2) * b * rate2
    log.xi0.E <- digamma(a0) + log(b) - log(1 + 2*b) * rate1
    log.xi1.E <- digamma(a1) + log(b) - log(1 + b) * rate1
    log.xi2.E <- digamma(a2) + log(b) - log(1 + b) * rate1
    result <- c(logdensity = log(den), xi0.E = xi0.E, xi1.E = xi1.E, xi2.E = xi2.E, 
                log.xi0.E = log.xi0.E, log.xi1.E = log.xi1.E, log.xi2.E = log.xi2.E, z.E = z.E, v.E = v.E)
    return(result)
  } else {
    xi0.E1 <- function(k, m) {x - k + y - m + a0}; xi0.E2 <- (1 + 2*b)/b
    xi1.E1 <- function(k) {k + a1}; xi1.E2 <- (1 + b)/b
    xi2.E1 <- function(m) {m + a2}; xi2.E2 <- (1 + b)/b
    log.xi0.E1 <- function(k, m) {digamma(x - k + y - m + a0)}; # log.xi0.E2 <- log(xi0.E2)
    log.xi1.E1 <- function(k) {digamma(k + a1)}; # log.xi1.E2 <- log(xi1.E2)
    log.xi2.E1 <- function(m) {digamma(m + a2)}; # log.xi2.E2 <- log(xi2.E2)
    
    l1 <- function(k, m) exp(lgamma(a1 + k) - lgamma(k+1) - lgamma(a1) + lgamma(x + y + a0 -m -k) - lgamma(x -k +1) - lgamma(a0 + y - m) + k *log(p1)
                             + lgamma(m + a2) - lgamma(m+1) - lgamma(a2) + lgamma(y +a0 -m) - lgamma(y -m +1) - lgamma(a0) + m *log(p1))
    l2 <- - (+x+y+a0)*log(1 + 2* b) + (x +y) * log(b) - (a1 + a2) * log(1 + b)
    l2 <- exp(l2)
    
    l1.mat <- sapply(0:x, function(k) sapply(0:y, l1, k=k))
    l1.mat <- (l1.mat * l2) #%>%print
    l1.sum <- sum(l1.mat) 
    
    xi0.mat <- sapply(0:x, function(k) sapply(0:y, xi0.E1, k=k))
    xi0.mat <- l1.mat * xi0.mat
    xi0.E <- sum(xi0.mat)/sum(l1.mat)/xi0.E2
    # if (xy0) {xi0.E <- (xi0.E * f00 * (1-pp) + a1 * b1 * pp) / (f00 *(1 - pp) + pp)}
    
    xi1.mat <- matrix(sapply(0:x, xi1.E1), x+1, y+1) #%>% print
    xi1.mat <- l1.mat * t(xi1.mat)
    xi1.E <- sum(xi1.mat)/sum(l1.mat)/xi1.E2
    
    xi2.mat <- matrix(sapply(0:y, xi2.E1), y+1, x+1) #%>% print
    xi2.mat <- l1.mat * xi2.mat
    xi2.E <- sum(xi2.mat)/sum(l1.mat)/xi2.E2
    
    log.xi0.mat <- sapply(0:x, function(k) sapply(0:y, log.xi0.E1, k=k))
    log.xi0.mat <- l1.mat * log.xi0.mat
    log.xi0.E <- sum(log.xi0.mat)/sum(l1.mat) - log(xi0.E2)
    
    log.xi1.mat <- matrix(sapply(0:x, log.xi1.E1), x+1, y+1) #%>% print
    log.xi1.mat <- l1.mat * t(log.xi1.mat)
    log.xi1.E <- sum(log.xi1.mat)/sum(l1.mat) - log(xi1.E2)
    
    log.xi2.mat <- matrix(sapply(0:y, log.xi2.E1), y+1, x+1) #%>% print
    log.xi2.mat <- l1.mat * log.xi2.mat
    log.xi2.E <- sum(log.xi2.mat)/sum(l1.mat) - log(xi2.E2)
    
    z.E <- 0
    v.E <- y
    
    ###### temporary #### (simpler but slower)
    ##  den <- dBvZINB2(x,y, a0, a1, a2, b, pp)
    ##  xi0.E <- dBvZINB2(x,y, a0+1, a1, a2, b, pp) *a0 *b/ den
    ##  xi1.E <- dBvZINB2(x,y, a0, a1+1, a2, b, pp) *a1 *b/ den
    ##  xi2.E <- dBvZINB2(x,y, a0, a1, a2+1, b, pp) *a2 *b/ den
    #####################
    
    result <- c(logdensity = log(l1.sum), xi0.E = xi0.E, xi1.E = xi1.E, xi2.E = xi2.E, 
                log.xi0.E = log.xi0.E, log.xi1.E = log.xi1.E, log.xi2.E = log.xi2.E, z.E = z.E, v.E = v.E)
    return(result)    
  }
}
dBvZINB2.Expt.vec <- Vectorize(dBvZINB2.Expt)

dBvZINB2.Expt.Z <- function(x, y, a0, a1, a2, b, pp) {
  f00 = (b + 1)^-(a1 + a2) * (2* b + 1)^-a0
  z.E <- pp / (f00 *(1-pp) + pp) * (x == 0 & y == 0)
  return(z.E)
}
dBvZINB2.Expt.Z.vec <- Vectorize(dBvZINB2.Expt.Z)

dBvZINB2.Expt2 <- function(x, y, a0, a1, a2, b, pp) {
  den <- dBvZINB2(x,y, a0, a1, a2, b, pp)
  xi0.E <- dBvZINB2(x,y, a0+1, a1, a2, b, pp) *a0 *b/ den
  xi1.E <- dBvZINB2(x,y, a0, a1+1, a2, b, pp) *a1 *b/ den
  xi2.E <- dBvZINB2(x,y, a0, a1, a2+1, b, pp) *a2 *b/ den
  result <- c(logdensity = log(den), xi0.E = xi0.E, xi1.E = xi1.E, xi2.E = xi2.E)
  return(result)    
}
dBvZINB2.Expt2.vec <- Vectorize(dBvZINB2.Expt2)

if (FALSE) {
  tmp <- dBvZINB2.Expt.vec(c(1,1,1),c(0,1,2),1,1,1,1,.5)
  tmp <- dBvZINB2.Expt.vec(c(0,1,1),c(0,1,2),1,1,1,1,.5)
  tmp <- dBvZINB2.Expt.vec(extractor(1),extractor(2),1,1,1,1,.5)
  t(tmp)[21:40,]
  apply(dBvZINB2.Expt.vec(c(10,1,2),c(10,1,1), 1.18971, 0.0001450095, 0.0004408919, 3.501135, 0),1,mean)
  dBvZINB2.Expt.vec(c(10,1,1),c(10,1,2), 1.193013282, 0.002745513, 0.003336139, 3.341625901, 0.000000000)
  
  # the above two should have the symmetric results!
  dBvZINB2.Expt(10, 1, 1,2,1,2,0)
  dBvZINB2.Expt(1, 10, 1,1,2,2,0)
}
#EM1 not converging
ML.BvZINB2 <- function (xvec, yvec, initial = NULL, tol=1e-8, showFlag=FALSE) {
  xy.reduced <- as.data.frame(table(xvec,yvec))
  names(xy.reduced) <- c("x", "y","freq")
  xy.reduced <- xy.reduced[xy.reduced$freq != 0,]
  xy.reduced$x <- as.numeric(as.character(xy.reduced$x))
  xy.reduced$y <- as.numeric(as.character(xy.reduced$y))
  xy.reduced$freq <- as.numeric(as.character(xy.reduced$freq))
  n <- sum(xy.reduced$freq)
  if (max(xvec)==0 & max(yvec)==0) {return(c(rep(1e-10,5),1))}
  #print(xy.reduced)
  
  # initial guess
  if (is.null(initial)) {
    xbar <- mean(xvec); ybar <- mean(yvec); xybar <- mean(c(xbar, ybar))
    s2.x <- var(xvec); s2.y <- var(yvec); if(is.na(s2.x)) {s2.x <- s2.y <- 1}
    s2.mean <- mean(s2.x, s2.y)
    cor.xy <- cor(xvec,yvec); if (is.na(cor.xy)) {cor.xy <- 0}
    zero <- sum(xvec == 0 & yvec == 0) / n
    initial <- rep(NA,5)
    initial[4] <- s2.mean /xybar
    initial[5] <- zero
    initial[2:3] <- c(xbar,ybar)/initial[4]
    initial[1] <- min(initial[2],initial[3]) * abs(cor.xy)
    initial[2:3] <-  initial[2:3] - initial[1]
    initial <- pmax(initial, 1e-5)
  }
  
  iter = 0
  param = initial
  if (showFlag) {print(c("iter", "a0", "a1", "a2", "b", "pi"))}
  repeat {
    iter = iter + 1
    #print(c(param))
    #print(lik(vec, pp=param[1], m0=param[2], m1=param[3], m2=param[4])) # debug
    param.old <- param # saving old parameters
    
    # updating
    expt <- dBvZINB2.Expt.vec(xy.reduced$x, xy.reduced$y, 
                              a0 = param[1], a1 = param[2], a2 = param[3], b = param[4], pp = param[5])
    expt <- as.vector(expt %*% xy.reduced$freq / n)
    # loglik = expt[1] * n
    # delta <- expt[9] / (expt[2] + expt[4])                   # delta = E(V) / (E(xi0 + xi2))
                                                        # pi = E(Z)
    
    if (FALSE) {
      # This is numerically instable: imprecise solution does not guarantee maximization
      b1.opt <- function (bb) {(sum(expt[2:4]) - sum (igamma(-log(bb) + expt[5:7])) * bb)} #igamma: inverse-digamma
      param[4]    <- multiroot(b1.opt, start = param[4])$root[1]            # b1
      # multiroot(b1.opt, start = param[4])$root %>% print  # b1 #####
      param[1:3]  <- igamma(-log(param[4]) + expt[5:7])                     # a1 ~ a3
      param[5] = expt[8]
    } else {
      #print(param)##
      opt.vec <- function(par.ab) {
        par.ab <- exp(par.ab)
        r1 <- sum(expt[2:4]) - sum(par.ab[1:3]) * par.ab[4]
        r2 <- expt[5:7] - digamma(par.ab[1:3]) - log(par.ab[4])
      # print(c(r1,r2)) ###
        return(c(r1,r2))
      }
      # param[1:4] <- multiroot(opt.vec, start=param[1:4], rtol=1e-20)$root %>%print
      #multiroot(opt.vec, start=param[1:4], rtol=1e-20)$f.root %>%print ####
      param.l <- log(param)
      result <- try(multiroot(opt.vec, start=param.l[1:4], rtol=1e-20)$root, silent=TRUE)
      if (class(result)=="try-error") {
        initial = rep(1,5)
        result <- multiroot(opt.vec, start=param[1:4], rtol=1e-20)$root
      }
      param[1:4] <- exp(result)
      param[5] = expt[8]
    }
    
  #print (expt) #####
    if (showFlag) {print(c(iter, round(param,5), expt[1] * n))} #lik: lik of previous iteration
    if (max(abs(param - param.old)) <= tol) {
      names(param) <- c("a0", "a1", "a2", "b", "pi")
      return(param)
      break
    }
  }
  #result <- data.frame(a0 = param[1], a1 = param[2], a2 = param[3], b1 = param[4], b2 = param[5], pi = param[6])
  #return(result)
}

#EM2 not converging: something must be wrong with the expectation step!!!
ML.BvZINB2.EM2 <- function (xvec, yvec, initial = NULL, tol=1e-8, showFlag=FALSE) {
  xy.reduced <- as.data.frame(table(xvec,yvec))
  names(xy.reduced) <- c("x", "y","freq")
  xy.reduced <- xy.reduced[xy.reduced$freq != 0,]
  xy.reduced$x <- as.numeric(as.character(xy.reduced$x))
  xy.reduced$y <- as.numeric(as.character(xy.reduced$y))
  xy.reduced$freq <- as.numeric(as.character(xy.reduced$freq))
  n <- sum(xy.reduced$freq)
  if (max(xvec)==0 & max(yvec)==0) {return(c(rep(1e-10,5),1))}
  #print(xy.reduced)
  
  # initial guess
  if (is.null(initial)) {
    xbar <- mean(xvec); ybar <- mean(yvec); xybar <- mean(c(xbar, ybar))
    s2.x <- var(xvec); s2.y <- var(yvec); if(is.na(s2.x)) {s2.x <- s2.y <- 1}
    s2.mean <- mean(s2.x, s2.y)
    cor.xy <- cor(xvec,yvec); if (is.na(cor.xy)) {cor.xy <- 0}
    zero <- sum(xvec == 0 & yvec == 0) / n
    initial <- rep(NA,5)
    initial[4] <- s2.mean /xybar
    initial[5] <- zero
    initial[2:3] <- c(xbar,ybar)/initial[4]
    initial[1] <- min(initial[2],initial[3]) * abs(cor.xy)
    initial[2:3] <-  initial[2:3] - initial[1]
    initial <- pmax(initial, 1e-5)
  }
  
  iter = 0
  param = initial
  if (showFlag) {print(c("iter", "a0", "a1", "a2", "b", "pi"))}
  repeat {
    iter = iter + 1
    #print(c(param))
    #print(lik(vec, pp=param[1], m0=param[2], m1=param[3], m2=param[4])) # debug
    param.old <- param # saving old parameters
    
    # updating
    expt <- dBvZINB2.Expt.vec(xy.reduced$x, xy.reduced$y, 
                              a0 = param[1], a1 = param[2], a2 = param[3], b = param[4], pp = param[5])
    expt <- as.vector(expt %*% xy.reduced$freq / n)
    # loglik = expt[1] * n
    # delta <- expt[9] / (expt[2] + expt[4])                   # delta = E(V) / (E(xi0 + xi2))
    # pi = E(Z)
    
    fn <- function(par.ab) {
      sum((par.ab[1:3]) * (expt[5:7] - log(par.ab[4]))) - sum(log(gamma(par.ab[1:3]))) - sum(expt[2:4]/par.ab[4])
    } # note argmax(pi) = EZ
    gr <- function(par.ab) {
      la <- expt[5:7] - log(par.ab[4]) - digamma(par.ab[1:3])
      lb <- sum(- par.ab[1:3] / par.ab[4] + expt[2:4] / par.ab[4]^2)
      return(c(la,lb))
    }
    
    fn.log = function (par.ab.l) { fn (exp(par.ab.l))}
    gr.log = function (par.ab.l) { gr (exp(par.ab.l))}
    param[1:4] = exp(optim(par=log(param[1:4]), fn = fn.log, gr = gr.log, control=list(fnscale=-1), method="BFGS")$par)
    param[5] = expt[8]
    density = dBvZINB2.vec(x = xy.reduced$x, y = xy.reduced$y, param[1], param[2], param[3], param[4], param[5], log=TRUE) %*% xy.reduced$freq 
    if (showFlag) {print(c(iter, round(param,5), expt[1] * n))} #lik: lik of previous iteration
    
    #  result <- try(multiroot(opt.vec, start=param.l[1:4], rtol=1e-20)$root, silent=TRUE)
    #  if (class(result)=="try-error") {
    #    initial = rep(1,5)
    #    result <- multiroot(opt.vec, start=param[1:4], rtol=1e-20)$root
    #  }
    
    #print (expt) #####
    if (max(abs(param - param.old)) <= tol) {
      names(param) <- c("a0", "a1", "a2", "b", "pi")
      return(param)
      break
    }
  }
  #result <- data.frame(a0 = param[1], a1 = param[2], a2 = param[3], b1 = param[4], b2 = param[5], pi = param[6])
  #return(result)
}


(a5[a5$deviance>1000,])
ML.BvZINB2(extractor(3), extractor(54), showFlag=TRUE)
ML.BvZINB2(extractor(8), extractor(58), showFlag=TRUE) #very slow:
# 373.00000     0.01645     0.16632     0.07736    49.25460     0.44795 -1899.10708



# direct maximization does not work!  !!! 
ML.BvZINB2.direct <- function (xvec, yvec, initial = NULL, tol=1e-8, showFlag=FALSE) { # does not work!!
  xy.reduced <- as.data.frame(table(xvec,yvec))
  names(xy.reduced) <- c("x", "y","freq")
  xy.reduced <- xy.reduced[xy.reduced$freq != 0,]
  xy.reduced$x <- as.numeric(as.character(xy.reduced$x))
  xy.reduced$y <- as.numeric(as.character(xy.reduced$y))
  xy.reduced$freq <- as.numeric(as.character(xy.reduced$freq))
  n <- sum(xy.reduced$freq)
  
  # initial guess
  if (is.null(initial)) {
    xbar <- mean(xvec); ybar <- mean(yvec); xybar <- mean(c(xbar, ybar))
    s2.x <- var(xvec); s2.y <- var(yvec); if(is.na(s2.x)) {s2.x <- s2.y <- 1}
    s2.mean <- mean(s2.x, s2.y)
    cor.xy <- cor(xvec,yvec); if (is.na(cor.xy)) {cor.xy <- 0}
    zero <- sum(xvec == 0 & yvec == 0) / n
    initial <- rep(NA,5)
    initial[4] <- s2.mean /xybar
    initial[5] <- zero
    initial[2:3] <- c(xbar,ybar)/initial[4]
    initial[1] <- min(initial[2],initial[3]) * abs(cor.xy)
    initial[2:3] <-  initial[2:3] - initial[1]
    initial <- pmax(initial, 1e-5)
  }
  lik.log = function (param.l) { lik.BvZINB2 (x = xvec, y = yvec, exp(param.l)) %>% print}
  param = exp(optim(par=log(initial), fn = lik.log, control=list(fnscale=-1))$par)
  return(param)
}
# direct + EM(for pi) : -> no use of EM for pi (slow. closed form solution exists!)
ML.BvZINB2.hybrid <- function (xvec, yvec, initial = NULL, tol=1e-8, showFlag=FALSE) { # does not work!!
  xy.reduced <- as.data.frame(table(xvec,yvec))
  names(xy.reduced) <- c("x", "y","freq")
  xy.reduced <- xy.reduced[xy.reduced$freq != 0,]
  xy.reduced$x <- as.numeric(as.character(xy.reduced$x))
  xy.reduced$y <- as.numeric(as.character(xy.reduced$y))
  xy.reduced$freq <- as.numeric(as.character(xy.reduced$freq))
  n <- sum(xy.reduced$freq)
  
  # initial guess
  if (is.null(initial)) {
    xbar <- mean(xvec); ybar <- mean(yvec); xybar <- mean(c(xbar, ybar))
    s2.x <- var(xvec); s2.y <- var(yvec); if(is.na(s2.x)) {s2.x <- s2.y <- 1}
    s2.mean <- mean(s2.x, s2.y)
    cor.xy <- cor(xvec,yvec); if (is.na(cor.xy)) {cor.xy <- 0}
    zero <- sum(xvec == 0 & yvec == 0) / n
    initial <- rep(NA,5)
    initial[4] <- s2.mean /xybar
    initial[5] <- zero
    initial[2:3] <- c(xbar,ybar)/initial[4]
    initial[1] <- min(initial[2],initial[3]) * abs(cor.xy)
    initial[2:3] <-  initial[2:3] - initial[1]
    initial <- pmax(initial, 1e-5)
  }
  
  iter = 0
  param = initial
  if (showFlag) {print(c("iter", "a0", "a1", "a2", "b", "pi"))}
  repeat {
    iter = iter + 1
    #print(c(param))
    #print(lik(vec, pp=param[1], m0=param[2], m1=param[3], m2=param[4])) # debug
    param.old <- param # saving old parameters
    
    # updating
    z.E <- dBvZINB2.Expt.Z.vec(xy.reduced$x, xy.reduced$y, 
                               a0 = param[1], a1 = param[2], a2 = param[3], b = param[4], pp = param[5])
    param[5] <- (z.E %*% xy.reduced$freq)/n
    lik.log = function (param.l) { lik.BvZINB2 (x = xvec, y = yvec, c(exp(param.l), pp = param[5]))}
    param = c(exp(optim(par=log(param[1:4]), fn = lik.log, control=list(fnscale=-1))$par),param[5])
    density = dBvZINB2.vec(x = xy.reduced$x, y = xy.reduced$y, param[1], param[2], param[3], param[4], param[5], log=TRUE) %*% xy.reduced$freq 
    if (showFlag) {print(c(iter, round(param,5), density))}
    if (max(abs(param - param.old)) <= tol) {
      names(param) <- c("a0", "a1", "a2", "b", "pi")
      return(param)
      break
    }
  }
  
  
}
# step-wise direct maximazation: pi -> alpha,beta -> pi -> ...
ML.BvZINB2.hybrid2 <- function (xvec, yvec, initial = NULL, tol=1e-8, showFlag=FALSE) { # does not work!!
  xy.reduced <- as.data.frame(table(xvec,yvec))
  names(xy.reduced) <- c("x", "y","freq")
  xy.reduced <- xy.reduced[xy.reduced$freq != 0,]
  xy.reduced$x <- as.numeric(as.character(xy.reduced$x))
  xy.reduced$y <- as.numeric(as.character(xy.reduced$y))
  xy.reduced$freq <- as.numeric(as.character(xy.reduced$freq))
  n <- sum(xy.reduced$freq)
  
  # initial guess
  if (is.null(initial)) {
    xbar <- mean(xvec); ybar <- mean(yvec); xybar <- mean(c(xbar, ybar))
    s2.x <- var(xvec); s2.y <- var(yvec); if(is.na(s2.x)) {s2.x <- s2.y <- 1}
    s2.mean <- mean(s2.x, s2.y)
    cor.xy <- cor(xvec,yvec); if (is.na(cor.xy)) {cor.xy <- 0}
    zero <- sum(xvec == 0 & yvec == 0) / n
    initial <- rep(NA,5)
    initial[4] <- s2.mean /xybar
    initial[5] <- zero
    initial[2:3] <- c(xbar,ybar)/initial[4]
    initial[1] <- min(initial[2],initial[3]) * abs(cor.xy)
    initial[2:3] <-  initial[2:3] - initial[1]
    initial <- pmax(initial, 1e-5)
  }
  
  iter = 0
  param = initial
  if (showFlag) {print(c("iter", "a0", "a1", "a2", "b", "pi"))}
  repeat {
    iter = iter + 1
    #print(c(param))
    #print(lik(vec, pp=param[1], m0=param[2], m1=param[3], m2=param[4])) # debug
    param.old <- param # saving old parameters
    
    lik.log = function (param.l) { lik.BvZINB2 (x = xvec, y = yvec, c(exp(param.l), pp = param[5]))}
    param = c(exp(optim(par=log(param[1:4]), fn = lik.log, control=list(fnscale=-1))$par),param[5])
    #lik.log2 = function (param.l2) { lik.BvZINB2 (x = xy.reduced$x, y = xy.reduced$y, c(pp = param[1:4], plogis(param.l2)))}
    #param = c(param[1:4],plogis(optim(par=qlogis(param[5]), fn = lik.log2, control=list(fnscale=-1))$par))
    lik.2 = function (param.2) { lik.BvZINB2 (x = xy.reduced$x, y = xy.reduced$y, c(pp = param[1:4], param.2))}
    param = c(param[1:4], optim(par=param[5], fn = lik.2, control=list(fnscale=-1), lower = 0, upper =1, method="Brent")$par)
    density = dBvZINB2.vec(x = xy.reduced$x, y = xy.reduced$y, param[1], param[2], param[3], param[4], param[5], log=TRUE) %*% xy.reduced$freq 
    if (showFlag) {print(c(iter, round(param,5), density))}
    if (max(abs(param - param.old)) <= tol) {
      names(param) <- c("a0", "a1", "a2", "b", "pi")
      return(param)
      break
    }
  }
  
  
}
# step-wise direct maximazation with closed form pi: pi -> alpha,beta -> pi -> ... !!!
ML.BvZINB2.hybrid3 <- function (xvec, yvec, initial = NULL, tol=1e-8, showFlag=FALSE) { # does not work!!
  xy.reduced <- as.data.frame(table(xvec,yvec))
  names(xy.reduced) <- c("x", "y","freq")
  xy.reduced <- xy.reduced[xy.reduced$freq != 0,]
  xy.reduced$x <- as.numeric(as.character(xy.reduced$x))
  xy.reduced$y <- as.numeric(as.character(xy.reduced$y))
  xy.reduced$freq <- as.numeric(as.character(xy.reduced$freq))
  n <- sum(xy.reduced$freq)
  n0 <- sum(xvec==0&yvec==0)
  
  # initial guess
  if (is.null(initial)) {
    xbar <- mean(xvec); ybar <- mean(yvec); xybar <- mean(c(xbar, ybar))
    s2.x <- var(xvec); s2.y <- var(yvec); if(is.na(s2.x)) {s2.x <- s2.y <- 1}
    s2.mean <- mean(s2.x, s2.y)
    cor.xy <- cor(xvec,yvec); if (is.na(cor.xy)) {cor.xy <- 0}
    zero <- sum(xvec == 0 & yvec == 0) / n
    initial <- rep(NA,5)
    initial[4] <- s2.mean /xybar
    # initial[5] <- zero
  initial[5] <- 0
    initial[2:3] <- c(xbar,ybar)/initial[4]
    initial[1] <- min(initial[2],initial[3]) * abs(cor.xy)
    initial[2:3] <-  initial[2:3] - initial[1]
    initial <- pmax(initial, 1e-5)
  }
  
  iter = 0
  param = initial
  if (showFlag) {print(c("iter", "a0", "a1", "a2", "b", "pi"))}
  repeat {
    iter = iter + 1
    #print(c(param))
    #print(lik(vec, pp=param[1], m0=param[2], m1=param[3], m2=param[4])) # debug
    param.old <- param # saving old parameters
    
    # 1. update alphas and betas
    lik.log = function (param.l) { lik.BvZINB2 (x = xvec, y = yvec, c(exp(param.l), pp = param[5]))}
    param = c(exp(optim(par=log(param[1:4]), fn = lik.log, control=list(fnscale=-1))$par),param[5])
    density = dBvZINB2.vec(x = xy.reduced$x, y = xy.reduced$y, param[1], param[2], param[3], param[4], param[5], log=TRUE) %*% xy.reduced$freq 
    if (showFlag) {print(c(iter, round(param,5), density))}
    
    # 2. update pi
    f00 = (param[4] + 1)^-sum(param[2:3]) * (2* param[4] + 1)^-param[1]
    pp = 1 - (n-n0)/n/(1-f00)
    pp = max(0,pp)
    param = c(param[1:4], pp)
    density = dBvZINB2.vec(x = xy.reduced$x, y = xy.reduced$y, param[1], param[2], param[3], param[4], param[5], log=TRUE) %*% xy.reduced$freq 
    if (showFlag) {print(c(iter, round(param,5), density))}
    if (max(abs(param - param.old)) <= tol) {
      names(param) <- c("a0", "a1", "a2", "b", "pi")
      return(param)
      break
    }
  }
  
  
}

# direct2 with gradient
dBvZINB2.gr <- function(x, y, a0, a1, a2, b, pp) {
  p = (2*b + 1) /(b + 1)
  gr.b.1 <- (x+y)/b - (x + y + a0) /(b + .5) - (a1 + a2) / (b + 1)
  gr.b <- function(k, m) {(k + m) * ((b + .5)^-1 - (b + 1)^-1) + gr.b.1}
  gr.a0.1 <- - digamma(a0) - log(1+2*b)
  gr.a0 <- function(k, m) {digamma(x +y +a0 -k -m) + gr.a0.1}
  {if (x==0 & y==0) {- log(1 + 2*b)} else  {digamma(x +y +a0 -k -m) + gr.a0.1}}
  gr.a1.1 <- - digamma(a1) - log(1 +b)
  gr.a1 <- function(k) {if (k==0) {- log(1 +b)} else  {digamma(a1 + k) + gr.a1.1}}
  gr.a2.1 <- - digamma(a2) - log(1 +b)
  gr.a2 <- function(m) {if (m==0) {- log(1 +b)} else  {digamma(a2 + m) + gr.a2.1}}
  
  l1 <- function(k, m) exp(lgamma(a1 + k) - lgamma(k+1) - lgamma(a1) + lgamma(x + y + a0 -m -k) - lgamma(x -k +1) - lgamma(a0 + y - m) + k *log(p)
                           + lgamma(m +a2) - lgamma(m+1) - lgamma(a2) + lgamma(y +a0 -m) - lgamma(y -m +1) - lgamma(a0) + m *log(p))
  l2 <- - (+x+y+a0)*log(1 + 2*b) + (x+y) * log(b) -(a1 + a2)*log(1+b)
  l2 <- exp(l2)
  l1.mat <- sapply(0:x, function(k) sapply(0:y, l1, k=k))
  l1.mat <- (l1.mat * l2) #%>%print
  l1.mat.zi <- l1.mat * (1-pp) + pp
  
  l1.sum <- sum(l1.mat)
  
  gr.b.mat <- sapply(0:x, function(k) sapply(0:y, gr.b, k=k))
  gr.b.mat <- l1.mat * gr.b.mat
  gr.b.sum <- sum(gr.b.mat)/sum(l1.mat.zi)
  
  gr.a0.mat <- sapply(0:x, function(k) sapply(0:y, gr.a0, k=k))
  gr.a0.mat <- l1.mat * gr.a0.mat
  gr.a0.sum <- sum(gr.a0.mat)/sum(l1.mat.zi)
  
  gr.a1.mat <- matrix(sapply(0:x, gr.a1),x+1, y+1) #%>% print
  gr.a1.mat <- l1.mat * t(gr.a1.mat)
  gr.a1.sum <- sum(gr.a1.mat)/sum(l1.mat.zi)
  
  gr.a2.mat <- matrix(sapply(0:y, gr.a2), y+1, x+1) #%>% print
  gr.a2.mat <- l1.mat * gr.a2.mat
  gr.a2.sum <- sum(gr.a2.mat)/sum(l1.mat.zi)
  
  result <- list(logdensity=log(l1.sum), gradient = c(gr.a0.sum, gr.a1.sum, gr.a2.sum, gr.b.sum))
  return(result)
} # DONE
dBvZINB2.gr.vec <- Vectorize(dBvZINB2.gr)
ML.BvZINB2.direct2 <- function (xvec, yvec, initial = NULL, tol=1e-8, showFlag=FALSE) { # does not work!!
  xy.reduced <- as.data.frame(table(xvec,yvec))
  names(xy.reduced) <- c("x", "y","freq")
  xy.reduced <- xy.reduced[xy.reduced$freq != 0,]
  xy.reduced$x <- as.numeric(as.character(xy.reduced$x))
  xy.reduced$y <- as.numeric(as.character(xy.reduced$y))
  xy.reduced$freq <- as.numeric(as.character(xy.reduced$freq))
  n <- sum(xy.reduced$freq)
  
  # initial guess
  if (is.null(initial)) {
    xbar <- mean(xvec); ybar <- mean(yvec); xybar <- mean(c(xbar, ybar))
    s2.x <- var(xvec); s2.y <- var(yvec); if(is.na(s2.x)) {s2.x <- s2.y <- 1}
    s2.mean <- mean(s2.x, s2.y)
    cor.xy <- cor(xvec,yvec); if (is.na(cor.xy)) {cor.xy <- 0}
    zero <- sum(xvec == 0 & yvec == 0) / n
    initial <- rep(NA,5)
    initial[4] <- s2.mean /xybar
    initial[5] <- zero
    initial[2:3] <- c(xbar,ybar)/initial[4]
    initial[1] <- min(initial[2],initial[3]) * abs(cor.xy)
    initial[2:3] <-  initial[2:3] - initial[1]
    initial <- pmax(initial, 1e-5)
  }
  
  gr.1 = function (param) {
    print(param)
    val <- dBvZINB2.gr.vec( x = xy.reduced$x, y = xy.reduced$y,  a0 = param[1], a1 = param[2], a2 = param[3], b = param[4], pp = param[5])
    gr <- sapply(val[2,],cbind) %*% xy.reduced$freq
    # print(c(param,gr)) ###
    return(gr)
  }
  
  #log-scaled params: param.l
  fn.log = function (param.l) { lik.BvZINB2 (x=xvec, y=yvec, c(exp(param.l[1:4]),plogis(param.l[5]))) %>% print}
  gr.log = function (param.l) { 
      param = c(exp(param.l[1:4]), plogis(param.l[5]))
      addendum = param *c(1,1,1,1,1-param[5])
      if (showFlag) {print(param)} #####
    return(as.vector(gr.1 (param)) * addendum) 
  }
  
  initial.l = c(log(initial[1:4]),qlogis(initial[5]))
  param.l = optim(par=initial.l, fn = fn.log, gr = gr.log, control=list(fnscale=-1), method="BFGS")$par
  param = c(exp(param.l[1:4]), plogis(param.l[5]))
  # print(c(param, lik.BvZINB2 (x=xvec, y=yvec, param)))
  return(param)
}
# direct3 with gradient for alpha and beta, and closed form for pi with repetition steps
ML.BvZINB2.direct3 <- function (xvec, yvec, initial = NULL, tol=1e-8, showFlag=FALSE) { # does not work!!
  xy.reduced <- as.data.frame(table(xvec,yvec))
  names(xy.reduced) <- c("x", "y","freq")
  xy.reduced <- xy.reduced[xy.reduced$freq != 0,]
  xy.reduced$x <- as.numeric(as.character(xy.reduced$x))
  xy.reduced$y <- as.numeric(as.character(xy.reduced$y))
  xy.reduced$freq <- as.numeric(as.character(xy.reduced$freq))
  n <- sum(xy.reduced$freq)
  n0 <- sum(xvec==0&yvec==0)
  
  # initial guess
  if (is.null(initial)) {
    xbar <- mean(xvec); ybar <- mean(yvec); xybar <- mean(c(xbar, ybar))
    s2.x <- var(xvec); s2.y <- var(yvec); if(is.na(s2.x)) {s2.x <- s2.y <- 1}
    s2.mean <- mean(s2.x, s2.y)
    cor.xy <- cor(xvec,yvec); if (is.na(cor.xy)) {cor.xy <- 0}
    zero <- sum(xvec == 0 & yvec == 0) / n
    initial <- rep(NA,5)
    initial[4] <- s2.mean /xybar
    initial[5] <- 0 #zero
    initial[2:3] <- c(xbar,ybar)/initial[4]
    initial[1] <- min(initial[2],initial[3]) * abs(cor.xy)
    initial[2:3] <-  initial[2:3] - initial[1]
  initial[5] <- 0 ##
    initial <- pmax(initial, 1e-5)
  }
  
  gr.1 = function (par.ab, pp) {
    val <- dBvZINB2.gr.vec( x = xy.reduced$x, y = xy.reduced$y,  a0 = par.ab[1], a1 = par.ab[2], a2 = par.ab[3], b = par.ab[4], pp = pp)
    gr <- sapply(val[2,],cbind) %*% xy.reduced$freq
    # print(c(param,gr)) ###
    return(gr[1:4]) # don't need gr for pi
  }
  
  #log-scaled params: param.l
  fn.log = function (par.ab, pp) { lik.BvZINB2 (x=xvec, y=yvec, c(exp(par.ab), pp))}
  gr.log = function (par.ab, pp) { 
    ## if (showFlag) {print(c(exp(par.ab), pp))} #####
    return(as.vector(gr.1 (exp(par.ab), pp)) * exp(par.ab)) 
  }
  
  iter = 0
  param = initial
  if (showFlag) {print(c("iter", "a0", "a1", "a2", "b", "pi"))}
  repeat {
    iter = iter + 1
    #print(c(param))
    #print(lik(vec, pp=param[1], m0=param[2], m1=param[3], m2=param[4])) # debug
    param.old <- param # saving old parameters
    
    # update ab
    par.ab.l = log(param[1:4])
    fn.log.ab <- function(s) fn.log(s, pp=param[5])
    gr.log.ab <- function(s) gr.log(s, pp=param[5])
    par.ab.l = optim(par=par.ab.l, fn = fn.log.ab, gr = gr.log.ab, control=list(fnscale=-1), method="BFGS")$par
    param = c(exp(par.ab.l), param[5])
    
    # update pi
    f00 = (param[4] + 1)^-sum(param[2:3]) * (2* param[4] + 1)^-param[1]
    pp = 1 - (n-n0)/n/(1-f00)
    pp = max(0,pp)
    param = c(param[1:4], pp)
    density = dBvZINB2.vec(x = xy.reduced$x, y = xy.reduced$y, param[1], param[2], param[3], param[4], param[5], log=TRUE) %*% xy.reduced$freq 
    
    if (showFlag) {print(c(iter, round(param,5), density))}
    if (max(abs(param - param.old)) <= tol) {
      names(param) <- c("a0", "a1", "a2", "b", "pi")
      return(param)
      break
    }
  }
  # print(c(param, lik.BvZINB2 (x=xvec, y=yvec, param)))
  return(param)
}


# simple tests
if (FALSE) {
  ML.BvZINB2(c(10,1,1),c(10,1,2), showFlag=TRUE) # c(1.274400, 0.002953337, 0.250098153, 2.951314982, 2.840646551, 0)
  ML.BvNB2(c(10,1,1),c(10,1,2), showFlag=TRUE)   # c(1.189742, 0.0001505849, 0.000292592, 3.361761, 3.641081,0)
  apply(dBvZINB3.Expt.vec(c(10,1,1),c(10,1,2), 1.189742, 0.0001505849, 0.000292592, 3.361761, 3.641081,0),1, mean)
  1.18971 0.0001450095 0.0004408919 3.501135
  (1.18971/(1+2*3.501135) + ( 0.0001450095 + 0.0004408919)/(1+3.501135)) / (1.18971+ 0.0001450095+ 0.0004408919)
  ML.BvZINB2(c(10,1,1),c(10,1,2), initial= c(1.18971, 0.0001450095, 0.0004408919, 3.501135, 0), showFlag=TRUE)
  
  ML.BvZINB2.direct(c(10,1,1),c(10,1,2))  
  # 1.188961e+00 1.911217e-08 1.354178e-14 3.504383e+00 3.364621e-07
  # -12.924
  lik.BvZINB2(c(10,1,1),c(10,1,2), c(1.188961e+00, 1.911217e-08, 1.354178e-14, 3.504383e+00, 3.364621e-07))
  
  lik.BvZINB2(c(10,1,1),c(10,1,2), c(1.18971, 0.0001450095, 0.0004408919, 3.501135, 0))
  lik.BvZINB2(c(10,1,1),c(10,1,2), c(1.1898696803, 0.0002000401, 0.0006001678, 3.5001799412, 0.0000000000 ))
    
  ML.BvZINB3(c(0,1,1),c(0,1,5), showFlag=TRUE)
  lik.BvZINB3(c(0,1,1),c(0,1,2), c(3.003010e-01,  2.348166e-03,  1.514485e-01,  1.895357e+00,  2.210670e+00, 4.940656e-324) )
  
  tt(1)
  ML.BvNB2(extractor(1), extractor(4), showFlag=TRUE)  #0.0004349035 0.009488825 0.003788559 68.25597 9.835188
  ML.BvZINB2(extractor(1), extractor(4), initial = c(0.0003058505, 0.01004748, 0.002164966, 57.80222,  0.36109), showFlag=TRUE)
  ML.BvNB2(extractor(1), extractor(3), showFlag=TRUE)
  ML.BvZINB2(extractor(1), extractor(3), showFlag=TRUE)
  
  ML.BvZINB2.direct(extractor(1), extractor(4))  
  lik.BvZINB2(extractor(1), extractor(4), c(2.215115e-05, 4.774164e-01, 5.562271e-02, 4.047793e+01, 1.349860e-04))
  ML.BvZINB2.direct(extractor(1), extractor(4), initial=c(0.00033 ,   0.01253,    0.00273   ,57.34690,    0.18943))
  lik.BvZINB2(extractor(1), extractor(4), c(4.613197e-12, 4.768649e-01, 5.551846e-02, 4.039589e+01, 3.821786e-06))
  
  lik.BvZINB3(extractor(1), extractor(4),c(0.0004349035, 0.009488825, 0.003788559, 68.25597, 9.835188, .95)) # -391.5657
  lik.BvNB3(extractor(1), extractor(4),c(0.0004349035, 0.009488825, 0.003788559, 68.25597, 9.835188)) # -308.7
  
  lik.BvNB2(extractor(1), extractor(3), c(1.033294e-05 ,0.008663334 ,0.05955079, 125.6262))
  [1] -1459.402
  
  
  
  
  {-3.356625e-01  5.371056e-02  1.736143e+00  6.363566e-01 -2.219790e+03 -9.442568e+01 -2.561689e+02
  [8]  8.643841e-01  1.053820e-01
  -3.356618e-01  7.034855e-02  2.092216e+00  7.615860e-01 -2.137599e+03 -9.131154e+01 -2.542944e+02
  [8]  8.284366e-01  1.279334e-01
  apply(dBvZINB3.Expt.vec(extractor(1), extractor(4), 0.0004349035, 0.009488825, 0.003788559, 68.25597, 9.835188, .95), 1, mean)
  igamma(-log(68.25597) + 3.227275e-02)
  digamma(0.0004349035)+ log(68.25597)}
  
  tt(2) # 31secs
  
  ML.BvNB3(extractor(2), extractor(5), method="BFGS", abstol=1e-30)
  ML.BvNB3(extractor(2), extractor(5), method="Nelder-Mead", abstol=1e-30)
  sum(sapply(dBvNB3.gr.vec(extractor(2), extractor(5), 0.0004670043, 0.002597285, 0.01597736, 10.48014, 35.76876)[1,],cbind))
  sum(sapply(dBvNB3.gr.vec(extractor(2), extractor(5), 0.0004686926, 0.002604057, 0.01598007, 10.45238, 35.75869)[1,],cbind))
  
  tt(1)
  ML.BvNB3(extractor(1), extractor(3), method="BFGS")
  tt(2) #8sec
  tt(1)
  ML.BvNB3(extractor(1), extractor(38), method="BFGS", showFlag=TRUE)
  ML.BvNB3(extractor(1), extractor(38), method="Nelder-Mead", showFlag=TRUE)
  tt(2) #31sec
  #lik.BvNB3(extractor(1), extractor(38), c(5.790158e-03, 4.300688e-03, 7.836757e-02, 7.586956e+01, 1.015767e+02))
  #lik.BvNB3(extractor(1), extractor(38), c())
  
}
# some deviance tests
if (FALSE) {
  tt(1)
  a <- ML.BvNB3(extractor(11), extractor(16))
  tt(2) #1.8 sec!
  lik.BvNB3(extractor(11), extractor(16),a) #lik = -805
  
  tt(1)
  a <- ML.BvNB3(extractor(11), extractor(16))
  tt(2) #1.9 sec!
  lik.BvNB3(extractor(11), extractor(16),a) #lik = -805
  
  tt(1)
  a <- ML.BvNB3(extractor(30), extractor(50))
  tt(2) #0.2 sec!
  lik.BvNB3(extractor(30), extractor(50),a) #lik = -113
  
  tt(1)
  a <- ML.BvNB3(extractor(8), extractor(11))
  tt(2) #17 sec!
  lik.BvNB3(extractor(8), extractor(11),a) #lik = -1938
  
  tt(1)
  a <- ML.BvNB3(extractor(1), extractor(3))
  tt(2) #17 sec!
  lik.BvNB3(extractor(1), extractor(3),a) #lik = -1458
  
}

tt(1)
a <- pairwise.MLE(data=data[iset.Mm.c2[[1]],], ML.fun = ML.BvNB3, showFlag=TRUE)
tt(2)

tt(1); ML.BvZINB2.direct(extractor(2), extractor(12), showFlag=TRUE); tt(2) #31sec lik: -49.65
tt(1); ML.BvZINB2.direct2(extractor(2), extractor(12), showFlag=TRUE); tt(2) #2.3sec lik: -49.68
tt(1); ML.BvZINB2.direct3(extractor(2), extractor(12), showFlag=TRUE); tt(2) #4.3sec lik: -49.68
tt(1); ML.BvZINB2.hybrid(extractor(2), extractor(12), showFlag=TRUE); tt(2)  #>1min lik: > -49.34
tt(1); ML.BvZINB2.hybrid2(extractor(2), extractor(12), showFlag=TRUE); tt(2) # 1.1min lik: -49.68
tt(1); ML.BvZINB2.hybrid3(extractor(2), extractor(12), showFlag=TRUE); tt(2) # 1.6min lik: -49.35

# initial guess of pi = n0/n
tt(1); ML.BvZINB2.direct(extractor(1), extractor(8), showFlag=TRUE); tt(2) # 1.3min lik: -1520
tt(1); ML.BvZINB2.direct2(extractor(1), extractor(8), showFlag=TRUE); tt(2) # 9sec lik: -1609
tt(1); ML.BvZINB2.direct3(extractor(1), extractor(8), showFlag=TRUE); tt(2) # 17sec lik: -1582
tt(1); ML.BvZINB2.hybrid(extractor(1), extractor(8), showFlag=TRUE); tt(2)  #>1min lik: > -
tt(1); ML.BvZINB2.hybrid2(extractor(1), extractor(8), showFlag=TRUE); tt(2) # 3.5min lik: -1523
tt(1); ML.BvZINB2.hybrid3(extractor(1), extractor(8), showFlag=TRUE); tt(2) # >20min lik: -1520

tt(1); ML.BvZINB2.EM2(extractor(1), extractor(8), showFlag=TRUE); tt(2) # >20min lik: -1520


 tt(1)
MLE.Geneset1$BvNB3 <- pairwise.MLE(data[iset.Mm.c2[[1]],], ML.fun = ML.BvNB3, showFlag=TRUE)  ## 1.8hrs
tt(2)

ML.BvNB2(0,155)

# not developed yet
if (FALSE) {

### 3. Deviance
dev.BvNB3 <- function(xvec, yvec, param = NULL, a0 = NULL, a1 = NULL, a2= NULL, b1 = NULL, b2 = NULL) {
  # If params = NULL, apply ML.BvNB3. if else, apply those params
  if (is.null (param)) { 
    if (is.null (a0) | is.null (a1) | is.null (a2) | is.null (b1) | is.null (b2)) {
      param = ML.BvNB3 (xvec = xvec, yvec = yvec)
    }
    else { param = c(a0, a1, a2, b1, b2)}
  }
  # Log-likelihood of the BvNB2 model
  lik.model <- lik.BvNB3 (x = xvec, y = yvec, param = param)
  
  # Reduced calculation
  xy.reduced <- as.data.frame(table(xvec,yvec))
  names(xy.reduced) <- c("x", "y","freq")
  xy.reduced <- xy.reduced[xy.reduced$freq != 0,]
  xy.reduced$x <- as.numeric(as.character(xy.reduced$x))
  xy.reduced$y <- as.numeric(as.character(xy.reduced$y))
  xy.reduced$freq <- as.numeric(as.character(xy.reduced$freq))
  
  # Saturated model BZIP params
  ML.BvNB3.vec <- Vectorize(ML.BvNB3)
  param.sat <- t(ML.BvNB3.vec(xy.reduced$x, xy.reduced$y)) # %>%print
  param.sat <- do.call(rbind, param.sat)  #"matrix of lists" into "a matrix"
  # tmp.a <<- param.sat ####
  lik.sat   <- sum(dBvNB3.vec(x= xy.reduced$x, y = xy.reduced$y, a0 = param.sat[,1], a1 = param.sat[,2], a2 = param.sat[,3], b1 = param.sat[,4], b2 = param.sat[,5],  log = TRUE) * xy.reduced$freq)
  
  return(data.frame(model.likelihood = lik.model, satrtd.likelihood = lik.sat, deviance = 2*(lik.sat - lik.model)))
}
dev.BvNB3(extractor(1), extractor(2), param=as.vector(a[,2]))  #455  vs BZIP.B:3150
dev.BvNB3(extractor(1), extractor(5), param=as.vector(a[,5]))  #849  vs BZIP.B:2877
dev.BvNB3(extractor(1), extractor(3))  #2028 vs BvNB2: 2030
dev.BvNB3(extractor(1), extractor(8))  #2066 vs BvNB2: 2077

## heatmap
a <- ML.BvNB3(extractor(1), extractor(3), method="BFGS")
tmp <- sapply(0:50, function(r) sapply(0:50, function(s) dBvNB3(s,r, a[1,1],a[1,2],a[1,3],a[1,4],a[1,5])))
plot_ly(z = tmp[1:2,2:10], type = "heatmap")
table(extractor(1), extractor(3))

MLE.Geneset1$BvNB3 ########
a <- cbind(1,2,dev.BvNB2(extractor(1), extractor(2), param = c(1,1,1,1,1,1,1)))
tt(1)
for (i in 1:dim(MLE.Geneset1$BvNB2)[1]) {
  #if (i ==300 | i == 349 ) {a[i,] <- c(a1, a2, NA,NA,NA)}
  #if (i >= 350) {
  tmp <- MLE.Geneset1$BvNB2[i,]
  a1 <- tmp[1]; a2 <- tmp[2]
  a3 <- tmp[4:7] #params
  a[i,] <- c(a1,a2,dev.BvNB2(extractor(as.numeric(a1),1), extractor(as.numeric(a2),1), param = a3))
  print(a[i,])  
}
tt(2) # 5.8 hours

MLE.Geneset1$BvNB2$dev <- a$deviance

sum(MLE.Geneset1$BvNB2$dev>863); sum(MLE.Geneset1$BvNB2$dev<863)
plot3.2.1 <-  ggplot(MLE.Geneset1$BvNB2, aes(non0.min, dev)) +
  geom_point(aes(color=1-non0.min)) + ggtitle(paste0("deviance of BvNB2 model - Geneset 1")) +
  ylim(c(0,25000)) +
  xlab("nonzero-count proprotion(min)") + ylab("Deviance of BvNB2") + geom_abline(intercept=863, slope=0, linetype, color="red", size=0.1)
png2(plot3.2.1, width = 720, height = 360)
}