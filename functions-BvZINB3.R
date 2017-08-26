### This is a zero-inflated version of BvNB3
# Needs dBvNB3 function
library(distr)

### 1. Density, likelihood
dBvZINB3 <- function(x, y, a0, a1, a2, b1, b2, pp, log=FALSE) {
  result <- dBvNB3(x=x, y=y, a0=a0, a1=a1, a2=a2, b1=b1, b2=b2, log=FALSE) * (1-pp) + ifelse(x + y == 0, pp, 0)
  return(ifelse(log, log(result), result))
}
if (FALSE) {
  dBvNB3(1,1,1,1,1,1,2)
  dBvZINB3(1,1,1,1,1,1,2,.5) ;dBvZINB3(0,0,1,1,1,1,2,.5) 
  tmp <- sapply(0:50, function(r) sapply(0:50, function(s) dBvZINB3(s,r,1,1,1,1,.5,.5)))
  sum(tmp) #1
}
dBvZINB3.vec <- Vectorize(dBvZINB3)
lik.BvZINB3 <- function(x, y, param) {
  sum(log(dBvZINB3.vec(x, y, param[1], param[2], param[3], param[4], param[5], param[6])))
}


### 2.EM

### nonzero cells: (1-pp) was not multiplied by!!! this caused decreasing likelihood in EM
dBvZINB3.Expt <- function(x, y, a0, a1, a2, b1, b2, pp) {
  p1 = (b1 + b2 + 1) /(b1 + 1); p2 = (b1 + b2 + 1) /(b2 + 1)
  f00 = (b1 + 1)^-a1 * (b2 + 1)^-a2 * (b1 + b2 + 1)^-a0
  if (x == 0 & y == 0) {
    den <- f00 *(1-pp) + pp
    rate1 <- f00 *(1-pp) / den
    rate2 <- pp / den
    xi0.E <- (a0 * b1) * (rate1 /(b1 + b2 + 1) + rate2)
    xi1.E <- (a1 * b1) * (rate1 /(b1 + 1) + rate2)
    xi2.E <- (a2 * b1) * (rate1 /(b2 + 1) + rate2)
    z.E   <- rate2
    v.E   <- (a0 + a2) * b2 * rate2
    log.xi0.E <- digamma(a0) + log(b1) - log(1 + b1 + b2) * rate1
    log.xi1.E <- digamma(a1) + log(b1) - log(1 + b1) * rate1
    log.xi2.E <- digamma(a2) + log(b1) - log(1 + b2) * rate1
    result <- c(logdensity = log(den), xi0.E = xi0.E, xi1.E = xi1.E, xi2.E = xi2.E, 
                log.xi0.E = log.xi0.E, log.xi1.E = log.xi1.E, log.xi2.E = log.xi2.E, z.E = z.E, v.E = v.E)
    return(result)
  } else {
    xi0.E1 <- function(k, m) {x - k + y - m + a0}; xi0.E2 <- (1 + b1 + b2)/b1
    xi1.E1 <- function(k) {k + a1}; xi1.E2 <- (1 + b1)/b1
    xi2.E1 <- function(m) {m + a2}; xi2.E2 <- (1 + b2)/b1
    log.xi0.E1 <- function(k, m) {digamma(x - k + y - m + a0)}; # log.xi0.E2 <- log(xi0.E2)
    log.xi1.E1 <- function(k) {digamma(k + a1)}; # log.xi1.E2 <- log(xi1.E2)
    log.xi2.E1 <- function(m) {digamma(m + a2)}; # log.xi2.E2 <- log(xi2.E2)
    
    l1 <- function(k, m) exp(lgamma(a1 + k) - lgamma(k+1) - lgamma(a1) + lgamma(x + y + a0 -m -k) - lgamma(x -k +1) - lgamma(a0 + y - m) + k *log(p1)
                             + lgamma(m + a2) - lgamma(m+1) - lgamma(a2) + lgamma(y +a0 -m) - lgamma(y -m +1) - lgamma(a0) + m *log(p2))
    l2 <- - (+x+y+a0)*log(1 + b1 + b2) + x * log(b1) + y * log(b2) - a1 * log(1 + b1) - a2 * log(1 + b2)
    l2 <- exp(l2) *(1-pp)
    
    l1.mat <- sapply(0:x, function(k) sapply(0:y, l1, k=k))
    l1.mat <- (l1.mat * l2) #%>%print
    l1.sum <- sum(l1.mat) 
    
    xi0.mat <- sapply(0:x, function(k) sapply(0:y, xi0.E1, k=k))
    xi0.mat <- l1.mat * xi0.mat
    xi0.E <- sum(xi0.mat)/sum(l1.mat)/xi0.E2
    # if (xy0) {xi0.E <- (xi0.E * f00 * (1-pp) + a1 * b1 * pp) / (f00 *(1 - pp) + pp)}
    
    xi1.mat <- matrix(sapply(0:x, xi1.E1),x+1, y+1) #%>% print
    xi1.mat <- l1.mat * t(xi1.mat)
    xi1.E <- sum(xi1.mat)/sum(l1.mat)/xi1.E2
    
    xi2.mat <- matrix(sapply(0:y, xi2.E1), y+1, x+1) #%>% print
    xi2.mat <- l1.mat * xi2.mat
    xi2.E <- sum(xi2.mat)/sum(l1.mat)/xi2.E2
    
    log.xi0.mat <- sapply(0:x, function(k) sapply(0:y, log.xi0.E1, k=k))
    log.xi0.mat <- l1.mat * log.xi0.mat
    log.xi0.E <- sum(log.xi0.mat)/sum(l1.mat) - log(xi0.E2)
    
    log.xi1.mat <- matrix(sapply(0:x, log.xi1.E1),x+1, y+1) #%>% print
    log.xi1.mat <- l1.mat * t(log.xi1.mat)
    log.xi1.E <- sum(log.xi1.mat)/sum(l1.mat) - log(xi1.E2)
    
    log.xi2.mat <- matrix(sapply(0:y, log.xi2.E1), y+1, x+1) #%>% print
    log.xi2.mat <- l1.mat * log.xi2.mat
    log.xi2.E <- sum(log.xi2.mat)/sum(l1.mat) - log(xi2.E2)
    
    z.E <- 0
    v.E <- y
    
    result <- c(logdensity = log(l1.sum), xi0.E = xi0.E, xi1.E = xi1.E, xi2.E = xi2.E, 
                log.xi0.E = log.xi0.E, log.xi1.E = log.xi1.E, log.xi2.E = log.xi2.E, z.E = z.E, v.E = v.E)
    return(result)    
  }
}
dBvZINB3.Expt.vec <- Vectorize(dBvZINB3.Expt)

if (FALSE) {
  tmp <- dBvZINB3.Expt.vec(c(1,1,1),c(0,1,2),1,1,1,1,2,.5)
  tmp <- dBvZINB3.Expt.vec(c(0,1,1),c(0,1,2),1,1,1,1,2,.5)
  tmp <- dBvZINB3.Expt.vec(extractor(1),extractor(2),1,1,1,1,2,.5)
  t(tmp)[21:40,]
  dBvZINB3.Expt.vec(c(10,1,2),c(10,1,1), 1.193013282, 0.003336139, 0.002745513, 3.618842924, 3.341625901, 0.000000000)
  dBvZINB3.Expt.vec(c(10,1,1),c(10,1,2), 1.193013282, 0.002745513, 0.003336139, 3.341625901, 3.618842924, 0.000000000)
  # the above two should have the symmetric results!
  dBvZINB3.Expt(10, 1, 1,2,1,2,1,0)
  dBvZINB3.Expt(1, 10, 1,1,2,1,2,0)
}

ML.BvZINB3 <- function (xvec, yvec, initial = NULL, tol=1e-8, showFlag=FALSE) {
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
    cor.xy <- cor(xvec,yvec); if (is.na(cor.xy)) {cor.xy <- 0}
    zero <- sum(xvec == 0 & yvec == 0) / n
    initial <- rep(NA,5)
    initial[4] <- s2.x /xbar
    initial[5] <- s2.y /ybar
    initial[6] <- zero
    initial[2:3] <- c(xbar,ybar)/initial[4:5]
    initial[1] <- min(initial[2:3]) * abs(cor.xy)
    initial[2:3] <-  initial[2:3] - initial[1]
    initial <- pmax(initial, 1e-5)
  print(initial) ###
  }
  
  iter = 0
  param = initial
  if (showFlag) {print(c("iter", "a0", "a1", "a2", "b1", "b2", "pi"))}
  repeat {
    iter = iter + 1
    #print(c(param))
    #print(lik(vec, pp=param[1], m0=param[2], m1=param[3], m2=param[4])) # debug
    param.old <- param # saving old parameters
    
    # updating
    expt <- dBvZINB3.Expt.vec(xy.reduced$x, xy.reduced$y, 
                              a0 = param[1], a1 = param[2], a2 = param[3], b1 = param[4], b2 = param[5], pp = param[6])
    expt <- as.vector(expt %*% xy.reduced$freq / n)
    # loglik = expt[1] * n
    delta <- expt[9] / (expt[2] + expt[4])                   # delta = E(V) / (E(xi0 + xi2))
    param[6] = expt[8]                                                    # pi = E(Z)
    
    opt.vec <- function(par.ab) {
      par.ab <- exp(par.ab)
      r1 <- sum(expt[2:4]) - sum(par.ab[1:3]) * par.ab[4]
      r2 <- expt[5:7] - digamma(par.ab[1:3]) - log(par.ab[4])
      # print(c(r1,r2)) ###
      return(c(r1,r2))
    }
    param.l <- log(param)
    result <- try(multiroot(opt.vec, start=param.l[1:4])$root, silent=TRUE)
    if (class(result)=="try-error") {
      initial = rep(1,5)
      result <- multiroot(opt.vec, start=initial[1:4], rtol=1e-20)$root
    }
    param[1:4] <- exp(result)
    param[5]    <- param[4] * delta                                       # b2
    param[6] = expt[8]
    
  #print (expt) #####
    if (showFlag) {print(c(iter, round(param,5), expt[1] * n))} #lik: lik of previous iteration
    if (max(abs(param - param.old)) <= tol) {
      names(param) <- c("a0", "a1", "a2", "b1", "b2", "pi")
      return(param)
      break
    }
  }
  #result <- data.frame(a0 = param[1], a1 = param[2], a2 = param[3], b1 = param[4], b2 = param[5], pi = param[6])
  #return(result)
}
# simple tests
if (FALSE) {
  ML.BvZINB3(c(10,1,1),c(10,1,2), showFlag=TRUE) # c(1.193014333 0.002745514 0.003336045 3.341621165 3.618839217 0.000000000 )
  ML.BvNB3(c(10,1,1),c(10,1,2), showFlag=TRUE)   # c(1.189742, 0.0001505849, 0.000292592, 3.361761, 3.641081,0)
  
  ML.BvZINB3(c(0,1,1),c(0,1,5), showFlag=TRUE)
  ML.BvNB3(c(0,1,1),c(0,1,5))
  tt(1)
  ML.BvNB3(extractor(1), extractor(4), showFlag=TRUE)  #0.0004349035 0.009488825 0.003788559 68.25597 9.835188
  ML.BvZINB3(extractor(1), extractor(4), initial = c(0.00080,    0.01505 ,   0.00621   ,67.41456   , 9.17832  ,  0.36109), showFlag=TRUE)
  ML.BvZINB3(extractor(1), extractor(3), showFlag=TRUE)
  # 0.000799916  0.015057420  0.006208375 67.414607790  9.180617081  0.361266622
  
  lik.BvZINB3(extractor(1), extractor(4),c(0.0004349035, 0.009488825, 0.003788559, 68.25597, 9.835188, .95)) # -391.5657
  lik.BvNB3(extractor(1), extractor(4),c(0.0004349035, 0.009488825, 0.003788559, 68.25597, 9.835188)) # -308.7
  
  
  3.404230e-05 9.740676e-03 5.435834e-03 7.059027e+01 6.627206e+00 9.500000e-01
  
  
  tt(2) # 31secs
  
  tt(1)
  ML.BvZINB3(extractor(1), extractor(3),showFlag=TRUE)
  ML.BvZINB3(extractor(1), extractor(3),initial=c(1.733055e-05, 0.009879464, 0.05864169, 69.22358, 134.6264,0),showFlag=TRUE)
  ML.BvNB3(extractor(1), extractor(3),showFlag=TRUE)
  lik.BvZINB3(extractor(1), extractor(3),c(1.733055e-05, 0.009879464, 0.05864169, 69.22358, 134.6264,0)) #1485.486
  tt(2) #8sec
  tt(1)
  ML.BvNB3(extractor(1), extractor(38), method="BFGS", showFlag=TRUE)
  ML.BvNB3(extractor(1), extractor(38), method="Nelder-Mead", showFlag=TRUE)
  tt(2) #31sec
  #lik.BvNB3(extractor(1), extractor(38), c(5.790158e-03, 4.300688e-03, 7.836757e-02, 7.586956e+01, 1.015767e+02))
  #lik.BvNB3(extractor(1), extractor(38), c())
  
}


# EM with booster
ML.BvZINB3.2 <- function (xvec, yvec, initial = NULL, tol=1e-8, showFlag=FALSE) {
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
    cor.xy <- cor(xvec,yvec); if (is.na(cor.xy)) {cor.xy <- 0}
    zero <- sum(xvec == 0 & yvec == 0) / n
    initial <- rep(NA,5)
    initial[4] <- s2.x /xbar
    initial[5] <- s2.y /ybar
    initial[6] <- zero
    initial[2:3] <- c(xbar,ybar)/initial[4:5]
    initial[1] <- min(initial[2:3]) * abs(cor.xy)
    initial[2:3] <-  initial[2:3] - initial[1]
    initial <- pmax(initial, 1e-5)
    print(initial) ###
  }
  
  booster <- function (param.matrix, xvec, yvec, n.cand = 10) {
    param.matrix[,6] <- qlogis(param.matrix[,6])  # logit transformation for probs
    param.matrix[,1:5] <- log(param.matrix[,1:5])  # log transformation for positives
    a <- param.matrix[1,]
    b <- param.matrix[5,]
    candidate <- matrix(b, byrow=TRUE, ncol=6, nrow = n.cand)
    index <- which((abs(b-a) > 1e-5) & is.finite(b) & is.finite(a))  # target param for grid search
    
    for (s in 1:n.cand) {
      candidate[s,index] <- b[index] + (b[index] - a[index]) * (s-1)
    }
    candidate[,6] <- plogis(candidate[,6])  # back-transformation
    candidate[,1:5] <- exp(candidate[,1:5])  # back-transformation for probs
    #print(candidate[,1:4])  #debug
    
    lik <- sapply(1:n.cand, function(s) {lik.BvZINB3(xvec, yvec, candidate[s,])})
    return(cbind(candidate,lik))
  }
  
  iter = 0
  boost = 0
  index = 1 # previous boosting index
  param = initial
  if (showFlag) {print(c("iter", "a0", "a1", "a2", "b1", "b2", "pi"))}
  repeat {
    iter = iter + 1
    #print(c(param))
    #print(lik(vec, pp=param[1], m0=param[2], m1=param[3], m2=param[4])) # debug
    param.old <- param # saving old parameters
    
    # updating
    expt <- dBvZINB3.Expt.vec(xy.reduced$x, xy.reduced$y, 
                              a0 = param[1], a1 = param[2], a2 = param[3], b1 = param[4], b2 = param[5], pp = param[6])
    expt <- as.vector(expt %*% xy.reduced$freq / n)
    # loglik = expt[1] * n
    delta <- expt[9] / (expt[2] + expt[4])                   # delta = E(V) / (E(xi0 + xi2))
    param[6] = expt[8]                                                    # pi = E(Z)
    
    opt.vec <- function(par.ab) {
      par.ab <- exp(par.ab)
      r1 <- sum(expt[2:4]) - sum(par.ab[1:3]) * par.ab[4]
      r2 <- expt[5:7] - digamma(par.ab[1:3]) - log(par.ab[4])
      # print(c(r1,r2)) ###
      return(c(r1,r2))
    }
    param.l <- log(param)
    result <- try(multiroot(opt.vec, start=param.l[1:4])$root, silent=TRUE)
    if (class(result)=="try-error") {
      initial = rep(1,5)
      result <- multiroot(opt.vec, start=initial[1:4], rtol=1e-20)$root
    }
    param[1:4] <- exp(result)
    param[5]    <- param[4] * delta                                       # b2
    param[6] = expt[8]
    
    
    # boosting
    # boosting
    if (iter == 6 + boost*5) {  # Creating an empty matrix
      param.boost <- matrix(NA, nrow = 5, ncol = 6)
    }
    if (iter >= 6 + boost*5 & iter <= 10 + boost*5 ) {  # Storing last ten params
      param.boost[iter - (5 + boost*5),] <- param
    }
    if (iter == 10 + boost*5) {
      param.boost <- booster(param.boost, xvec, yvec, n.cand = max(5, index * 2))
      index <- which.max(param.boost[,7])
      if (showFlag) {print(param.boost)}
      if (showFlag) {print(paste0("Jump to the ",index, "th parameter"))}
      param <- param.boost[index,1:6]
      boost <- boost + 1
    }
    
    
    #print (expt) #####
    if (showFlag) {print(c(iter, round(param,5), expt[1] * n))} #lik: lik of previous iteration
    if (max(abs(param - param.old)) <= tol) {
      names(param) <- c("a0", "a1", "a2", "b1", "b2", "pi")
      return(param)
      break
    }
  }
  #result <- data.frame(a0 = param[1], a1 = param[2], a2 = param[3], b1 = param[4], b2 = param[5], pi = param[6])
  #return(result)
}
ML.BvZINB3.2(extractor(11), extractor(16),showFlag=TRUE) #>-804.44 >10mins
ML.BvZINB3.2(extractor(1), extractor(5),showFlag=TRUE)




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
ML.BvNB3(extractor(1), extractor(12), showFlag=TRUE)
ML.BvNB3(extractor(51), extractor(12), showFlag=TRUE)

tt(1)
MLE.Geneset1$BvNB3 <- pairwise.MLE(data[iset.Mm.c2[[1]],], ML.fun = ML.BvNB3, showFlag=TRUE)  ## 1.8hrs
tt(2)

ML.BvNB2(0,155)



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
