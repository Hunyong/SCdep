##########################################################################################
## 0. Common functions
##########################################################################################
# remove global assign later !!!! <<- 
pairwise.MLE <- function(data, ML.fun, nondata.col = 1, p.nonzero = TRUE, rm.id.vec = TRUE, ...) {
  # data: a dataframe of which pairs are to be analyzed
  # nondata.col: the column number to be ignored for measurement eg. first col (gene names)
  # name: short name for the distribution
  data <- as.matrix(data[,-nondata.col])
  dim.p <- dim(data)[1]
  comb <- expand.grid(1:dim.p,1:dim.p)
  if (rm.id.vec)  {rw <- which(comb[,1] > comb[,2])} else {rw <- which(comb[,1] >= comb[,2])}
  comb <- data.frame(comb[rw,c(2,1)])
  # [lower.tri(matrix(1:dim.p^2,dim.p),TRUE),][,c(2,1)]
  # comb <- data.frame(t(combn(1:dim.p,2)))
  comb$pair <- apply(comb,1,function(x) paste(x,collapse="-"))
  
  # empty.result <- ML.fun(xvec = rep(0,3), yvec = rep(0,3))
  # comb <- cbind(comb, matrix(ncol = dim(empty.result)[2]))
  # names(comb) <- c(1,2,"pair", names(empty.result))
  names(comb) <- c(1,2,"pair")
  MLE <- apply(comb[,1:2], 1, function(s) {
    # if (s[1] <= 6 | s[2] <=23) {return(data.frame(matrix(NA,1,4)))} # debug #7,24 has problem
    x <- data[s[1],]
    y <- data[s[2],]
    print(s)   # debug
    tmp.sss <<- s   # debug
    tmp.aaa <<- cbind(x,y)   # debug
    return(ML.fun(xvec = x, yvec = y))
  })
  MLE <- do.call(rbind, MLE)
  comb <- cbind(comb, MLE)
  
  if (p.nonzero == TRUE) {
    n <- dim(data)[2]
    comb$non0.1 <- sapply(1:dim(comb)[1], function(i) {sum(data[comb[i,1],] != 0) / n})
    comb$non0.2 <- sapply(1:dim(comb)[1], function(i) {sum(data[comb[i,2],] != 0) / n})
    comb$non0.min <- pmin(comb$non0.1,comb$non0.2)
  }
  
  return(comb)
}
if (FALSE) { # example
  head(pairwise.MLE(data=data[iset.Mm.c2[[1]],], ML.fun = ML.BP))
}


entropy.generic <- function(..., densityfn, summation = 100) {
  if (sum(is.na(c(...)))>0) {return(data.frame(H1 = NA, H2 = NA, H12 = NA, MI = NA))} else {
    joint <- sapply(0:summation, function(x) {sapply(0:summation, function(y) {
      densityfn (x, y, ...)})})
    marginal1 <- apply(joint,1,sum)
    marginal2 <- apply(joint,2,sum)
    
    H1 <- - marginal1 * log(marginal1); H1[!is.finite(H1)] <- 0; H1 <- sum(H1)
    H2 <- - marginal2 * log(marginal2); H2[!is.finite(H2)] <- 0; H2 <- sum(H2)
    H12 <- - joint * log(joint);        H12[!is.finite(H12)] <- 0; H12 <- sum(H12)
    M12 <- max(- H12 + H1 + H2,0)
    return(data.frame(H1 = H1, H2 = H2, H12 = H12, MI = M12)) 
  }
}  
if (FALSE) {
  # two fns give the same results
  entropy.BZIP(.02,1,1,1)
  entropy.generic(.02,1,1,1,densityfn=dBZIP)
  entropy.generic(c(.97,.01,.01,.01,1,1,1),densityfn=dBZIP.B)
}

# binary profile function: for ML.BZIP.B
bin.profile <- function(xvec, yvec) {
  xvec[xvec != 0] = 1
  yvec[yvec != 0] = 1
  vec <- cbind(xvec,yvec)
  a <- rep(0,4)
  a[1] <- sum(apply(vec,1,prod)) # 1,1
  a[2] <- sum(xvec) - a[1]       # 1,0
  a[3] <- sum(yvec) - a[1]       # 0,1
  a[4] <- length(xvec) - sum(a)  # 0,0
  return(a)
}



##########################################################################################
## 1. Bivariate Poisson
##########################################################################################

dBP <- function(x, y = NULL, m0, m1, m2, log = FALSE, max = 500) {
  # max = 500: when counts exceed the max, p = 0
  if (length(x) == 2) {
    y <- x[2]
    x <- x[1]
  }
  if (m1 * m2 == 0) {  # previous code doesn't consider when mu1, mu2, or both is 0
    if ((x - y) %*% (m1 - m2) >= 0) {
      m <- min(x,y); s <- x + y - 2*m; mm <- max(m1, m2)
      result <- dpois (m, m0) * dpois (s, mm)
      # print(c("m", m, "s", s, "mm", mm, "result", result))
    } else { result <- 0 }
    if (log) {result <- log(result)}
    return(result)
  }
  f1 <- dpois(x, m1, log = TRUE)
  f2 <- dpois(y, m2, log = TRUE)
  m <- min(x,y); p <- m0/m1/m2; if (m == 0) {p <- 0} # in order not to make p NaN
  fun.a <- function(x, y, s, p, adj) {
    ifelse (p == 0, ifelse(s==0, 1, 0), exp(lchoose(x,s)+lchoose(y,s)+lfactorial(s)+ s*log(p) - adj)) 
  }
  if (max(x,y) > 100) {adj = 500} else {adj = 0}  # Handle numerical error for large numbers
  if (max(x,y) > max) {result <- ifelse(log, -Inf, 0)} else {
    f3 <- log(sum(sapply(0:m, function(s) fun.a(x=x, y=y, s=s, p=p, adj=adj)))) + adj
    result <- f1 + f2 - m0 + f3
    if (!log) {result <- exp(result)}
  }
  if (!is.finite(ifelse(log,exp(result),result))) {result <- 0}
  return(result)
}
dBP.vec <- Vectorize(dBP)
# sample data
x <- as.numeric(data[iset.Mm.c2[[1]][1],-1])
y <- as.numeric(data[iset.Mm.c2[[1]][27],-1])  # cor(x,y) = -0.004 (p=0.9)
z <- as.numeric(data[iset.Mm.c2[[1]][2],-1])   # cor(x,z) = 0.107 (p=0.002)
z2 <- as.numeric(data[iset.Mm.c2[[1]][10],-1])   # cor(x,z2) = 0.46 (p=0.000)
if (FALSE) {# example
  dBP(2, 2, 0, 1, 1)
  dBP(100, 100, 0, 1, 1)
  dBP(200, 200, 0, 1, 1)
  dBP(300, 300, 0, 1, 1)
  
  # check if marginal of BP = poi
  sum(sapply(0:100, function(a) dBP(a,2,0,2,2)))
  dpois(2,2)
  
  sum(sapply(0:100, function(a) dBP(a,2,1,2,2)))
  dpois(2,3)
  
  # joint(BP) > pois^2 (around diagonal)
  dBP(2,2,1,2,2)
  dpois(2,3)^2
  
  # joint(BP) < pois*pois (off-diagonal)
  dBP(1,3,1,2,2)
  dpois(1,3)*dpois(3,3)
}

ML.BP <- function(xvec, yvec, tol = 1e-6) {
  xvec <- as.numeric(xvec); yvec <- as.numeric(yvec)
  len <- length(xvec)
  # xbar <- mean(xvec)
  # ybar <- mean(yvec)
  vec <- cbind(xvec,yvec)
  bar <- apply(vec,2,mean)
  m <- min(bar)
  lik.BP <- function(a) (-sum(apply(vec, 1, dBP, m0 = a, m1 = bar[1] - a, m2 = bar[2] - a, log = TRUE)))
  # result <- nlminb(m, lik, lower = 0)
  if (m == 0) {  # when either is all zero, then mu0 is automatically 0.
    result = data.frame(par = 0, value = NA)  # lik not actually NA but can be obtained if necessary
  } else {
    result <- optim(m, lik.BP, lower = 0, upper = m, method="Brent")
  }
  result <- data.frame(mu0 = result$par, mu1 = bar[1] - result$par, mu2 = bar[2] - result$par, likelihood = -result$value)
  rownames(result) <- NULL
  return(result)
}
if (FALSE) { # example
  ML.BP(x,y)
  ML.BP(x,z)
  set.seed(1000); a1 <- rpois(20,1); a2 <- rpois(20,2); a3 <- rpois(20,3)
  ML.BP(a1+a2, a1+a3)
}

# Numerical calcualtion of MI for bivariate Poisson
# instead of MI.BP, use entropy.BP!!
MI.BP <- function(m0, m1, m2, summation = 100) {
  if (summation < 2* max(m0 + m1, m0 + m2)) {
    summation <- 2* max(m0 + m1, m0 + m2)
  }
  sum(
    sapply(0:summation, function(y)
      sum(
        sapply(0:summation, function(x) {
          joint <- dBP(x, y, m0, m1, m2) #%>% print
          marginal1 <- dpois(x, m0 + m1) #%>% print
          marginal2 <- dpois(y, m0 + m2)
          if (joint == 0) {return(0)} else {
            return(joint * (log(joint) - log(marginal1) - log(marginal2)))
          }
        }))))
}

if (FALSE) {# example
  MI.BP(0.1, 1, 1, summation = 100)
  MI.BP(1e-7, .1, .1, summation = 100)
  
  # example
  ML.BP(x,z2)
  MI.BP(0.2311328, 0.4588672, 1.252617, summation = 100)  #0.039445
  mi.empirical.2(x,z2)  #0.1082
  0.2311328 / sqrt(0.2311328 + 0.4588672) / sqrt(0.2311328 + 1.252617)  # cor 0.2284
  cor(x,z2)  #0.4648
}

MI.ML.BP <- function(xvec, yvec, tol = 1e-6, summation = 100) {
  MLE <- ML.BP (xvec = xvec, yvec = yvec, tol = tol)
  MI <- MI.BP (m0 = MLE$mu0, m1 = MLE$mu1, m2 = MLE$mu2, summation = summation)
  return(MI)
}
cor.ML.BP <- function(xvec, yvec, cor = FALSE, tol = 1e-6, summation = 100) {
  MLE <- ML.BP (xvec = xvec, yvec = yvec, tol = tol)
  return(with(MLE, {mu0 / sqrt((mu0 + mu1)*(mu0 + mu2))}))
}

entropy.BP <- function(m0, m1, m2, summation = 100) {
  if (summation < 2* max(m0 + m1, m0 + m2)) {
    summation <- 2* max(m0 + m1, m0 + m2)
  }
  marginal1 <- sapply(0:summation, function(s) {dpois(s, m0 + m1)})
  marginal2 <- sapply(0:summation, function(s) {dpois(s, m0 + m2)})
  joint <- sapply(0:summation, function(x) {sapply(0:summation, function(y) {dBP(x = x, y = y, m0 = m0, m1 = m1, m2 = m2)})})
  
  H1 <- - marginal1 * log(marginal1); H1[!is.finite(H1)] <- 0; H1 <- sum(H1)
  H2 <- - marginal2 * log(marginal2); H2[!is.finite(H2)] <- 0; H2 <- sum(H2)
  H12 <- - joint * log(joint);        H12[!is.finite(H12)] <- 0; H12 <- sum(H12)
  M12 <- max(- H12 + H1 + H2,0)
  return(data.frame(H1 = H1, H2 = H2, H12 = H12, MI = M12))
}  
entropy.BP.vec <- Vectorize(entropy.BP)
if (FALSE) {  # entropy.BP is faster than MI.BP
  a <- Sys.time(); entropy.BP(1,1,1); b <- Sys.time(); MI.BP(1,1,1); Sys.time()- b; b-a
}

cor.BP <- function(m0, m1, m2) {
  return(m0 / sqrt((m0 + m1)*(m0 + m2)))
}
cor.BP.vec <- Vectorize(cor.BP)





##########################################################################################
## 2. Bivariate ZIP
##########################################################################################

## 2.1 basic functions for BvZIP
dBZIP <- function(x, y = NULL, pp, m0, m1, m2, log = FALSE) {
  fxy <- (1-pp) * dBP (x=x, y=y, m0 = m0, m1 = m1, m2 = m2) + ifelse((x == 0 & y == 0), pp,0)
  if (log) {fxy <- log(fxy)}
  return(fxy)
}
dBZIP.vec <- Vectorize(dBZIP)
lik.BZIP <- function(x, y, param, pp = NULL, m0  = NULL, m1  = NULL, m2  = NULL){
  if (is.null(pp)|is.null(m0)|is.null(m1)|is.null(m2)) {
    pp = param[1]; m0 = param[2]; m1 = param[3]; m2 = param[4]
  } 
  sum(dBZIP.vec(x, y, pp = pp, m0 = m0, m1 = m1, m2 = m2, log=TRUE))
}
rBZIP <- function(n, pp, m0, m1, m2) {
  z <- rbinom(n,1,pp)
  common.tmp <- rpois(n,m0)
  x.tmp <- rpois(n,m1)
  y.tmp <- rpois(n,m2)
  x <- (common.tmp + x.tmp)*(1-z)
  y <- (common.tmp + y.tmp)*(1-z)
  return(data.frame(x=x, y=y))
}
rBZIP.vec <- Vectorize(rBZIP)

## 2.2 param estimation
# formal EM algorithm
ML.BZIP <- function(xvec, yvec, tol = 1e-6, initial = NULL, showFlag = FALSE) { #MLE based on score equations : fail (not convex)
  # counter, likelihood, param.prev for recursive function
  len <- length(xvec)          # n
  vec <- data.frame(xvec=xvec, yvec=yvec)
  sum.x.y <- apply(vec,2,sum)  # mu0+mu1, mu0+mu2
  if (sum(sum.x.y) ==0 ) { return(data.frame(pp = 0, mu0 = 0, mu1 = 0, mu2 = 0))} # everything is zero ==> nonestimable, set pi = 0
  
  # E-step
  fun.cond.exp <- function(x, y, pp, m0, m1, m2) {
    if (x+y == 0) {
      cond.prob <- pp/(pp+(1-pp)*exp(-m0-m1-m2))
      cond.expt <- c(1, m0, m1, m2) * cond.prob
    } else {
      m = min(x,y)
      prob <- sapply(0:m, function(u) {dpois(x = u, lambda = m0) * dpois(x = x - u, lambda = m1) * dpois(x = y - u, lambda = m2)})
      # print(prob); print(prob==0) #debug
      if (sum(prob) == 0) {
        # using the ratio of two probs
        if (m==0) {prob <- 1} else {
          prob <- sapply(0:(m-1), function(u) (m0/m1/m2*(x-u)*(y-u)/(u+1)))
          prob <- cumprod(c(1,prob))
          #prob <- sapply(0:m, function(u) {dBP(x = u, y=0, m0 = 0, m1 = m0, m2 = 0) * dBP(x = x-u, y=0, m0 = 0, m1 = m1, m2 = 0) * 
          #    dBP(x = 0, y= y-u, m0 = 0, m1 = 0, m2 = m2)*10^100}) 
          # used dBP instead of simply dpois to handle large numbers (when mu=1, x = 172, dpois=0 not small number)
          # and adjust by 10^100 (num/den cancel out adj.)
        }
        #print(c(1,prob)) # debug
      #print(prob)
      }
      EU <- sum((0:m)*prob)/sum(prob)
      cond.expt <- c(0, 0, x, y) + EU *c(0, 1, -1, -1)
    }
    return(cond.expt)
  }
  fun.cond.exp <- Vectorize(fun.cond.exp)
  
  # M-step 
  param.update <- function(x, y, pp, m0, m1, m2) {
    result <- fun.cond.exp(x = x, y = y, pp = pp, m0 = m0, m1 = m1, m2 = m2)
    return(apply(result,1,mean))
  }
  
  # initial guess
  if (is.null(initial)) { # when initial(starting point) is not provided
    initial <- c(0.5,1,1,1)
  }
  
  # Repeat
  iter = 0
  param = initial
  if (showFlag) {print(c("iter", "pi", paste0("mu",0:2)))}
  repeat {
    iter = iter + 1
    #print(c(param))
    #print(lik(vec, pp=param[1], m0=param[2], m1=param[3], m2=param[4])) # debug
    param.old <- param # saving old parameters
    param <- param.update (x = xvec, y = yvec, pp = param[1], m0 = param[2], m1 = param[3], m2 = param[4])
    if (showFlag) {print(c(iter, round(param,5), lik.BZIP(xvec, yvec, param)))}
    if (max(abs(param - param.old)) <= tol) {
      param <- data.frame(matrix(param,1,4))
      names(param) <- c("pp", "mu0", "mu1", "mu2")
      return(param)
      break
    }
  }
  return(param)
}

  # ML.BZIP examples
  if (FALSE) {
    ML.BZIP(vec[,1], vec[,2], showFlag=TRUE, initial = rep(0,4), tol = 1e-10)
    ML.BZIP(vec[,1], vec[,2], showFlag=TRUE, tol = 1e-8)
    ML.BZIP(extractor(1), extractor(2), initial = rep(0,4), showFlag=TRUE)
    ML.BZIP(rep(0,100), rep(0,100), initial = rep(0,4), showFlag=TRUE)
    # For all 0 pairs, identifiability issue: same likelihood for (pi=1) and (pi=0)
    ML.BZIP(c(1,rep(0,99)), rep(0,100), initial = rep(0,4), showFlag=TRUE)
    ML.BZIP(c(1,rep(0,799)), rep(0,800), initial = rep(0,4), showFlag=TRUE)
    ML.BZIP(c(1,rep(0,799)), rep(0,800), initial = c(799/800,0,0,0), showFlag=TRUE) # not converging. initial should be all zero
  }

  # unstable BZIP ML estimation
  if (FALSE) {
    ML.BZIP(tmp.1[,1],tmp.1[,2],showFlag=T) # 0.085183 0.3744193 0.04096429 0.08468888
    lik.BZIP(tmp.1[,1],tmp.1[,2], pp=.08511, m0 = .3744, m1 = .04096, m2=.0847) #lik = -55
    lik.BZIP(tmp.1[,1],tmp.1[,2], pp=.5, m0 = 1, m1 = .1, m2=.1) #lik = -59
    # far away from truth
    
  }

  # simulation to check if ML works okay
  if (FALSE) {
    set.seed(1); tt(1)
    result <- data.frame(pp= NA, mu0 = NA, mu1 = NA, mu2 = NA)
    for (i in 1:30) {
      print (c(i, "/ 30"))
      tmp.1 <- rBZIP(50, .5, 1, 1, 1)
      result[i,] <- ML.BZIP(tmp.1[,1], tmp.1[,2])
      print(result[i,])
    }
    result; tt(2)
    apply(result,2,mean) # 0.503, 0.977, 1.023, 1.044
    
    set.seed(1); tt(1)
    result2 <- result3 <- data.frame(pp= NA, mu0 = NA, mu1 = NA, mu2 = NA)
    for (i in 1:30) {
      print (c(i, "/ 30"))
      tmp.1 <- rBZIP(50, .5, 1, .1, .1)
      result2[i,] <- ML.BZIP(tmp.1[,1], tmp.1[,2])
      result3[i,] <- ML.BZIP(tmp.1[,1], tmp.1[,2], initial = c(.5,1,.1,.1))
      print(result2[i,])
      print(result3[i,])
    }
    result2;result3; tt(2) #1.7 mins
    apply(result2,2,mean,na.rm=T) # 0.4727396 0.9906099 0.1171386 0.1099863 
    apply(result3,2,mean,na.rm=T) # 0.4727392 0.9906087 0.1171385 0.1099862  okay!
    
    lik.BZIP(tmp.1[,1], tmp.1[,2], c(0.62369159, 0.04185936, 1.45443309, 1.44852183 ))
    lik.BZIP(tmp.1[,1], tmp.1[,2], c(0.62369159, 1, .1, .1 ))
    ML.BZIP(tmp.1[,1],tmp.1[,2],showFlag=T)
  }

  # Numerical calcualtion of MI for bivariate Poisson 
  ### !!! Something's wrong with this code !!!
  MI.BZIP <- function(pp, m0, m1, m2, summation = 100) {
    sum(
      sapply(0:summation, function(y)
        sum(
          sapply(0:summation, function(x) {
            joint <- dBZIP(x, y, pp, m0, m1, m2) #%>% print
            marginal1 <- dpois(x, m0 + m1) * (1-pp) #%>% print
            marginal2 <- dpois(y, m0 + m2) * (1-pp)
            if (x==0 & y==0) {
              marginal1 <- marginal1 + pp  # ZIP model
              marginal2 <- marginal2 + pp  # ZIP model
            }
            if (joint == 0) {return(0)} else {
              return(joint * (log(joint) - log(marginal1) - log(marginal2)))
            }
          }))))
  }
  ### !!! This code is also suspended until MI.BZIP is corrected !!!
  MI.ML.BZIP <- function(xvec, yvec, tol = 1e-6, summation = 100,...) {
    MLE <- ML.BZIP (xvec = xvec, yvec = yvec,...)
    MI <- MI.BZIP (pp = MLE[1], m0 = MLE[2], m1 = MLE[3], m2 = MLE[4], summation = summation)
    return(MI)
  }
  # MI.ML.BZIP examples
  if (FALSE) {# example
    MI.ML.BZIP(vec[,1], vec[,2])
    
    a <- Sys.time()
    MI.ML.BZIP(vec[,1], vec[,2]) #0.15511
    Sys.time() - a  # 12 sec
    a <- Sys.time()
    MI.ML.BP(vec[,1], vec[,2]) #0.2.839e-09
    Sys.time() - a # 5.3 sec
  }
  
  
moment.BZIP <- function(pp, m0, m1, m2) {
  MEAN.BP <- (m0 + c(m1, m2))
  VAR.BP <- diag(MEAN.BP); VAR.BP[c(2,3)] <- m0
  MEAN <- (1-pp)*MEAN.BP
  VAR <- (1-pp)*VAR.BP + pp*(1-pp) * MEAN.BP %o% MEAN.BP
  COR <- VAR[2] / sqrt(prod(diag(VAR)))
  return(list(mean = MEAN, var = VAR, cor = COR))
}


## 2.3 MI estimation
entropy.BZIP <- function(pp, m0, m1, m2, summation = 100) {
  if (is.na(pp) |is.na (m0) | is.na(m1) | is.na(m2)) {return(data.frame(H1 = NA, H2 = NA, H12 = NA, MI = NA))} else {
    if (summation < 2* max(m0 + m1, m0 + m2)) {
      summation <- 2* max(m0 + m1, m0 + m2)
    }
    joint <- sapply(0:summation, function(x) {sapply(0:summation, function(y) {dBP(x = x, y = y, m0 = m0, m1 = m1, m2 = m2)})})
    joint <- joint *(1-pp)
    joint[1,1] <- joint[1,1] + pp
    marginal1 <- apply(joint,1,sum)
    marginal2 <- apply(joint,2,sum)
    
    H1 <- - marginal1 * log(marginal1); H1[!is.finite(H1)] <- 0; H1 <- sum(H1)
    H2 <- - marginal2 * log(marginal2); H2[!is.finite(H2)] <- 0; H2 <- sum(H2)
    H12 <- - joint * log(joint);        H12[!is.finite(H12)] <- 0; H12 <- sum(H12)
    M12 <- max(- H12 + H1 + H2,0)
    return(data.frame(H1 = H1, H2 = H2, H12 = H12, MI = M12)) 
  }
}  
entropy.BZIP.vec <- Vectorize(entropy.BZIP)

## 2.4 deviance
# Finding saturated parameters (ML for individual lik)
BZIP.saturated <- function(x, y) {
  m = min(x, y)
  return(param = c((x==0)*(y==0), m, x-m, y-m))
}
BZIP.saturated <- Vectorize(BZIP.saturated)
dev.BZIP <- function(xvec, yvec, param = NULL, pp = NULL, m0  = NULL, m1  = NULL, m2  = NULL) {
  # If params = NULL, apply ML.BZIP. if else, apply those params
  if (is.null (param)) { 
    if (is.null (pp) | is.null (m0) | is.null (m1) | is.null (m2) ) {
      param = ML.BZIP (xvec = xvec, yvec = yvec)
    }
    else { param = c(pp, m0, m1, m2)}
  }
  # Log-likelihood of the BZIP model
  lik.model <- lik.BZIP (x = xvec, y = yvec, param = param)
  
  # Saturated model BZIP params
  param.sat <- t(BZIP.saturated(xvec, yvec))
  lik.sat   <- sum(dBZIP.vec(x= xvec, y = yvec, pp = param.sat[,1], m0 = param.sat[,2], m1 = param.sat[,3], m2 = param.sat[,4], log = TRUE))
  
  return(data.frame(model.likelihood = lik.model, satrtd.likelihood = lik.sat, deviance = 2*(lik.sat - lik.model)))
}
dev.BP <- function(xvec, yvec, param = NULL, m0  = NULL, m1  = NULL, m2  = NULL) {
  # If params = NULL, apply ML.BZIP. if else, apply those params
  if (is.null (param)) { 
    if (is.null (m0) | is.null (m1) | is.null (m2) ) {
      param = ML.BP (xvec = xvec, yvec = yvec)
    }
    else { param = c(m0, m1, m2)}
  }
  # Log-likelihood of the BP model
  lik.model <- sum(dBP.vec(x= xvec, y = yvec, m0 = param[1], m1 = param[2], m2 = param[3], log = TRUE))
  
  # Saturated model BP params
  param.sat <- t(BZIP.saturated(xvec, yvec))[,-1] #Same as (BZIP.sat without pi)
  lik.sat   <- sum(dBP.vec(x= xvec, y = yvec, m0 = param.sat[,1], m1 = param.sat[,2], m2 = param.sat[,3], log = TRUE))
  
  return(data.frame(model.likelihood = lik.model, satrtd.likelihood = lik.sat, deviance = 2*(lik.sat - lik.model)))
}



##########################################################################################
## 3. Bivariate ZIP.B: General BZIP1 with marginal ZIP condition (6params)
##########################################################################################

## 3.1 basic functions for BvZIP.B
dBZIP.B <- function(x, y = NULL, p1, p2, p3, p4 = 1 - p1 - p2 - p3, m0, m1, m2, log = FALSE) {
  fxy <- p1 * dBP (x=x, y=y, m0 = m0, m1 = m1, m2 = m2) +
    p2 * {if (y == 0) dBP (x=x, y=y, m0 = 0, m1 = m0 + m1, m2 = 0) else 0} +
    p3 * {if (x == 0) dBP (x=x, y=y, m0 = 0, m1 = 0, m2 = m0 + m2) else 0} +
    p4 * {if (x + y == 0) 1 else 0}
  if (log) {fxy <- log(fxy)}
  return(fxy)
}
dBZIP.B.vec <- Vectorize(dBZIP.B)
lik.BZIP.B <- function(x, y, param, p1 = NULL, p2 = NULL, p3 = NULL, p4 = NULL, m0  = NULL, m1  = NULL, m2  = NULL){
  if (is.null(p1)|is.null(p2)|is.null(p3)|is.null(p4)|is.null(m0)|is.null(m1)|is.null(m2)) {
    p1 = param[1]; p2 = param[2]; p3 = param[3]; p4 = param[4]; 
    m0 = param[5]; m1 = param[6]; m2 = param[7]
  } 
  sum(dBZIP.B.vec(x, y, p1 = p1, p2 = p2, p3 = p3, p4 = p4, 
                m0 = m0, m1 = m1, m2 = m2, log=TRUE))
}
rBZIP.B <- function(n, p1, p2, p3, p4, m0, m1, m2) {
  z <- rmultinom(n,1,c(p1,p2,p3,p4))
  u.tmp <- rpois(n,m0)
  v.tmp <- rpois(n,m1)
  w.tmp <- rpois(n,m2)
  #(Xi, Yi) = ((Ei1 + Ei2)(Ui + Vi), (Ei1 + Ei3)(Ui + Wi)).
  x <- (z[1,] + z[2,]) * (u.tmp + v.tmp)
  y <- (z[1,] + z[3,]) * (u.tmp + w.tmp)
  return(data.frame(x=x, y=y))
}

## 3.2 parameter estimation
# Method of moment, for starting values of MLE
MME.BZIP.B <- function(xvec, yvec) {
  m.x <- mean(xvec); m.x2 <- mean(xvec^2)
  m.y <- mean(yvec); m.y2 <- mean(yvec^2)
  m.x2.y <- mean(xvec^2*yvec); m.y2.x <- mean(yvec^2*xvec)
  m.x.y <- mean(xvec*yvec)
  w <- (m.x2.y + m.y2.x) / (2 * m.x.y)
  lambda1 <- m.x2 / m.x - 1
  lambda2 <- m.y2 / m.y - 1
  mu0 <- (lambda1*lambda2) * (w - (2 + lambda1 + lambda2)/2) / (lambda1 + lambda2 + 1 - w)
  #print(c(w,lambda1, lambda2, mu0))
  mu0 <- max(1e-10, min(mu0, m.x, m.y))
  mu1 <- m.x - mu0
  mu2 <- m.y - mu0
  
  pi1 <- m.x.y / (lambda1*lambda2 + mu0)
  pi2 <- m.x / lambda1 - pi1
  pi3 <- m.y / lambda2 - pi1
  pi4 <- 1 - pi1 - pi2 - pi3
  pi <- c(pi1, pi2, pi3, pi4)
  pi <- pmax(0, pmin(pi,1))
  pi <- pi/sum(pi)
  return(data.frame(pi1 = pi[1], pi2 = pi[2], pi3 = pi[3], pi4 = pi[4], mu0 = mu0, mu1 = mu1, mu2 = mu2))
}

ML.BZIP.B <- function(xvec, yvec, tol = 1e-6, initial = NULL, showFlag = FALSE) { #MLE based on score equations : fail (not convex)
  # counter, likelihood, param.prev for recursive function
  len <- length(xvec)          # n
  vec <- data.frame(xvec=xvec, yvec=yvec)
  sum.x.y <- apply(vec,2,sum)  # mu0+mu1, mu0+mu2
  if (sum(sum.x.y) ==0 ) { return(data.frame(p1 = 0, p2 = 0, p3 = 0, p4 = 1,  mu0 = 0, mu1 = 0, mu2 = 0))} # everything is zero ==> nonestimable, set pi = 0
  
  # E-step
  fun.cond.exp <- function(x, y, p1, p2, p3, p4, m0, m1, m2) {
    PrA <- p1 * dBP(x=x, y=y, m0 = m0, m1 = m1, m2 = m2)
    PrB <- p2 * (y==0) * dBP(x=x, y=0, m0 =0, m1 = m0 + m1, m2 = 0)
    PrC <- p3 * (x==0) * dBP(x=0, y=y, m0 =0, m1 = 0, m2 = m0 + m2)
    PrD <- p4 * (x + y ==0)
    PrAgr <- c(PrA, PrB, PrC, PrD)
    PrSum <- sum(PrAgr)
    ### Conditional expectations (EU, EV, EW) for each profile case of (x, y)
    if (x + y ==0) {
      EU <- PrD/PrSum * m0
      EV <- (PrC + PrD)/PrSum * m1
      EW <- (PrB + PrD)/PrSum * m2
    } else
      
    if (x > 0 & y == 0) {
      prob <- sapply(0:x, function(u) {u * dpois(x = u, lambda = m0) * dpois(x = x - u, lambda = m1)})
      EU <- sum(prob) * p2 / PrSum
      EV <- x - EU
      EW <- PrB/PrSum * m2
    } else
      
    if (x == 0 & y > 0) {
      prob <- sapply(0:y, function(u) {u * dpois(x = u, lambda = m0) * dpois(x = y - u, lambda = m2)})
      EU <- sum(prob) * p3 / PrSum
      EV <- PrC/PrSum * m1
      EW <- y - EU
    } else
      
    if (x * y > 0) {
      m = min(x,y)
      prob <- sapply(0:m, function(u) {u * dpois(x = u, lambda = m0) * dpois(x = x - u, lambda = m1) * dpois(x = y - u, lambda = m2)})
      EU <- sum(prob) * p1 / PrSum
      EV <- x - EU
      EW <- y - EU
    }
    
    ### Conditional expectations (EE) for all profile cases of (x, y)
    EE <- PrAgr/PrSum
    return(c(EE,EU,EV,EW))
  }
  fun.cond.exp <- Vectorize(fun.cond.exp)
  
  # M-step 
  param.update <- function(x, y, p1, p2, p3, p4, m0, m1, m2) {
    result <- fun.cond.exp(x, y, p1, p2, p3, p4, m0, m1, m2)
    result[result < 0] <- 0
    return(apply(result,1,mean))
  }
  
  # initial guess
  if (is.null(initial)) { # when initial(starting point) is not provided
    initial <- MME.BZIP.B(xvec=xvec, yvec=yvec)
    if (sum(is.na(initial))>0) {
      initial <- bin.profile(xvec, yvec)   # freq of each zero-nonzero profile
      initial <- initial/sum(initial)      # relative freq
      initial[5:7] <- c(0.001, sum.x.y/len/(1-initial[4]))
    }
  }
  
  # Repeat
  iter = 0
  param = initial
  if (showFlag) {print(c("iter", paste0("p",1:4), paste0("mu",0:2)))}
  repeat {
    iter = iter + 1
    # print(lik(vec, pp=param[1:4], m0=param[5], m1=param[6], m2=param[7])) # debug
    param.old <- param # saving old parameters
    param <- param.update (x = xvec, y = yvec, p1 = param[1], p2 = param[2], p3 = param[3], p4 = param[4], m0=param[5], m1=param[6], m2=param[7])
    if (showFlag) {print(c(iter, round(param,5), lik.BZIP.B(xvec,yvec,param=param) ))} #lik.BZIP(xvec, yvec, param)
    if (max(abs(param - param.old)) <= tol) {
      param <- data.frame(matrix(param,1,7))
      names(param) <- c("p1", "p2", "p3", "p4", "mu0", "mu1", "mu2")
      return(param)
      break
    }
  }
  return(param)
}

# ML.BZIP.B2: Simpler EM + boosting + initial guess manipulation (for small counts pi4=0)
ML.BZIP.B2 <- function(xvec, yvec, tol = 1e-6, initial = NULL, showFlag = FALSE, maxiter=200) { #MLE based on score equations : fail (not convex)
  len <- length(xvec)          # n
  vec <- data.frame(xvec=xvec, yvec=yvec)
  sum.x.y <- apply(vec,2,sum)  # mu0+mu1, mu0+mu2
  if (sum(sum.x.y) ==0 ) { return(data.frame(p1 = 0, p2 = 0, p3 = 0, p4 = 1,  mu0 = 0, mu1 = 0, mu2 = 0))} # everything is zero ==> nonestimable, set pi = 0
  
  # E-step
  fun.cond.exp.a <- function(x, y, p1, p2, p3, p4, m0, m1, m2) {
    PrA <- p1 * dBP(x=x, y=y, m0 = m0, m1 = m1, m2 = m2)
    PrB <- p2 * (y==0) * dBP(x=x, y=0, m0 =0, m1 = m0 + m1, m2 = 0)
    PrC <- p3 * (x==0) * dBP(x=0, y=y, m0 =0, m1 = 0, m2 = m0 + m2)
    PrD <- p4 * (x + y ==0)
    PrAgr <- c(PrA, PrB, PrC, PrD)
    PrSum <- sum(PrAgr)
    
    ### Conditional expectations (EE) for all profile cases of (x, y)
    EE <- (PrAgr/PrSum)
    
    ### Conditional expectations (EkU) for all profile cases of (x, y)
    m = min(x,y)
    num <- sum(sapply(0:m, function(u) {u * dpois(x = u, lambda = m0) * dpois(x = x - u, lambda = m1) * dpois(x = y - u, lambda = m2)}))
    den <- dBP(x=x, y=y, m0 = m0, m1 = m1, m2 = m2)
    EE1U <- EE[1] * num/den
    
    num <- sum(sapply(0:x, function(u) {u * dpois(x = u, lambda = m0) * dpois(x = x - u, lambda = m1)}))
    den <- dBP(x=x, y=0, m0 = 0, m1 = m0 + m1, m2 = 0)
    EE2U <- EE[2] * num/den
    
    num <- sum(sapply(0:y, function(u) {u * dpois(x = u, lambda = m0) * dpois(x = y - u, lambda = m2)}))
    den <- dBP(x=0, y=y, m0 = 0, m1 = 0, m2 = m0 + m2)
    EE3U <- EE[3] * num/den
    
    EE4U <- EE[4] * m0
    
    m0.den <- 1
    m0.num <- (EE1U + EE2U + EE3U + EE4U)
    m1.den <- sum(EE[1:2])
    m1.num <- m1.den * x - EE1U - EE2U
    m2.den <- sum(EE[c(1,3)])
    m2.num <- m2.den * y - EE1U - EE3U
    
    #print(c(EE,EU,EV,EW))   ####debug
    return(c(EE, m0.den, m0.num, m1.den, m1.num, m2.den, m2.num))
  }
  fun.cond.exp <- Vectorize(fun.cond.exp.a)
  
  # M-step 
  param.update <- function(x, y, p1, p2, p3, p4, m0, m1, m2) {
    result <- fun.cond.exp(x, y, p1, p2, p3, p4, m0, m1, m2)
    #print(result)  ####debug
    result[result < 0] <- 0
    result <- apply(result,1,mean)
    result2 <- result[1:4]
    result2[5] <- result[6]/result[5]
    result2[6] <- result[8]/result[7]
    result2[7] <- result[10]/result[9]
    result2[is.na(result2)] <- as.numeric(c(p1, p2, p3, p4, m0, m1, m2))[is.na(result2)]
    return(result2)
  }
  
  # initial guess
  if (is.null(initial)) { # when initial(starting point) is not provided
    initial <- MME.BZIP.B(xvec=xvec, yvec=yvec)
    if (sum(is.na(initial))>0) {
      initial <- bin.profile(xvec, yvec)   # freq of each zero-nonzero profile
      initial <- initial/sum(initial)      # relative freq
      initial[5:7] <- c(0.001, sum.x.y/len/(1-initial[4]))
    }
    if (min(sum.x.y) < 5 & min(sum.x.y) > 0) {
      initial[2:3] <- (initial[4] - 1e-10) * sum.x.y/sum(sum.x.y)
      initial[4] <- 1e-10}
    #print(initial)
  }
  
  booster <- function (param.matrix, xvec, yvec, n.cand = 10) {
    param.matrix[,1:4] <- qlogis(param.matrix[,1:4])  # logit transformation for probs
    a <- param.matrix[1,]
    b <- param.matrix[5,]
    candidate <- matrix(b, byrow=TRUE, ncol=7, nrow = n.cand)
    index <- which((abs(b-a) > 1e-5) & is.finite(b) & is.finite(a))  # target param for grid search
    
    for (s in 1:n.cand) {
      candidate[s,index] <- b[index] + (b[index] - a[index]) * 2^(s-1)
    }
    candidate[,1:4] <- plogis(candidate[,1:4])  # back-transformation for probs
    #print(candidate[,1:4])  #debug
    
    # adjusting sum-to-zero
    index.p1 <- which(1:4 %in% index)     # changing part in pi vector
    index.p2 <- which(!(1:4 %in% index))  # unchanging part
    if (length(index.p1) > 1) {
      part1.sum <- apply(candidate[, index.p1],1,sum)
      if (length(index.p2) == 0) {part2.sum <- 0 } else 
        if (length(index.p2) == 1) {part2.sum <-candidate[, index.p2]} else 
        {part2.sum <- apply(candidate[, index.p2],1,sum)}
      adj.ratio <- (1-part2.sum)/part1.sum
      #print(c(index.p1,index.p2))  #debug
      #print(cbind(part1.sum, part2.sum))  #debug
      #print(adj.ratio)  #debug
      candidate[,index.p1] <- candidate[,index.p1]*adj.ratio   # make sum-to-zero      
    } 
    
    lik <- sapply(1:n.cand, function(s) {lik.BZIP.B(xvec, yvec, candidate[s,])})
    return(cbind(candidate,lik))
  }
  
  # Repeat
  iter = 0
  boost = 0
  index = 1 # previous boosting index
  param = initial
  if (showFlag) {print(c("iter", paste0("p",1:4), paste0("mu",0:2)))}
  repeat {
    if (iter >= maxiter) { warning("EM exceeded maximum number of iterations")
      param <- data.frame(matrix(NA,1,7))
      names(param) <- c("p1", "p2", "p3", "p4", "mu0", "mu1", "mu2")
      return(param)}
    iter = iter + 1
    # print(lik(vec, pp=param[1:4], m0=param[5], m1=param[6], m2=param[7])) # debug
    param.old <- param # saving old parameters
    param <- param.update (x = xvec, y = yvec, p1 = param[1], p2 = param[2], p3 = param[3], p4 = param[4], m0=param[5], m1=param[6], m2=param[7])
    # boosting
    if (iter == 6 + boost*5) {  # Creating an empty matrix
      param.boost <- matrix(NA, nrow = 5, ncol = 7)
    }
    if (iter >= 6 + boost*5 & iter <= 10 + boost*5 ) {  # Storing last ten params
      param.boost[iter - (5 + boost*5),] <- param
    }
    if (iter == 10 + boost*5) {
      param.boost <- booster(param.boost, xvec, yvec, n.cand = max(5, index * 2))
      index <- which.max(param.boost[,8])
      if (showFlag) {print(param.boost)}
      if (showFlag) {print(paste0("Jump to the ",index, "th parameter"))}
      param <- param.boost[index,1:7]
      boost <- boost + 1
    }
    if (showFlag) {print(c(iter, round(param,5), lik.BZIP.B(xvec,yvec,param=param) ))}
    if (max(abs(param - param.old)) <= tol) {
      param <- data.frame(matrix(param,1,7))
      names(param) <- c("p1", "p2", "p3", "p4", "mu0", "mu1", "mu2")
      return(param)
      break
    }
    
  }
  return(param)
}
  # Comparing algorithms
  if (FALSE) {
    tt(1)
    ML.BZIP.B(extractor(1),extractor(2), showFlag=TRUE)   #EM (BZIP.B)
    # p1         p2          p3        p4       mu0     mu1     mu2
    # 0.002569217 0.03868079 0.006423043 0.952327 3.272441e-09 16.72727 3.61411
    # 1285 iterations, 2.67 mins (lik=-487.8919)
    # lik.BZIP.B(extractor(1),extractor(2), c(0.002569217, 0.03868079, 0.006423043, 0.952327, 3.272441e-09, 16.72727, 3.61411)) #-487.8919    
    tt(2)
    
    tt(1)
    ML.BZIP.B2(extractor(1),extractor(2), showFlag=TRUE)   #EM2 (BZIP.B)
    # p1         p2          p3        p4       mu0     mu1     mu2
    # 0.002569209 0.03868079 0.006423021 0.952327 7.816525e-09 16.72727 3.614232
    # 10 iterations, 2.45 mins  (lik=-487.8919)
    #lik.BZIP.B(extractor(1),extractor(2), c(0.002569209, 0.03868079, 0.006423021, 0.952327, 7.816525e-09, 16.72727, 3.614232)) #-487.8919
    tt(2)
    
    tt(1)
    ML.BZIP(extractor(1),extractor(2), showFlag=TRUE)      #EM (BZIP)
    #   pp          mu0      mu1       mu2
    # 0.9525 0.0001888937 14.52612 0.6840215
    # 1122 iterations, 58 sec
    tt(2)
    
    tt(1)
    ML.BZIP.old5(extractor(1),extractor(2), showFlag=TRUE)   #Direct maximization
    # pp          mu0      mu1       mu2
    # 0.9525 7.699164e-09 14.52631 0.6842103
    # 3 iterations, 0.53 sec
    tt(2)  
  }

## 3.3 entropy
entropy.BZIP.B <- function(p1, p2, p3, p4, m0, m1, m2, summation = 100) {
  if (is.na(p1) | is.na(p2) | is.na(p3) | is.na(p4) |is.na (m0) | is.na(m1) | is.na(m2)) {
    return(data.frame(H1 = NA, H2 = NA, H12 = NA, MI = NA))} else {
      if (summation < 2* max(m0 + m1, m0 + m2)) {
        summation <- 2* max(m0 + m1, m0 + m2)
      }
      joint <- sapply(0:summation, function(x) {sapply(0:summation, function(y) {
        dBZIP.B(x = x, y = y, p1 = p1, p2 = p2, p3 = p3, p4 = p4, m0 = m0, m1 = m1, m2 = m2)})})
      marginal1 <- apply(joint,1,sum)
      marginal2 <- apply(joint,2,sum)
      
      H1 <- - marginal1 * log(marginal1); H1[!is.finite(H1)] <- 0; H1 <- sum(H1)
      H2 <- - marginal2 * log(marginal2); H2[!is.finite(H2)] <- 0; H2 <- sum(H2)
      H12 <- - joint * log(joint);        H12[!is.finite(H12)] <- 0; H12 <- sum(H12)
      M12 <- max(- H12 + H1 + H2,0)
      return(data.frame(H1 = H1, H2 = H2, H12 = H12, MI = M12)) 
    }
} 
entropy.BZIP.B.vec <- Vectorize(entropy.BZIP.B)


## 3.4 deviance
# Finding saturated parameters (ML for individual lik)
BZIP.B.saturated <- function(x, y) {
  m = min(x, y)
  return(param = c((x!=0)*(y!=0),(x!=0)*(y==0),(x==0)*(y!=0),(x==0)*(y==0), m, x-m, y-m))
}
BZIP.B.saturated <- Vectorize(BZIP.B.saturated)
dev.BZIP.B <- function(xvec, yvec, param = NULL, p1 = NULL, p2 = NULL, p3 = NULL, p4 = NULL, 
                     m0  = NULL, m1  = NULL, m2  = NULL) {
  # If params = NULL, apply ML.BZIP. if else, apply those params
  if (is.null (param)) { 
    if (is.null (p1) | is.null (p2) | is.null (p3) | is.null (p4) | 
        is.null (m0) | is.null (m1) | is.null (m2) ) {
      param = ML.BZIP.B2 (xvec = xvec, yvec = yvec)
    }
    else { param = c(p1, p2, p3, p4, m0, m1, m2)}
  }
  # Log-likelihood of the BZIP model
  lik.model <- lik.BZIP.B (x = xvec, y = yvec, param = param)
  
  # Saturated model BZIP params
  param.sat <- t(BZIP.B.saturated(xvec, yvec))
  lik.sat   <- sum(dBZIP.B.vec(x= xvec, y = yvec, p1 = param.sat[,1], p2 = param.sat[,2], p3 = param.sat[,3], p4 = param.sat[,4],
                               m0 = param.sat[,5], m1 = param.sat[,6], m2 = param.sat[,7], log = TRUE))
  
  return(data.frame(model.likelihood = lik.model, satrtd.likelihood = lik.sat, deviance = 2*(lik.sat - lik.model)))
}

##########################################################################################
## 9. Under construction
##########################################################################################


##### need to update parametric correlation
cor.ML.BZIP <- function(xvec, yvec, ...) {
  MLE <- ML.BZIP (xvec = xvec, yvec = yvec,...)
  COR <- moment.BZIP (pp = MLE[1], m0 = MLE[2], m1 = MLE[3], m2 = MLE[4])$cor
  return(COR)
}

# test of correlation function
if (FALSE) {
  vec2 <- matrix(c(rep(0,100),1,1,1, rep(0,100),3,3,3), 103)
  cor(vec2[,1],vec2[,2])          # 1.00
  cor.ML.BZIP(vec2[,1],vec2[,2])  # 0.60
  cor.ML.BP(vec2[,1],vec2[,2])    # 0.58
  
  vec2[98,] <- c(2,0); vec2[99,] <- c(0,3)
  cor(vec2[,1],vec2[,2])          # 0.55
  cor.ML.BZIP(vec2[,1],vec2[,2])  # 0.83
  cor.ML.BP(vec2[,1],vec2[,2])    # 0.38
  
  vec2 <- vec2[90:103, ]
  cor(vec2[,1],vec2[,2])          # 0.41
  cor.ML.BZIP(vec2[,1],vec2[,2])  # 0.78
  cor.ML.BP(vec2[,1],vec2[,2])    # 0.31
  
  vec2 <- vec2[8:14, ]
  cor(vec2[,1],vec2[,2])          # 0.06
  cor.ML.BZIP(vec2[,1],vec2[,2])  # 0.70
  cor.ML.BP(vec2[,1],vec2[,2])    # 0.08
  table(vec2[,1], vec2[,2])
  
  vec2 <- matrix(c(0,0,1,1,1,3,3,0,0,0),5)
  cor(vec2[,1],vec2[,2])          # -1.00
  cor.ML.BZIP(vec2[,1],vec2[,2])  # 0.71   # not fitting well
  cor.ML.BP(vec2[,1],vec2[,2])    # 0.00
  table(vec2[,1], vec2[,2])
}

