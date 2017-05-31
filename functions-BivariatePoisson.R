##########################################################################################
## 0. Common functions
##########################################################################################

pairwise.MLE <- function(data, ML.fun, nondata.col = 1, rm.id.vec = TRUE, ...) {
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
    x <- data[s[1],-1]
    y <- data[s[2],-1]
    return(ML.fun(xvec = x, yvec = y))
  })
  MLE <- do.call(rbind, MLE)
  
  comb <- cbind(comb, MLE)
  return(comb)
}
if (FALSE) { # example
  head(pairwise.MLE(data=data[iset.Mm.c2[[1]],], ML.fun = ML.BP))
}

##########################################################################################
## 1. Bivariate Poisson
##########################################################################################

dBP <- function(x, y = NULL, m0, m1, m2, log = FALSE) {
  if (length(x) == 2) {
    y <- x[2]
    x <- x[1]
  }
  f1 <- dpois(x, m1, log = TRUE)
  f2 <- dpois(y, m2, log = TRUE)
  m <- min(x,y); p <- m0/m1/m2
  fun.a <- function(x, y, s, p, adj) {
    ifelse (p == 0, ifelse(s==0, 1, 0), exp(lchoose(x,s)+lchoose(y,s)+lfactorial(s)+ s*log(p) - adj)) 
  }
  if (max(x,y) > 100) {adj = 300} else {adj = 0}  # Handle numerical error for large numbers
  if (max(x,y) > 200) {result <- ifelse(log, -Inf, 0)} else {
    f3 <- log(sum(sapply(0:m, function(s) fun.a(x=x, y=y, s=s, p=p, adj=adj)))) + adj
    result <- f1 + f2 - m0 + f3
    if (!log) {result <- exp(result)}
  }
  if (!is.finite(ifelse(log,exp(result),result))) {result <- 0}
  return(result)
}
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
  lik <- function(a) (-sum(apply(vec, 1, dBP, m0 = a, m1 = bar[1] - a, m2 = bar[2] - a, log = TRUE)))
  # result <- nlminb(m, lik, lower = 0)
  if (m == 0) {  # when either is all zero, then mu0 is automatically 0.
    result = data.frame(par = 0, value = NA)  # lik not actually NA but can be obtained if necessary
  } else {
    result <- optim(m, lik, lower = 0, upper = m, method="Brent")
  }
  result <- data.frame(mu0 = result$par, mu1 = bar[1] - result$par, mu2 = bar[2] - result$par, likelihood = result$value)
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
MI.BP <- function(m0, m1, m2, summation = 100) {
  sum(
    sapply(1:summation, function(y)
      sum(
        sapply(1:summation, function(x) {
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

# old version of MLE
ML.BP.old <- function(xvec, yvec, tol = 1e-6) {
  vec <- cbind(xvec,yvec)
  bar <- apply(vec,2,mean)
  m <- min(bar)
  maxBP <- function(start, end, length.out = 10) {
    candidate <- seq(start, end, length.out = length.out)
    lik <- sapply(candidate, 
                  function(a) (sum(apply(vec, 1, dBP, m0 = a, m1 = bar[1] - a, m2 = bar[2] - a, log = TRUE))))
    amax <- candidate[(lik == max(lik))][1]
    vmax <- max(lik)
    return(data.frame(argmax = amax, maxvalue = vmax))
  }
  iter = 0
  start = 0; end = m - 1e-10
  repeat {
    iter = iter + 1
    a.new <- maxBP(start, end, length.out = 10)
    start <- a.new$argmax * 2/3
    end <- a.new$argmax * 4/3
    if (abs(a.new$argmax) < 0.1) {
      end <- a.new$argmax + 10^(-iter)  
    }
    print(a.new); print(start); print(end)
    if (abs(start - end) < tol) {break}
    if (iter >= 10) {break}
  }
  return(data.frame(argmax = a.new$argmax, likelihood = a.new$maxvalue))
}


##########################################################################################
## 2. Bivariate ZIP
##########################################################################################

dBZIP <- function(x, y = NULL, pp, m0, m1, m2, log = FALSE) {
  fxy <- (1-pp) * dBP (x=x, y=y, m0 = m0, m1 = m1, m2 = m2) + ifelse((x == 0 & y == 0), pp,0)
  if (log) {fxy <- log(fxy)}
  return(fxy)
}


# throw this away later
ML.BZIP.old <- function(xvec, yvec, tol = 1e-6) { #MLE based on score equations : fail (not convex)
  len <- length(xvec)          # n
  vec <- cbind(xvec,yvec)
  vec.s <- vec[apply(vec,1,sum) != 0,]  # nonzero pair (at least one is nonzero)
  len.s <- dim(vec.s)[2]       # n - n* (number of nonzero pair)
  
  
  fun.a <- function(x, y, s, p) {
    if (p == 0) {
      return(data.frame(num = 0, den = ifelse(s==0, 1, 0)))
    } else {
      den = exp(lchoose(x,s)+lchoose(y,s)+lfactorial(s)+ s*log(p))
      num = den*s
      return(data.frame(num = num, den = den))
    }
  }
  fun.b <- function(vec, p) {
    if (length(x) == 2) {
      y <- x[2]
      x <- x[1]
    }
    m <- min(x, y)
    a <- as.data.frame(sapply(0:m, function(s) fun.a(x, y, s, p)))
    a[,] <- unlist(a)
    b <- apply(a, 1, sum)
    b <- b[1] / b[2]
    return(b)
  }
  fun.c <- function(vec, p) {    # getting B
    return(sum(apply(vec, 1, fun.b, p = p)))
  }
  
  
  fun.d <- function(m.vec, vec.s, len, len.s = dim(vec.s)[1]) {
    p <- exp(sum(log(m.vec)*c(1,-1,-1)))
    f00 <- exp(-sum(m.vec)) #%>% print
    pp <- 1 - (len-len.s)/len/(1-f00)  # pi = 1- n*/(n (1-f00))
    RHS <- len*(1-pp)   # n (1-pi)   
    vec.sum <- c(0,apply(vec.s,2,sum))
    B <- fun.c(vec.s, p)
    #print(c(B, f00))
    m.vec.new <- (vec.sum + c(1,-1,-1) * B) / len/(1-pp)
    return(list(sum(abs((m.vec.new - m.vec))), m.vec.new))
  }
  
  fun.e <- function(m.vec, vec.s, len, len.s = dim(vec.s)[1], tol = 1e-6) {
    if (is.null(m.vec)) {m.vec <- c(0,0,0)}
    iter = 0
    repeat {
      iter = iter + 1
      result <- fun.d(m.vec = m.vec, vec.s = vec.s, len = len, len.s = len.s)
      print(c("iter:",iter, round(result[[2]],6), round(result[[1]],6)))
      if (result[[1]] <= tol) {
        return(result)
        break
      } else if (iter >= 100) {
        return(NA)
        break
      }
      m.vec <- result[[2]]
    }
  }
  
  result <- optim(m, lik, lower = 0, upper = m, method="Brent")
  require(optimr)
  optimr(par = c(0.000001, 0.1, 0.2), fun.d, lower = c(0,0,0), vec.s = vec.s, len = 800)
  aa <- .9999999; fun.d(c(1e-20, 552/(34+800*aa),1/(34+800*aa)), vec.s, 800)
  
  # result <- nlminb(m, lik, lower = 0)
  result <- data.frame(mu0 = result$par, mu1 = bar[1] - result$par, mu2 = bar[2] - result$par, likelihood = result$value)
  rownames(result) <- NULL
  return(result)
}


ML.BZIP <- function(xvec, yvec, tol = 1e-8, initial = NULL, showFlag = FALSE) { #MLE based on score equations : fail (not convex)
  len <- length(xvec)          # n
  vec <- cbind(xvec, yvec)
  vec.s <- vec[apply(vec,1,sum) != 0,]  # nonzero pair (at least one is nonzero)
  len.s <- dim(vec.s)[1]       # n - n* (number of nonzero pair)
  
  fun.a <- function(x, y, s, p) {
    if (p == 0) {
      return(data.frame(num = 0, den = ifelse(s==0, 1, 0)))
    } else {
      den = exp(lchoose(x,s)+lchoose(y,s)+lfactorial(s)+ s*log(p))
      num = den*s
      return(data.frame(num = num, den = den))
    }
  }
  fun.b <- function(vec, p) {
    if (length(x) == 2) {
      y <- x[2]
      x <- x[1]
    }
    m <- min(x, y)
    a <- as.data.frame(sapply(0:m, function(s) fun.a(x, y, s, p)))
    a[,] <- unlist(a)
    b <- apply(a, 1, sum)
    b <- b[1] / b[2]
    return(b)
  }
  fun.c <- function(vec, p) {    # getting B
    return(sum(apply(vec, 1, fun.b, p = p)))
  }
  lik <- function(vec, pp, m0, m1, m2) {
    return(sum(apply(vec, 1, function(x) {dBZIP(x=x[1], y=x[2], pp=pp, m0 = m0, m1 = m1, m2 = m2, log = TRUE)})))    
  }
  lik.profile <- function(vec, pp, m0, xi1, xi2) {
    return(lik(vec=vec, pp = pp, m0 = m0, m1 = xi1 - m0, m2 = xi2 -m0))
  }  
  fun.d <- function(param, vec.s, len, len.s = dim(vec.s)[1]) {
    # xi1 = mu0 + mu1, xi2 = mu0 + mu2
    sum.x.y <- apply(vec.s,2,sum)
    if (sum(sum.x.y) ==0 ) { return(data.frame(pp = 1, mu0 = 0, mu1 = 0, mu2 = 0))} # everything is zero ==> nonestimable, set pi = 0
    # EM algorithm
    iter = 0
    if (showFlag) {print(c("iter", "pi", paste0("mu",0:2)))}
    if (showFlag) {print(c(iter, round(param,5)))}
    repeat {
      iter = iter + 1
      param.old <- param # saving old parameters
      xi <- sum.x.y / len / (1-param[1])
      if (min(xi) == 0 ) {result <- data.frame(par = 0)} else {
        result <- optim(par = param[2], lik.profile, vec = vec, pp = param[1], xi1 = xi[1], xi2 = xi[2], lower = 0, upper = min(xi), method="Brent")  
      }
      param[2] <- result$par
      param[c(3,4)] <- xi - param[2]
      p <- exp(sum(log(param[2:4])*c(1,-1,-1)))
      f00 <- exp(-sum(param[2:4])) #%>% print
      if (len == len.s | len.s/len/(1-f00) > 1) { param[1] <- 0} else {
        param[1] <- 1 - len.s/len/(1-f00)  
      }
      if (len == len.s | param[1] < 0) { param[1] <- 0}
      if (showFlag) {print(c(iter, round(param,5)))}
      if (max(abs(param - param.old)) <= tol) {
        param <- data.frame(matrix(param,1,4))
        names(param) <- c("pp", "mu0", "mu1", "mu2")
        return(param)
        break
      } else if (iter >= 100) {
        param <- data.frame(matrix(NA,1,4))
        names(param) <- c("pp", "mu0", "mu1", "mu2")
        return(param)
        break
      }  
    }
  }
  if (is.null(initial)) { # when initial(starting point) is not provided
    initial <- rep(0,4)
    initial[1] <- 1 - len.s/len # inital guess for pi
    initial[3:4] <- apply(vec.s, 2, sum)/len.s
    initial[2] <- min(initial[3:4])
    initial[3:4] <- initial[3:4] - initial[2]
  }
  result <- fun.d(param = initial, vec.s = vec.s, len = len, len.s = len.s)
  return(result)
}
if (FALSE) {# example
  ML.BZIP(vec[,1], vec[,2], showFlag=TRUE, initial = rep(0,4), tol = 1e-10)
  ML.BZIP(vec[,1], vec[,2], showFlag=TRUE, tol = 1e-8)
  ML.BZIP(extractor(1), extractor(2), initial = rep(0,4), showFlag=TRUE)
}

# Numerical calcualtion of MI for bivariate Poisson
MI.BZIP <- function(pp, m0, m1, m2, summation = 100) {
  sum(
    sapply(1:summation, function(y)
      sum(
        sapply(1:summation, function(x) {
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

moment.BZIP <- function(pp, m0, m1, m2) {
  MEAN.BP <- (m0 + c(m1, m2))
  VAR.BP <- diag(MEAN.BP); VAR.BP[c(2,3)] <- m0
  MEAN <- (1-pp)*MEAN.BP
  VAR <- (1-pp)*VAR.BP + pp*(1-pp) * MEAN.BP %o% MEAN.BP
  COR <- VAR[2] / sqrt(prod(diag(VAR)))
  return(list(mean = MEAN, var = VAR, cor = COR))
}

MI.ML.BZIP <- function(xvec, yvec, tol = 1e-6, summation = 100,...) {
  MLE <- ML.BZIP (xvec = xvec, yvec = yvec,...)
  MI <- MI.BZIP (pp = MLE[1], m0 = MLE[2], m1 = MLE[3], m2 = MLE[4], summation = summation)
  return(MI)
}
if (FALSE) {# example
  MI.ML.BZIP(vec[,1], vec[,2])
  
  a <- Sys.time()
  MI.ML.BZIP(vec[,1], vec[,2]) #0.15511
  Sys.time() - a  # 12 sec
  a <- Sys.time()
  MI.ML.BP(vec[,1], vec[,2]) #0.2.839e-09
  Sys.time() - a # 5.3 sec
}

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

