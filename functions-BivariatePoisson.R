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
  m <- min(x,y); p <- m0/m1/m2; if (m == 0) {p <- 0} # in order not to make p NaN
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
#dBP <- Vectorize(dBP)
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

# basic functions for BvZIP
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


ML.BZIP.old2 <- function(xvec, yvec, tol = 1e-8, initial = NULL, showFlag = FALSE) { #MLE based on score equations : fail (not convex)
  len <- length(xvec)          # n
  vec <- data.frame(xvec=xvec, yvec=yvec)
  vec.s <- vec[apply(vec,1,sum) != 0,]  # nonzero pair (at least one is nonzero)
  len.s <- dim(vec.s)[1]       # n* (number of nonzero pair)
  
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
      # print(lik(vec, pp=param[1], m0=param[2], m1=param[3], m2=param[4])) # debug
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
      } else if (iter >= 50) {
        print("iteration exceeded")
        param <- data.frame(matrix(NA,1,4))
        names(param) <- c("pp", "mu0", "mu1", "mu2")
        return(param)
        break
      }  
    }
  }
  if (len.s == 1) { # when there is only one nonzero pair, pi should be 0
    result.a <- (ML.BP(xvec=xvec, yvec=yvec)[,1:3])
    return(data.frame(pp= 0, mu0 = result.a$mu0, mu1 = result.a$mu1, mu2 = result.a$mu2))
  } else {
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
}

# formal EM algorithm
ML.BZIP.old3 <- function(xvec, yvec, tol = 1e-8, initial = NULL, showFlag = FALSE) { #MLE based on score equations : fail (not convex)
  len <- length(xvec)          # n
  vec <- data.frame(xvec=xvec, yvec=yvec)
  vec.s <- vec[apply(vec,1,sum) != 0,]  # nonzero pair (at least one is nonzero)
  len.s <- dim(vec.s)[1]       # n* (number of nonzero pair)
  len.0 <- len - len.s         # n0
  sum.x.y <- apply(vec,2,sum)  # mu0+mu1, mu0+mu2
  
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
    if (length(vec) == 2) {
      y <- as.numeric(vec[2])
      x <- as.numeric(vec[1])
    }
    m <- min(x, y)
    # a <- as.data.frame(sapply(0:m, function(s) {print(c(x,y,s,p)); return(fun.a(x, y, s, p))})) %>% print
    a <- as.data.frame(sapply(0:m, function(s) fun.a(x, y, s, p)))
    a[,] <- unlist(a)
    b <- apply(a, 1, sum)
    b <- b[1] / b[2]
    return(b)
  }
  fun.c <- function(vec, p) {    # getting B
    return(sum(apply(vec, 1, fun.b, p = p)))
  }
  fun.c.profile <- function(vec, mu0, pp, mu.x.y) {    # solve for mu0
    p = mu0 / (mu.x.y[1] - mu0) / (mu.x.y[2] - mu0)
    return((fun.c(vec, p) / len / (1-pp) - mu0)^2)
  }
  xi.k <- function(param) {   #update xi = E(Delta_i | obs'd, kth param)
    return(param[1] / (param[1] + (1-param[1]) *exp(- param %*% c(0,1,1,1) )))
  }
  
  fun.d <- function(param, vec.s, len, len.s = dim(vec.s)[1]) {
    if (sum(sum.x.y) ==0 ) { return(data.frame(pp = 0, mu0 = 0, mu1 = 0, mu2 = 0))} # everything is zero ==> nonestimable, set pi = 0
    if (min(sum.x.y) ==0 ) { return(data.frame(pp = 0, mu0 = 0, mu1 = sum.x.y[1]/len, mu2 = sum.x.y[2]/len))} # everything is zero ==> nonestimable, set pi = 0
    # EM algorithm
    iter = 0
    if (showFlag) {print(c("iter", "pi", paste0("mu",0:2)))}
    if (showFlag) {print(c(iter, round(param,5)))}
    repeat {
      iter = iter + 1
      # print(lik(vec, pp=param[1], m0=param[2], m1=param[3], m2=param[4])) # debug
      param.old <- param # saving old parameters
      
    if (min(sum.x.y) == 0 ) {result <- data.frame(par = 0)} else {
      result <- optim(par = param[2], lik.profile, vec = vec, pp = param[1], xi1 = xi[1], xi2 = xi[2], lower = 0, upper = min(xi), method="Brent")
    }
      
      xi <- xi.k(param)
      param[1] <- len.0/len*xi
      mu.x.y <- sum.x.y / len / (1 - param[1])
      param[2] <- optim(param[2], fun.c.profile, lower = 0, upper = min(mu.x.y), vec = vec.s, pp = param[1], mu.x.y = mu.x.y, method="Brent")$par
      param[3:4] <- mu.x.y - param[2]
      
      # if (len == len.s | param[1] < 0) { param[1] <- 0}
      if (showFlag) {print(c(iter, round(param,5)))}
      if (max(abs(param - param.old)) <= tol) {
        param <- data.frame(matrix(param,1,4))
        names(param) <- c("pp", "mu0", "mu1", "mu2")
        return(param)
        break
      } else if (iter >= 50) {
        print("iteration exceeded")
        param <- data.frame(matrix(NA,1,4))
        names(param) <- c("pp", "mu0", "mu1", "mu2")
        return(param)
        break
      }  
    }
  }
  #if (len.s == 1) { # when there is only one nonzero pair, pi should be 0
  #  result.a <- (ML.BP(xvec=xvec, yvec=yvec)[,1:3])
  #  return(data.frame(pp= 0, mu0 = result.a$mu0, mu1 = result.a$mu1, mu2 = result.a$mu2))
  #} else {
  if (is.null(initial)) { # when initial(starting point) is not provided
    initial <- rep(0,4)
    initial[1] <- 1 - len.s/len # inital guess for pi
    initial[3:4] <- apply(vec.s, 2, sum)/len.s
    initial[2] <- min(initial[3:4])
    initial[3:4] <- initial[3:4] - initial[2]
  }
  if (initial[1] == 0) {
    initial[1] <- 0.01
  }
  result <- fun.d(param = initial, vec.s = vec.s, len = len, len.s = len.s)
  return(result)
}

# formal EM algorithm
ML.BZIP <- function(xvec, yvec, tol = 1e-6, initial = NULL, showFlag = FALSE) { #MLE based on score equations : fail (not convex)
  len <- length(xvec)          # n
  vec <- data.frame(xvec=xvec, yvec=yvec)
  vec.s <- vec[apply(vec,1,sum) != 0,]  # nonzero pair (at least one is nonzero)
  len.s <- dim(vec.s)[1]       # n* (number of nonzero pair)
  len.0 <- len - len.s         # n0
  sum.x.y <- apply(vec,2,sum)  # mu0+mu1, mu0+mu2
  
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
    if (length(vec) == 2) {
      y <- as.numeric(vec[2])
      x <- as.numeric(vec[1])
    }
    m <- min(x, y)
    # a <- as.data.frame(sapply(0:m, function(s) {print(c(x,y,s,p)); return(fun.a(x, y, s, p))})) %>% print
    a <- as.data.frame(sapply(0:m, function(s) fun.a(x, y, s, p)))
    a[,] <- unlist(a)
    b <- apply(a, 1, sum)
    b <- b[1] / b[2]
    return(b)
  }
  fun.c <- function(vec, p) {    # getting B
    return(sum(apply(vec, 1, fun.b, p = p)))
  }
  fun.c.profile <- function(vec, mu0, pp, mu.x.y) {    # solve for mu0
    return(-lik.BZIP(x = vec[,1], y= vec[,2], pp=pp, m0 = mu0, m1 = mu.x.y[1] - mu0, m2 = mu.x.y[2] -mu0))
  }
  xi.k <- function(param) {   #update xi = E(Delta_i | obs'd, kth param)
    return(param[1] / (param[1] + (1-param[1]) *exp(- param %*% c(0,1,1,1) )))
  }
  
  fun.d <- function(param, vec.s, len, len.s = dim(vec.s)[1]) {
    if (sum(sum.x.y) ==0 ) { return(data.frame(pp = 0, mu0 = 0, mu1 = 0, mu2 = 0))} # everything is zero ==> nonestimable, set pi = 0
    if (min(sum.x.y) ==0 ) { return(data.frame(pp = 0, mu0 = 0, mu1 = sum.x.y[1]/len, mu2 = sum.x.y[2]/len))} # everything is zero ==> nonestimable, set pi = 0
    # EM algorithm
    iter = 0
    if (showFlag) {print(c("iter", "pi", paste0("mu",0:2), "loglik"))}
    if (showFlag) {print(c(iter, round(param,5), lik.BZIP(xvec,yvec, param)))}
    repeat {
      iter = iter + 1
      # print(lik(vec, pp=param[1], m0=param[2], m1=param[3], m2=param[4])) # debug
      param.old <- param # saving old parameters
      
      #if (min(sum.x.y) == 0 ) {result <- data.frame(par = 0)} else {
      #  result <- optim(par = param[2], lik.BZIP, x = vec[,1], y = vec[,2], pp = param[1], xi1 = xi[1], xi2 = xi[2], lower = 0, upper = min(xi), method="Brent")
      #}
      
      xi <- xi.k(param)
      param[1] <- len.0/len*xi
      mu.x.y <- sum.x.y / len / (1 - param[1])
      param[2] <- optim(param[2], fun.c.profile, lower = 0, upper = min(mu.x.y), vec = vec.s, pp = param[1], mu.x.y = mu.x.y, method="Brent")$par
      param[3:4] <- mu.x.y - param[2]
      
      # if (len == len.s | param[1] < 0) { param[1] <- 0}
      if (showFlag) {print(c(iter, round(param,5), lik.BZIP(xvec, yvec, param)))}
      if (max(abs(param - param.old)) <= tol) {
        param <- data.frame(matrix(param,1,4))
        names(param) <- c("pp", "mu0", "mu1", "mu2")
        return(param)
        break
      } else if (iter >= 1000) {
        print("iteration exceeded")
        param <- data.frame(matrix(NA,1,4))
        names(param) <- c("pp", "mu0", "mu1", "mu2")
        return(param)
        break
      }  
    }
  }
  #if (len.s == 1) { # when there is only one nonzero pair, pi should be 0
  #  result.a <- (ML.BP(xvec=xvec, yvec=yvec)[,1:3])
  #  return(data.frame(pp= 0, mu0 = result.a$mu0, mu1 = result.a$mu1, mu2 = result.a$mu2))
  #} else {
  if (is.null(initial)) { # when initial(starting point) is not provided
    initial <- rep(0,4)
    initial[1] <- 1 - len.s/len # inital guess for pi
    initial[3:4] <- apply(vec.s, 2, sum)/len.s
    initial[2] <- min(initial[3:4])
    initial[3:4] <- initial[3:4] - initial[2]
  }
  if (initial[1] == 0) {
    initial[1] <- 0.01
  }
  result <- fun.d(param = initial, vec.s = vec.s, len = len, len.s = len.s)
  return(result)
}
if (FALSE) {
  ML.BZIP(tmp.1[,1],tmp.1[,2],showFlag=T) # 0.085183 0.3744193 0.04096429 0.08468888
  lik.BZIP(tmp.1[,1],tmp.1[,2], pp=.08511, m0 = .3744, m1 = .04096, m2=.0847) #lik = -55
  lik.BZIP(tmp.1[,1],tmp.1[,2], pp=.5, m0 = 1, m1 = .1, m2=.1) #lik = -59
  # far away from truth
  
}



if (FALSE) {# example
  ML.BZIP(vec[,1], vec[,2], showFlag=TRUE, initial = rep(0,4), tol = 1e-10)
  ML.BZIP(vec[,1], vec[,2], showFlag=TRUE, tol = 1e-8)
  ML.BZIP(extractor(1), extractor(2), initial = rep(0,4), showFlag=TRUE)
  ML.BZIP(rep(0,100), rep(0,100), initial = rep(0,4), showFlag=TRUE)
  # For all 0 pairs, identifiability issue: same likelihood for (pi=1) and (pi=0)
  ML.BZIP(c(1,rep(0,99)), rep(0,100), initial = rep(0,4), showFlag=TRUE)
  ML.BZIP(c(1,rep(0,799)), rep(0,800), initial = rep(0,4), showFlag=TRUE)
  ML.BZIP(c(1,rep(0,799)), rep(0,800), initial = c(799/800,0,0,0), showFlag=TRUE) # not converging. initial should be all zero
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

#### simulation to check if ML works okay
set.seed(1); tt(1)
result <- data.frame(pp= NA, mu0 = NA, mu1 = NA, mu2 = NA)
for (i in 1:30) {
  print (c(i, "/ 30"))
  tmp.1 <- rBZIP(50, .5, 1, 1, 1)
  result[i,] <- ML.BZIP(tmp.1[,1], tmp.1[,2])
  print(result[i,])
}
result; tt(2)
apply(result,2,mean) # 0.507, 0.865, 1.148, 1.169

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

entropy.BZIP <- function(pp, m0, m1, m2, summation = 100) {
  if (summation < 2* max(m0 + m1, m0 + m2)) {
    summation <- 2* max(m0 + m1, m0 + m2)
  }
  joint <- sapply(0:summation, function(x) {sapply(0:summation, function(y) {dBP(x = x, y = y, m0 = m0, m1 = m1, m2 = m2)})})
  joint <- joint *(1-pp)
  joint[1,1] <- joint[1,1] + pp
  marginal1 <- apply(joint,1,sum)
  marginal1 <- apply(joint,2,sum)
  
  H1 <- - marginal1 * log(marginal1); H1[!is.finite(H1)] <- 0; H1 <- sum(H1)
  H2 <- - marginal2 * log(marginal2); H2[!is.finite(H2)] <- 0; H2 <- sum(H2)
  H12 <- - joint * log(joint);        H12[!is.finite(H12)] <- 0; H12 <- sum(H12)
  M12 <- max(- H12 + H1 + H2,0)
  return(data.frame(H1 = H1, H2 = H2, H12 = H12, MI = M12))
}  
entropy.BZIP.vec <- Vectorize(entropy.BZIP)


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

