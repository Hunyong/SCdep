
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
ML.BZIP.old4 <- function(xvec, yvec, tol = 1e-6, initial = NULL, showFlag = FALSE, counter=0) { #MLE based on score equations : fail (not convex)
  # counter for recursive function
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
      } else if (iter >= 100) {
        if (counter + iter > 1000) {
          print("iteration exceeded")
          param <- data.frame(matrix(NA,1,4))
          names(param) <- c("pp", "mu0", "mu1", "mu2")
          return(param)  
          break
        }
        param[1] <- param[1]*.1
        return(ML.BZIP(xvec=xvec, yvec=yvec, tol = max(1e-8,tol*.1), initial = param, showFlag = showFlag, counter= counter + iter))
      }  
    }
  }
  #if (len.s == 1) { # when there is only one nonzero pair, pi should be 0
  #  result.a <- (ML.BP(xvec=xvec, yvec=yvec)[,1:3])
  #  return(data.frame(pp= 0, mu0 = result.a$mu0, mu1 = result.a$mu1, mu2 = result.a$mu2))
  #} else {
  if (is.null(initial)) { # when initial(starting point) is not provided
    initial <- rep(0,4)
    initial[1] <- (1 - len.s/len - exp(-sum(sum.x.y)/len.s))/(1 - exp(-sum(sum.x.y)/len.s)) # inital guess for pi
    initial[3:4] <- apply(vec.s, 2, sum)/len.s
    initial[2] <- min(initial[3:4])
    initial[3:4] <- initial[3:4] - initial[2]
  }
  if (initial[1] == 0) {
    initial[1] <- 10^-10
  }
  result <- fun.d(param = initial, vec.s = vec.s, len = len, len.s = len.s)
  return(result)
}

ML.BZIP.old5 <- function(xvec, yvec, tol = 1e-6, initial = NULL, showFlag = FALSE, counter=0, likelihood = -Inf, param.prev = rep(0,4)) { #MLE based on score equations : fail (not convex)
  # counter, likelihood, param.prev for recursive function
  len <- length(xvec)          # n
  vec <- data.frame(xvec=xvec, yvec=yvec)
  vec.s <- vec[apply(vec,1,sum) != 0,]  # nonzero pair (at least one is nonzero)
  len.s <- dim(vec.s)[1]       # n* (number of nonzero pair)
  len.0 <- len - len.s         # n0
  sum.x.y <- apply(vec,2,sum)  # mu0+mu1, mu0+mu2
  if (sum(sum.x.y) ==0 ) { return(data.frame(pp = 0, mu0 = 0, mu1 = 0, mu2 = 0))} # everything is zero ==> nonestimable, set pi = 0
  
  fun.c.profile <- function(vec, mu0, pp, mu.x.y) {    # solve for mu0
    return(-lik.BZIP(x = vec[,1], y= vec[,2], pp=pp, m0 = mu0, m1 = mu.x.y[1] - mu0, m2 = mu.x.y[2] -mu0))
  }
  xi.k <- function(param) {   #update xi = E(Delta_i | obs'd, kth param)
    return(param[1] / (param[1] + (1-param[1]) *exp(- param %*% c(0,1,1,1) )))
  }
  param.update <- function(param, vec, pp = NULL) {
    xi <- xi.k(param)
    if (is.null(pp)) {param[1] <- len.0/len*xi} else {
      param[1] <- pp
    }
    mu.x.y <- sum.x.y / len / (1 - param[1])
    if (min(mu.x.y) == 0) { param[2] <- 0} else {
      param[2] <- optim(param[2], fun.c.profile, lower = 0, upper = min(mu.x.y), vec = vec, pp = param[1], mu.x.y = mu.x.y, method="Brent")$par
    }
    param[3:4] <- mu.x.y - param[2]
    return(param)
  } 
  fun.d <- function(param, vec, vec.s, len, len.s = dim(vec.s)[1]) {
    # EM algorithm
    iter = 0
    likelihood <- lik.BZIP(xvec,yvec, param)
    if (showFlag) {print(c("iter", "pi", paste0("mu",0:2), "loglik"))}
    if (showFlag) {print(c(iter+counter, round(param,5), likelihood))}
    repeat {
      iter = iter + 1
      # print(lik(vec, pp=param[1], m0=param[2], m1=param[3], m2=param[4])) # debug
      param.old <- param # saving old parameters
      likelihood.old <- likelihood
      
      # for fast computing use nonzero cells only "vec.s" (resulting likelihood is not guaranteed to exceed previous lik)
      param <- param.update (param = param, vec = vec.s)
      likelihood <- lik.BZIP(xvec, yvec, param)
      
      if (likelihood < likelihood.old) { # likelihood check. if not better, then use whole vector "vec"
        param <- param.old
        param <- param.update (param = param, vec = vec)
        likelihood <- lik.BZIP(xvec, yvec, param)
      }
      
      # if (len == len.s | param[1] < 0) { param[1] <- 0}
      if (showFlag) {print(c(iter+counter, round(param,5), lik.BZIP(xvec, yvec, param)))}
      if (max(abs(param - param.old)) <= tol) {
        param <- data.frame(matrix(param,1,4))
        names(param) <- c("pp", "mu0", "mu1", "mu2")
        return(param)
        break
      } else if (iter >= 15) {
        if (counter + iter > 1000) {
          print("iteration exceeded")
          param <- data.frame(matrix(NA,1,4))
          names(param) <- c("pp", "mu0", "mu1", "mu2")
          return(param)  
          break
        }
        param.jump <- param.update (param = param, vec = vec, pp=param[1]*.1)
        param <- ML.BZIP(xvec=xvec, yvec=yvec, tol = max(1e-8,tol*.1), initial = param.jump, showFlag = showFlag, counter= counter + iter, likelihood = likelihood, param.prev = param)
        return(param)
      }  
    }
  }
  #if (len.s == 1) { # when there is only one nonzero pair, pi should be 0
  #  result.a <- (ML.BP(xvec=xvec, yvec=yvec)[,1:3])
  #  return(data.frame(pp= 0, mu0 = result.a$mu0, mu1 = result.a$mu1, mu2 = result.a$mu2))
  #} else {
  if (is.null(initial)) { # when initial(starting point) is not provided
    initial <- rep(0,4)
    initial[1] <- (1 - len.s/len - exp(-sum(sum.x.y)/len.s))/(1 - exp(-sum(sum.x.y)/len.s)) # inital guess for pi
    initial[3:4] <- apply(vec.s, 2, sum)/len.s
    initial[2] <- min(initial[3:4])
    initial[3:4] <- initial[3:4] - initial[2]
  }
  if (initial[1] == 0) {
    initial[1] <- 10^-10
  }
  if (lik.BZIP(xvec,yvec, initial) < likelihood) {
    initial <- param.prev
    a <- sapply(c(1,2,2.7,3), function(k) {
      initial.candidate <- param.update(param = initial, vec=vec, pp = initial[1]*k/3)
      likelihood.candidate <- lik.BZIP(xvec,yvec, initial.candidate)
      return(c(likelihood.candidate, initial.candidate))
    })
    initial <- a[2:5, which(a[1,]==max(a[1,]))[1]]
  }
  result <- fun.d(param = initial, vec = vec, vec.s=vec.s, len = len, len.s = len.s)
  return(result)
}


ML.BZIP.B2.old1 <- function(xvec, yvec, tol = 1e-6, initial = NULL, showFlag = FALSE) { #MLE based on score equations : fail (not convex)
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
    
    result <- c(EE, EE1U, EE2U, EE3U)
    
    m0.den <- sum(EE[1:3])
    m0.num <- (EE1U + EE2U + EE3U)
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
    print(initial)
  }
  
  booster <- function (param.matrix, xvec, yvec, n.cand = 10) {
    param.matrix[,1:4] <- qlogis(param.matrix[,1:4])  # logit transformation for probs
    a <- param.matrix[1,]
    b <- param.matrix[10,]
    candidate <- matrix(b, byrow=TRUE, ncol=7, nrow = n.cand)
    index <- which((abs(b-a) > 1e-6) & is.finite(b) & is.finite(a))  # target param for grid search
    for (s in 1:n.cand) {
      candidate[s,index] <- b[index] + (b[index] - a[index]) * 2^(s-1)
    }
    candidate[,1:4] <- plogis(candidate[,1:4])  # back-transformation for probs
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
    iter = iter + 1
    # print(lik(vec, pp=param[1:4], m0=param[5], m1=param[6], m2=param[7])) # debug
    param.old <- param # saving old parameters
    param <- param.update (x = xvec, y = yvec, p1 = param[1], p2 = param[2], p3 = param[3], p4 = param[4], m0=param[5], m1=param[6], m2=param[7])
    # boosting
    if (iter == 21 + boost*30) {  # Creating an empty matrix
      param.boost <- matrix(NA, nrow = 10, ncol = 7)
      boost <- boost + 1
    }
    if (iter >= 21 + (boost-1)*30 & iter <= 30 + (boost-1)*30 ) {  # Storing last ten params
      param.boost[iter - (20 + (boost-1)*30),] <- param
    }
    if (iter == 30 + (boost-1)*30) {
      param.boost <- booster(param.boost, xvec, yvec, n.cand = ifelse(index > 8, 20, 10))
      index <- which.max(param.boost[,8])
      print(param.boost)
      print(paste0("Jump to the ",index, "th parameter"))
      param <- param.boost[index,1:7]
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
