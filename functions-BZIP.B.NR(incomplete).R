ML.BZIP.B.NR <- function(xvec, yvec, tol = 1e-6, initial = NULL, showFlag = FALSE) { #MLE based on score equations : fail (not convex)
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
  fun.cond.exp <- Vectorize(fun.cond.exp)
  
  # M-step 
  param.update <- function(x, y, p1, p2, p3, p4, m0, m1, m2) {
    result <- fun.cond.exp(x, y, p1, p2, p3, p4, m0, m1, m2)
    #print(result)  ####debug
    result[result < 0] <- 0
    result <- apply(result,1,mean)
    Score <- vector("numeric",6)
    Score[1] <- (result[1] / p1  - (1-sum(result[1:3])) / p4)  #average score vector
    Score[2] <- (result[2] / p2  - (1-sum(result[1:3])) / p4)
    Score[3] <- (result[3] / p3  - (1-sum(result[1:3])) / p4)
    Score[4] <- (result[6] / m0  - result[5])
    Score[5] <- (result[8] / m1  - result[7])
    Score[6] <- (result[10] / m2  - result[9])
    Info <- matrix(0,6,6)
    diag(Info) <- c(result[1:3] / c(p1,p2,p3)^2 , result[c(6,8,10)] / c(m0,m1,m2)) #average Info
    Info[1,c(2,3)] <- Info[c(2,3),1] <- Info[2,3] <- Info[3,2] <- - (1-sum(result[1:3])) / p4^2
    #print(Score)
    #print(Info)
    result2 <- solve(Info) %*% Score
    result2 <- c(result[1:3], 1 - sum(result[1:3]), result[4:6])
    result2 <- as.numeric(c(p1, p2, p3, p4, m0, m1, m2)) -  result2
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
