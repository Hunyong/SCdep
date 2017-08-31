# BvZINB4: BvZINB3 + varying zero inflation parameters
library(distr) # for igamma (inverse digamma)

dBvZINB4 <- function(x, y, a0, a1, a2, b1, b2, p1, p2, p3, p4, log=FALSE) {
  dxy <- dBvNB3(x=x, y=y, a0=a0, a1=a1, a2=a2, b1=b1, b2=b2, log=FALSE)
  dx <- dnbinom(x=x, a0+a1, 1/(1+b1))
  dy <- dnbinom(x=y, a0+a1, 1/(1+b2))
  result <- dxy * p1 + dx * ifelse(y==0,p2,0) + dy * ifelse(x==0,p3,0) + ifelse(x+y==0,p4,0)
  return(ifelse(log, log(result), result))
}
if (FALSE) {
  dBvZINB4(1,1,1,1,1,1,.5,.25,.25,.25,.25)
  tmp <- sapply(0:50, function(r) sapply(0:50, function(s) dBvZINB4(s,r,1,1,1,1,.5,.25,.25,.25,.25)))
  sum(tmp)
}

dBvZINB4.vec <- Vectorize(dBvZINB4)
lik.BvZINB4 <- function(x, y, param) {
  sum(log(dBvZINB4.vec(x, y, param[1], param[2], param[3], param[4], param[5], param[6], param[7], param[8], param[9])))
}

!! ## modify code so that the likelihood doesn't exceed 1.(Some very small number is calculated as a large positive number)

### 2.EM

### nonzero cells: (1-pp) was not multiplied by!!! this caused decreasing likelihood in EM
dBvZINB4.Expt <- function(x, y, a0, a1, a2, b1, b2, p1, p2, p3, p4) {
  # Base density
    t1 = (b1 + b2 + 1) /(b1 + 1); t2 = (b1 + b2 + 1) /(b2 + 1)
    l1 <- function(k, m) exp(lgamma(a1 + k) - lgamma(k+1) - lgamma(a1) + lgamma(x + y + a0 -m -k) - lgamma(x -k +1) - lgamma(a0 + y - m) + k *log(t1)
                             + lgamma(m + a2) - lgamma(m+1) - lgamma(a2) + lgamma(y +a0 -m) - lgamma(y -m +1) - lgamma(a0) + m *log(t2))
    l1.B <- - (+x+y+a0)*log(1 + b1 + b2) + x * log(b1) + y * log(b2) - a1 * log(1 + b1) - a2 * log(1 + b2)
    l1.B <- exp(l1.B) * p1
    l2.B <- exp(- (x + a0 + a1)*log(1 + b1) + x * log(b1)) * p2 * ifelse(y==0, 1, 0)
    l3.B <- exp(- (y + a0 + a2)*log(1 + b2) + y * log(b2)) * p3 * ifelse(x==0, 1, 0)
    l4.B <- p4 * ifelse(x + y == 0, 1, 0)
    
    l.mat <- sapply(0:x, function(k) sapply(0:y, l1, k=k))  # %>% print
    l.sum <- sum(l.mat * (l1.B + l2.B +  l3.B +  l4.B))     # %>% print

  # expectation components
    R0.E1 <- function(k, m) {x - k + y - m + a0}
    log.R0.E1 <- function(k, m) {digamma(x - k + y - m + a0)}
    R0.E1.B <- b1/(1 + b1 + b2)
    R0.E2.B <- b1/(1 + b1)
    R0.E3.B <- b1/(1 + b2)
    R0.E4.B <- b1
    
    R1.E1 <- function(k) {k + a1}
    log.R1.E1 <- function(k) {digamma(k + a1)}
    R1.E1.B <- b1/(1 + b1)
    R1.E2.B <- b1/(1 + b1)
    R1.E3.B <- b1
    R1.E4.B <- b1
    
    R2.E1 <- function(m) {m + a2}
    log.R2.E1 <- function(m) {digamma(m + a2)}
    R2.E1.B <- b1/(1 + b2)
    R2.E2.B <- b1
    R2.E3.B <- b1/(1 + b2)
    R2.E4.B <- b1
    
    R0.mat <- sapply(0:x, function(k) sapply(0:y, R0.E1, k=k))
    R0.mat <- R0.mat * l.mat
    R0.E <- sum(R0.mat*(l1.B * R0.E1.B + l2.B * R0.E2.B 
                         +  l3.B * R0.E3.B +  l4.B * R0.E4.B))/l.sum
    
    R1.mat <- t(matrix(sapply(0:x, R1.E1), x+1, y+1))
    R1.mat <- R1.mat * l.mat
    R1.E <- sum(R1.mat*(l1.B * R1.E1.B + l2.B * R1.E2.B 
                        +  l3.B * R1.E3.B +  l4.B * R1.E4.B))/l.sum
    
    R2.mat <- matrix(sapply(0:y, R2.E1), y+1, x+1) #%>% print
    R2.mat <- R2.mat * l.mat
    R2.E <- sum(R2.mat*(l1.B * R2.E1.B + l2.B * R2.E2.B 
                        +  l3.B * R2.E3.B +  l4.B * R2.E4.B))/l.sum
    
    log.R0.mat <- sapply(0:x, function(k) sapply(0:y, log.R0.E1, k=k))
    log.R0.mat <- log.R0.mat * l.mat
    log.R0.E <- sum(log.R0.mat*(l1.B + l2.B +  l3.B +  l4.B))/l.sum
    log.R0.E <- log.R0.E +
      sum (c(l1.B, l2.B, l3.B, l4.B) * log(c(R0.E1.B, R0.E2.B, R0.E3.B, R0.E4.B))) / sum (l1.B, l2.B, l3.B, l4.B)
    
    log.R1.mat <- t(matrix(sapply(0:x, log.R1.E1), x+1, y+1))
    log.R1.mat <- log.R1.mat * l.mat
    log.R1.E <- sum(log.R1.mat*(l1.B + l2.B +  l3.B +  l4.B))/l.sum
    log.R1.E <- log.R1.E +
      sum (c(l1.B, l2.B, l3.B, l4.B) * log(c(R1.E1.B, R1.E2.B, R1.E3.B, R1.E4.B))) / sum (l1.B, l2.B, l3.B, l4.B)
    
    
    log.R2.mat <- matrix(sapply(0:y, log.R2.E1), y+1, x+1)
    log.R2.mat <- log.R2.mat * l.mat
    log.R2.E <- sum(log.R2.mat*(l1.B + l2.B +  l3.B +  l4.B))/l.sum
    log.R2.E <- log.R2.E +
      sum (c(l1.B, l2.B, l3.B, l4.B) * log(c(R2.E1.B, R2.E2.B, R2.E3.B, R2.E4.B))) / sum (l1.B, l2.B, l3.B, l4.B)
    
    E.E <- c(l1.B, l2.B, l3.B, l4.B)
    E.E <- E.E/sum(E.E)
    
    v.E <- y * ifelse(y == 0, 0, 1) + (a0 + a2) * b2 * sum(E.E[3:4])
    
    result <- c(log(l.sum), R0.E, R1.E, R2.E, log.R0.E, log.R1.E, log.R2.E, E.E, v.E)
    names(result) <- c("logdensity", paste0("R", 0:2, ".E"), paste0("log.R", 0:2, ".E"), paste0("E",1:4,".E"), "v.E")
    return(result)
}
dBvZINB4.Expt.vec <- Vectorize(dBvZINB4.Expt)

if (FALSE) {
  tmp <- dBvZINB4.Expt.vec(c(1,1,1),c(0,1,2),1,1,1,1,2,.25,.25,.25,.25)
  tmp <- dBvZINB4.Expt.vec(c(0,1,1),c(0,1,2),1,1,1,1,2,.25,.25,.25,.25)
  tmp <- dBvZINB4.Expt.vec(extractor(1),extractor(2),1,1,1,1,2,.25,.25,.25,.25)
  t(tmp)[21:40,]
  dBvZINB4.Expt.vec(c(10,1,2),c(10,1,1), 1.193013282, 0.003336139, 0.002745513, 3.618842924, 3.341625901, .25,.25,.25,.25)
}

ML.BvZINB4 <- function (xvec, yvec, initial = NULL, tol=1e-8, showFlag=FALSE) {
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
    initial <- rep(NA,9)
    initial[4] <- s2.x /xbar
    initial[5] <- s2.y /ybar
    initial[2:3] <- c(xbar,ybar)/initial[4:5]
    initial[1] <- min(initial[2:3]) * abs(cor.xy)
    initial[2:3] <-  initial[2:3] - initial[1]
    
    initial[6:9] <- bin.profile(xvec, yvec)   # freq of each zero-nonzero profile
    initial[6:9] <- initial[6:9]/sum(initial[6:9])      # relative freq
    initial <- pmax(initial, 1e-5)
    
    print(initial) ###
  }
  
  
  iter = 0
  param = initial
  if (showFlag) {print(c("iter", "a0", "a1", "a2", "b1", "b2", paste0("p",1:4), "lik"))}
  repeat {
    iter = iter + 1
    #print(c(param))
    #print(lik(vec, pp=param[1], m0=param[2], m1=param[3], m2=param[4])) # debug
    param.old <- param # saving old parameters
    
    # updating
    expt <- dBvZINB4.Expt.vec(xy.reduced$x, xy.reduced$y, 
                              a0 = param[1], a1 = param[2], a2 = param[3], b1 = param[4], b2 = param[5], p1 = param[6], p2 = param[7], p3 = param[8], p4 = param[9])
    expt <- as.vector(expt %*% xy.reduced$freq / n) #%>% print
    # loglik = expt[1] * n
    delta <- expt[12] / (expt[2] + expt[4])                   # delta = E(V) / (E(xi0 + xi2))
    param[6:9] = expt[8:11]                                                    # pi = E(Z)
    
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
      initial = rep(1,4)
      result <- multiroot(opt.vec, start=initial[1:4], rtol=1e-20)$root
    }
    param[1:4] <- exp(result)
    param[5]   <- param[4] * delta                                       # b2
    param[6:9] <- expt[8:11]
    
    #print (expt) #####
    if (showFlag) {print(c(iter, round(param,4), expt[1] * n))} #lik: lik of previous iteration
    if (max(abs(param - param.old)) <= tol) {
      names(param) <- c("a0", "a1", "a2", "b1", "b2", paste0("p",1:4))
      return(param)
      break
    }
  }
  #result <- data.frame(a0 = param[1], a1 = param[2], a2 = param[3], b1 = param[4], b2 = param[5], pi = param[6])
  #return(result)
}
# simple tests
if (FALSE) {
  ML.BvZINB4(c(10,1,1),c(10,1,2), showFlag=TRUE) # c(1.193014333 0.002745514 0.003336045 3.341621165 3.618839217 0.000000000 )
  # lik.BvZINB4(c(10,1,1),c(10,1,2),c(0.7186211, 0.4954254, 0.5652637, 2.9788157, 3.0834235, 1,0,0,0))
  # [1] -13.82585
  # lik.BvZINB4(c(10,1,1),c(10,1,2),c(1.193014333, 0.002745514, 0.003336045, 3.341621165, 3.618839217, 1,0,0,0))
  # [1] -12.90997
  ML.BvZINB4(c(0,1,1),c(0,1,5), showFlag=TRUE)
  tt(1)
  ML.BvZINB4(extractor(1), extractor(4), showFlag=TRUE)
  ML.BvZINB4(extractor(5), extractor(6), showFlag=TRUE)
  tt(2)
  ML.BvZINB4(extractor(1), extractor(3), showFlag=TRUE)
  # 0.000799916  0.015057420  0.006208375 67.414607790  9.180617081  0.361266622
  
  lik.BvZINB3(extractor(1), extractor(4),c(0.0004349035, 0.009488825, 0.003788559, 68.25597, 9.835188, .95)) # -391.5657
  lik.BvNB3(extractor(1), extractor(4),c(0.0004349035, 0.009488825, 0.003788559, 68.25597, 9.835188)) # -308.7
  
  
  3.404230e-05 9.740676e-03 5.435834e-03 7.059027e+01 6.627206e+00 9.500000e-01
  
  
  tt(2) # 31secs
  
  tt(1)
  ML.BvZINB4(extractor(1), extractor(3),showFlag=TRUE)
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
      candidate[s,index] <- b[index] + (b[index] - a[index]) * 3^(s-1)
    }
    candidate[,6] <- plogis(candidate[,6])  # back-transformation
    candidate[,1:5] <- exp(candidate[,1:5])  # back-transformation for probs
    #print(candidate[,1:4])  #debug
    
    lik <- sapply(1:n.cand, function(s) {lik.BvZINB3(xvec, yvec, candidate[s,])})
    lik <- ifelse(is.infinite(lik), -Inf, lik)  # sometimes likelihood is inf which is nonsense. force it to -Inf
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
      param.boost <- booster(param.boost, xvec, yvec, n.cand = min(max(5, index * 2),20))
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

