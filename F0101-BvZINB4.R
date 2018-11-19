# BvZINB4: BvZINB3 + varying zero inflation parameters
source("F0101-BvNB3.R")
source("F0101-BvZINB4-supp.R")

dBvZINB4 <- function(x, y, a0, a1, a2, b1, b2, p1, p2, p3, p4, log=FALSE) {
  dxy <- dBvNB3(x=x, y=y, a0=a0, a1=a1, a2=a2, b1=b1, b2=b2, log=FALSE)
  dx <- dnbinom(x=x, a0+a1, 1/(1+b1))
  dy <- dnbinom(x=y, a0+a2, 1/(1+b2))
  result <- dxy * p1 + dx * ifelse(y==0,p2,0) + dy * ifelse(x==0,p3,0) + ifelse(x+y==0,p4,0)
  return(ifelse(log, log(result), result))
}
if (FALSE) {
  dBvZINB4(1,1,1,1,1,1,.5,.25,.25,.25,.25)
  tmp <- sapply(0:50, function(r) sapply (0:50, function(s) dBvZINB4(s,r,1,1,1,1,.5,.25,.25,.25,.25)))
  sum(tmp)
}

dBvZINB4.vec <- Vectorize(dBvZINB4)
lik.BvZINB4 <- function(x, y, param) {
  sum(log(dBvZINB4.vec(x, y, param[1], param[2], param[3], param[4], param[5], param[6], param[7], param[8], param[9])))
}

rBvZINB4 <- function(n, a0, a1, a2, b1, b2, p1, p2, p3, p4, param=NULL) {
  if (!is.null(param)) {a0 = param[1]; a1 = param[2]; a2 = param[3]; b1 = param[4]; b2 = param[5]
                        p1 = param[6]; p2 = param[7]; p3 = param[8]; p4 = param[9]
  }
  
  
  rmat <- matrix(rgamma(n*3, shape = c(a0, a1, a2), rate = 1/b1), n, 3, byrow=TRUE)
  rmat2 <- rmat
  rmat2[,3] <- rmat2[,1] + rmat2[,3]
  rmat2[,2] <- rmat2[,1] + rmat2[,2]
  rmat2 <- rmat2[,2:3]
  rmat2[,2] <- rmat2[,2]*b2/b1
  uv <- matrix(rpois(n*2, rmat2), n, 2)
  
  E <- t(rmultinom(n, 1, c(p1, p2, p3, p4)))
  z <- cbind(E[,1]+E[,2], E[,1]+E[,3])
  
  xy <- uv * z
  colnames(xy) <- c("x", "y")
  
  return(xy)
}


### 2.EM

### nonzero cells: (1-pp) was not multiplied by!!! this caused decreasing likelihood in EM
dBvZINB4.Expt <- function(x, y, a0, a1, a2, b1, b2, p1, p2, p3, p4, debug = FALSE) {
  # Base density
  t1 = (b1 + b2 + 1) /(b1 + 1); t2 = (b1 + b2 + 1) /(b2 + 1)
  adj.A <- adj.B1 <- adj.C <- adj.sum <- 0
  l1 <- function(k, m, adjj=0) exp(lgamma(a1 + k) - lgamma(k+1) - lgamma(a1) + lgamma(x + y + a0 -m -k) - lgamma(x -k +1) - lgamma(a0 + y - m) 
                           + lgamma(m + a2) - lgamma(m+1) - lgamma(a2) + lgamma(y +a0 -m) - lgamma(y -m +1) - lgamma(a0) - adjj)
  l1.C <- function(k, m, adjj=0) exp(k *log(t1) + m *log(t2) - adjj)
  l1.B <- - (+x+y+a0)*log(1 + b1 + b2) + x * log(b1) + y * log(b2) - a1 * log(1 + b1) - a2 * log(1 + b2)
  # l1.B to be updated several lines later depending on l2.B ~ l4.B
  l2.B <- exp(- (x + a0 + a1)*log(1 + b1) + x * log(b1) + adj.B1) * p2 * ifelse(y==0, 1, 0)
  l3.B <- exp(- (y + a0 + a2)*log(1 + b2) + y * log(b2) + adj.B1) * p3 * ifelse(x==0, 1, 0)
  l4.B <- p4 * ifelse(x + y == 0, 1, 0) * exp(adj.B1)
  
  # l1.AC For numerical stability use only.
  l1.AC <- function(k, m, adjj=0) exp(lgamma(a1 + k) - lgamma(k+1) - lgamma(a1) + lgamma(x + y + a0 -m -k) - lgamma(x -k +1) - lgamma(a0 + y - m) 
                                   + lgamma(m + a2) - lgamma(m+1) - lgamma(a2) + lgamma(y +a0 -m) - lgamma(y -m +1) - lgamma(a0) + k *log(t1) + m *log(t2) - adjj)
  
  
  # cat("l1.B ", l1.B,"\n")  
  if (l1.B < - 200 & log(l2.B + l3.B + l4.B) < 0) {
if (debug) cat("adjustment activated for l1.B\n")     
    adj.B1 = ((-l1.B - 200) %/% 100) * 100 # prevent exp(l1.B) from being 0
    l1.B = l1.B + adj.B1
  }
  l1.B <- exp(l1.B) * p1
if (debug)  cat("l1.B ", l1.B,"\n")  
  
    
  
  l.A.mat <- sapply(0:x, function(k) sapply(0:y, l1, k = k, adjj = adj.A))  # %>% print
  l.C.mat <- sapply(0:x, function(k) sapply(0:y, l1.C, k = k, adjj = adj.C))  # %>% print
  while (log(sum( l.A.mat)) > 250) {
if (debug)  cat("adjustment activated for A.mat\n")     
    adj.A = adj.A + 200
    l.A.mat <- sapply(0:x, function(k) sapply(0:y, l1, k = k, adjj = adj.A))  # %>% print
  }
  while (log(sum( l.C.mat)) > 250) {
if (debug)  cat("adjustment activated for C.mat\n")     
    adj.C = adj.C + 200
    l.C.mat <- sapply(0:x, function(k) sapply(0:y, l1.C, k = k, adjj = adj.C))  # %>% print
  }
  
# print(l.C.mat)
  
#   if (is.infinite(sum( l.A.mat))) {
# cat("activated once")    
#     adj.A = 200
#     l.A.mat <- sapply(0:x, function(k) sapply(0:y, function(m) {l1(k =k, m = m) *exp(-adj.A)}))
#     if (is.infinite(sum( l.A.mat))) { ## added for further adjustment
# cat("activated twice")      
#       adj.A = 500
#       l.A.mat <- sapply(0:x, function(k) sapply(0:y, function(m) {l1(k =k, m = m) *exp(-adj.A)}))
#     }  
#   } #%>%print
  #adjustment is cancelled out for each Expectation, so can be ignored. But for the final likelihood it should be adjusted at the end.
  sum.AC <- sum(l.A.mat * l.C.mat)

  if (is.infinite(sum.AC)| log(sum.AC) > 200) {
    if (debug) cat("adjustment activated for AC.mat (too large)\n")
    adj.A = adj.A + 100
    adj.C = adj.C + 100
    l.A.mat <- sapply(0:x, function(k) sapply(0:y, l1, k = k, adjj = adj.A))  # %>% print
    l.C.mat <- sapply(0:x, function(k) sapply(0:y, l1.C, k = k, adjj = adj.C))  # %>% print
    sum.AC <- sum(l.A.mat * l.C.mat)
  } else if (log(sum.AC) < - 100) {
    if (debug) cat("adjustment activated for AC.mat (too small)\n")
    adj.A = adj.A - 200  # floor(log(sum(l.A.mat)/x/y)*2/3)
    adj.C = adj.C - 200
    l.A.mat <- sapply(0:x, function(k) sapply(0:y, l1, k = k, adjj = adj.A))  # %>% print
    l.C.mat <- sapply(0:x, function(k) sapply(0:y, l1.C, k = k, adjj = adj.C))  # %>% print
    l.AC.mat <- sapply(0:x, function(k) sapply(0:y, l1.AC, k = k, adjj = adj.C + adj.A))
    sum.AC <- sum(l.AC.mat)
    # abcde.1 <<- l.A.mat
    # abcde.2 <<- l.C.mat
    # abcde.3 <<- l.AC.mat
  }
  sum.A <- sum(l.A.mat)
  l.sum <- sum.AC * l1.B + sum.A * sum (l2.B +  l3.B +  l4.B) * exp(-adj.C)
  if (l.sum == 0) {
    adj.sum = -floor(log(sum.AC)*2/3 + log(l1.B)*2/3)
    if (debug) cat("adjustment activated for l.sum (adj = ", adj.sum, ")\n")
    l.sum <- sum.AC * exp(adj.sum) * l1.B + sum.A * (exp(adj.sum) * sum (l2.B +  l3.B +  l4.B)) * exp(-adj.C)
    # abcde.4 <<- c(l.sum = l.sum, sum.AC = sum.AC, l1.B = l1.B, sum.A = sum.A,  l2.B = l2.B, l3.B = l3.B, l4.B = l4.B, adj.C = adj.C)
    ## paranthesis matters. sum.A = some number, exp(adj.sum) = almost inf, sum(l2.B +  l3.B +  l4.B) = 0, ...
    # Then without paranthesis, Inf * 0 = NaN,
    # But with paranthesis, c * (large number * 0) = c * 0 = 0
  }

if (debug) {
  cat("sum.AC", sum.AC,"\n\n")
  cat("sum.A", sum.A,"\n\n")
  cat("sum(l.C.mat)", sum(l.C.mat),"\n\n")
  cat("l1.B", l1.B,"\n\n")
  cat("l2.B", l2.B,"\n\n")
  cat("l3.B", l3.B,"\n\n")
  cat("l4.B", l4.B,"\n\n")
  cat("l.sum ", l.sum, "\n")
}
  # print(c(l.sum, log(l.sum))); print(l.A.mat); print(l.C.mat); print(c(l1.B, l2.B, l3.B, l4.B, adj.A)) #####
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
  
  R0.mat <- sapply (0:x, function(k) sapply(0:y, R0.E1, k=k))
  R0.mat <- R0.mat * l.A.mat
  R0.E <- sum(R0.mat * l.C.mat * exp(adj.sum) * l1.B * R0.E1.B)  / l.sum + 
          sum(R0.mat*(  l2.B * R0.E2.B
                      + l3.B * R0.E3.B
                      + l4.B * R0.E4.B)*exp(-adj.C + adj.sum)) / l.sum
  
  # cat("R0.E ", R0.E, "\n")  
  R1.mat <- t(matrix(sapply(0:x, R1.E1), x+1, y+1))
  R1.mat <- R1.mat * l.A.mat
  R1.E <- sum(R1.mat * l.C.mat * exp(adj.sum) * l1.B * R1.E1.B) / l.sum +  
          sum(R1.mat*(  l2.B * R1.E2.B 
                      + l3.B * R1.E3.B 
                      + l4.B * R1.E4.B)*exp(-adj.C + adj.sum)) / l.sum
  # cat("R1.E ", R1.E, "\n")  

  R2.mat <- matrix(sapply(0:y, R2.E1), y+1, x+1) #%>% print
  R2.mat <- R2.mat * l.A.mat
  R2.E <- sum(R2.mat * l.C.mat * exp(adj.sum) * l1.B * R2.E1.B) / l.sum +  
          sum(R2.mat*(  l2.B * R2.E2.B 
                      + l3.B * R2.E3.B 
                      + l4.B * R2.E4.B)*exp(-adj.C + adj.sum)) / l.sum
  # cat("R2.E ", R2.E, "\n")  
  
  log.R0.mat <- sapply(0:x, function(k) sapply(0:y, log.R0.E1, k=k))
  log.R0.mat <- log.R0.mat * l.A.mat
  log.R0.E <- sum(log.R0.mat * l.C.mat) * exp(adj.sum) * l1.B + sum(log.R0.mat) * c(l2.B +  l3.B +  l4.B)*exp(-adj.C + adj.sum)
  log.R0.E <- log.R0.E +
              sum.AC * exp(adj.sum) * l1.B             *  log (R0.E1.B)                             +  
              sum.A  * c(l2.B, l3.B, l4.B) %*% log (c(R0.E2.B, R0.E3.B, R0.E4.B)) *exp(-adj.C + adj.sum)
  log.R0.E <- log.R0.E / l.sum
  
  # cat("log.R0.E ", log.R0.E, "\n")  
  log.R1.mat <- t(matrix(sapply(0:x, log.R1.E1), x+1, y+1))
  log.R1.mat <- log.R1.mat * l.A.mat
  log.R1.E <- sum(log.R1.mat * l.C.mat) * exp(adj.sum) * l1.B + sum(log.R1.mat) * c(l2.B +  l3.B +  l4.B)*exp(-adj.C + adj.sum)
  log.R1.E <- log.R1.E +
              sum.AC * exp(adj.sum) * l1.B             *  log (R1.E1.B)                             + 
              sum.A  * c(l2.B, l3.B, l4.B) %*% log (c(R1.E2.B, R1.E3.B, R1.E4.B))*exp(-adj.C + adj.sum)
  log.R1.E <- log.R1.E / l.sum
  
  # cat("log.R1.E ", log.R1.E, "\n")  
  log.R2.mat <- matrix(sapply(0:y, log.R2.E1), y+1, x+1)
  log.R2.mat <- log.R2.mat * l.A.mat
  log.R2.E <- sum(log.R2.mat * l.C.mat) * exp(adj.sum) * l1.B + sum(log.R2.mat) * c(l2.B +  l3.B +  l4.B)*exp(-adj.C + adj.sum)
  log.R2.E <- log.R2.E +
              sum.AC * exp(adj.sum) * l1.B                 *  log (R2.E1.B)                       + 
              sum.A  * c(l2.B, l3.B, l4.B) %*% log (c(R2.E2.B, R2.E3.B, R2.E4.B))*exp(-adj.C + adj.sum)
  log.R2.E <- log.R2.E / l.sum
  
  # cat("log.R2.E ", log.R2.E, "\n")  
  E.E <- c(sum.AC * exp(adj.sum) * l1.B, sum.A * c(l2.B, l3.B, l4.B)*exp(-adj.C + adj.sum))
  E.E <- E.E/sum(E.E)
  
  # cat("E.E ", E.E, "\n")  
  #v.E <- ifelse(y == 0, 0, y) + (a0 + a2) * b2 * sum(E.E[c(2,4)])
  v.E <- (sum.AC * exp(adj.sum) * l1.B * y + 
          sum.A * l2.B * a2 * b2*exp(-adj.C + adj.sum) +
          dnbinom(x, a0 + a1 + 1, b1/(1+b1)) * exp(-adj.A - adj.C + adj.sum) * a0 * b2 * p2 * ifelse(y==0, 1, 0) +
          sum.A * l3.B * y *exp(-adj.C + adj.sum) +
          sum.A * l4.B * (a0 + a2) * b2 *exp(-adj.C + adj.sum)) / l.sum

  result <- c(log(l.sum) + adj.A -adj.B1 + adj.C - adj.sum, R0.E, R1.E, R2.E, log.R0.E, log.R1.E, log.R2.E, E.E, v.E) #%>%print
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

# maxiter control added, output =param + lik + #iter
ML.BvZINB4 <- function (xvec, yvec, initial = NULL, tol=1e-8, maxiter = 200, showFlag=FALSE, showPlot = FALSE) {
  xy.reduced <- as.data.frame(table(xvec,yvec))
  names(xy.reduced) <- c("x", "y","freq")
  xy.reduced <- xy.reduced[xy.reduced$freq != 0,]
  xy.reduced$x <- as.numeric(as.character(xy.reduced$x))
  xy.reduced$y <- as.numeric(as.character(xy.reduced$y))
  xy.reduced$freq <- as.numeric(as.character(xy.reduced$freq))
  n <- sum(xy.reduced$freq)
  if (max(xvec)==0 & max(yvec)==0) {return(c(rep(1e-10,5),1,0,0,0, 0, 1, 0))} # 9 params, lik, iter, pureCor
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
    #print(initial) ###
  }

  cor.trace <<- data.frame(iter=1, pureCor=1)
  iter = 0
  param = initial
  if (showFlag) {print(c("iter", "a0", "a1", "a2", "b1", "b2", paste0("p",1:4), "lik", "pureCor"))}
  if (showPlot) {
    par(mfrow=c(2,1))
    par(mar=c(2,4,1,4))
    }
  repeat {
    iter = iter + 1
    # print(c(param))
    # print(lik(vec, pp=param[1], m0=param[2], m1=param[3], m2=param[4])) # debug
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
      result <- multiroot(opt.vec, start = initial[1:4], rtol=1e-20)$root
    }
    param[1:4] <- exp(result)
    param[5]   <- param[4] * delta                                       # b2
    
    pureCor <- stat.BvZINB4(param = param, measure = "pureCor")
    cor.trace[iter,] <<- c(iter,pureCor)
    if (showPlot & (iter %% 20 == 0)) {
      span <- min(max(iter-200+1,1),101):iter
      span2 <- max(iter-100+1,1):iter
      yspan <- c(min(0.2, min(cor.trace[span,2]-0.05)),max (max(cor.trace[span,2])+0.05,0.4))
      yspan2 <- c(min(max(cor.trace[span2,2]) - 0.001, min(cor.trace[span2,2]-0.001)),max (max(cor.trace[span2,2])+0.001,0.4))
      plot(cor.trace[span,"iter"], cor.trace[span,"pureCor"], xlab="iteration", ylab="pureCorrelation", pch=".", col="blue", ylim = yspan)
      plot(cor.trace[span2,"iter"], cor.trace[span2,"pureCor"], xlab="iteration", ylab="pureCorrelation", pch=20, col="red")
    }
    #print (expt) #####
    if (showFlag) {print(c(iter, round(param,4), expt[1] * n, pureCor))} #lik: lik of previous iteration
    if (maxiter <= iter) {
      lik <- lik.BvZINB4(xvec, yvec, param = param)
      result <- c(param, lik, iter, pureCor)
      names(result) <- c("a0", "a1", "a2", "b1", "b2", paste0("p",1:4), "lik","iter", "pureCor")
      return(result)
    }
    if (max(abs(param - param.old)) <= tol) {
      lik <- lik.BvZINB4(xvec, yvec, param = param)
      result <- c(param, lik, iter, pureCor)
      names(result) <- c("a0", "a1", "a2", "b1", "b2", paste0("p",1:4), "lik","iter", "pureCor")
      return(result)
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
  #1,8 -1522.5424        -483.66650 2077.7510
  
  #7,17 BvNB2 2352 -> 860
  # 1 13 1008.3853 -> not much
  #  8 53 1720.4281 -> not much
  # 8 36  2733.8509 ->  not much < 2670?
  # 6 38       -1632.6652        -544.85949 2175.6113
  # 4 44  3977.0302  -> 1200
  # 3 23       -1581.2832        -481.11177 2200.3428
  # 2 58  3660.0371  -> 1200
  # 5 36  1475.9486  -> not much
  # 5 38  2399.4086  -> not much
  #  9 28  3055.5890  -> 1060
  #  10 16        -808.4026        -226.31388 1164.1775  not much
  # 11 18   -1859.2748        -588.25860 2542.0324 not much
  # 17 18       -1729.3567        -546.02464 2366.6640 not much
  
  
  #3.404230e-05 9.740676e-03 5.435834e-03 7.059027e+01 6.627206e+00 9.500000e-01
  
  
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
# maxiter control added, output =param + lik + #iter
# Mar 15, 2018: Print pureCor instead of cor
ML.BvZINB4.2 <- function (xvec, yvec, initial = NULL, tol=1e-8, maxiter=200, showFlag=FALSE, showPlot=FALSE, cor.conv = FALSE, boosting=TRUE, debug = FALSE) {
  if (debug) {showFlag=TRUE}
  xy.reduced <- as.data.frame(table(xvec,yvec))
  names(xy.reduced) <- c("x", "y","freq")
  xy.reduced <- xy.reduced[xy.reduced$freq != 0,]
  xy.reduced$x <- as.numeric(as.character(xy.reduced$x))
  xy.reduced$y <- as.numeric(as.character(xy.reduced$y))
  xy.reduced$freq <- as.numeric(as.character(xy.reduced$freq))
  n <- sum(xy.reduced$freq)
  if (max(xvec)==0 & max(yvec)==0) {return(c(rep(1e-10,5),1,0,0,0, 0, 1, 0))} # 9 params, lik, iter, pureCor
  #print(xy.reduced)
  
  # initial guess
  if (is.null(initial)) {
    xbar <- mean(xvec); ybar <- mean(yvec); xybar <- mean(c(xbar, ybar))
    s2.x <- var(xvec); s2.y <- var(yvec); if(is.na(s2.x)|is.na(s2.y)) {s2.x <- s2.y <- 1}
    cor.xy <- cor(xvec,yvec); if (is.na(cor.xy)) {cor.xy <- 0}
    zero <- sum(xvec == 0 & yvec == 0) / n
    
    initial <- rep(NA,9)
    initial[4] <- s2.x /ifelse(xbar==0,1e-4, xbar) #%>% print
    initial[5] <- s2.y /ifelse(ybar==0,1e-4, ybar) #%>% print
    initial[2:3] <- c(xbar,ybar)/pmax(initial[4:5], c(0.1,0.1)) #%>% print
    initial[1] <- min(initial[2:3]) * abs(cor.xy) #%>% print
    initial[2:3] <-  initial[2:3] - initial[1] #%>% print
    
    initial[6:9] <- bin.profile(xvec, yvec)   # freq of each zero-nonzero profile
    initial[6:9] <- initial[6:9]/sum(initial[6:9])      # relative freq
    initial <- pmax(initial, 1e-5)
    if(is.na(sum(initial))) { initial[is.na(initial)] <- 1}
  # print(initial) ###
  }
# print(initial)
  booster <- function (param.matrix, xvec, yvec, n.cand = 10) {
    param.matrix[,6:9] <- qlogis(param.matrix[,6:9])  # logit transformation for probs
    param.matrix[,1:5] <- log(param.matrix[,1:5])  # log transformation for positives
    a <- param.matrix[1,]
    b <- param.matrix[5,]
    candidate <- matrix(b, byrow=TRUE, ncol=9, nrow = n.cand)
    index <- which((abs(b-a) > 1e-5) & is.finite(b) & is.finite(a))  # target param for grid search
    
    for (s in 1:n.cand) {
      candidate[s,index] <- b[index] + (b[index] - a[index]) * 3^(s-1)
    }
    candidate[,6:9] <- plogis(candidate[,6:9])  # back-transformation
    candidate[,6:9] <- candidate[,6:9]/ apply(candidate[,6:9],1,sum) # normalize
    candidate[,1:5] <- exp(candidate[,1:5])  # back-transformation for probs
    #print(candidate[,1:4])  #debug
    
    lik <- sapply(1:n.cand, function(s) {lik.BvZINB4(xvec, yvec, candidate[s,])})
    lik <- ifelse(is.infinite(lik), -Inf, lik)  # sometimes likelihood is inf which is nonsense. force it to -Inf
    if (sum(!is.finite(lik)) > 0) {
      return(cbind(candidate,lik)[1:max(min(which(!is.finite(lik)))-1,1),])
    } else {return(cbind(candidate,lik))}
  }
  
  cor.trace <<- data.frame(iter=1, pureCor=1)
  iter = 0
  param = initial
  lik = Inf
  pureCor = 0
  boost = 0
  index = 1 # previous boosting index
  if (showPlot) {
    par(mfrow=c(2,1))
    par(mar=c(2,4,1,4))
  }
# cat(442)
  repeat {
    iter = iter + 1
    param.old <- param # saving old parameters
abcd.old <<- param.old
    if (debug) {lik.old <- lik} #debugging
    pureCor.old <- pureCor
    # updating
# cat(449)
    expt <- dBvZINB4.Expt.vec(xy.reduced$x, xy.reduced$y, 
                              a0 = param[1], a1 = param[2], a2 = param[3], b1 = param[4], b2 = param[5], p1 = param[6], p2 = param[7], p3 = param[8], p4 = param[9])
abc <<- expt
    expt <- as.vector(expt %*% xy.reduced$freq / n)
# cat(453)
    # loglik = expt[1] * n
    delta <- expt[12] / (expt[2] + expt[4])                   # delta = E(V) / (E(xi0 + xi2))
    param[6:9] = expt[8:11]                                                    # pi = E(Z)
# cat(457)
abcd <<- param
    opt.vec <- function(par.ab) {
      par.ab <- exp(par.ab)
      r1 <- sum(expt[2:4]) - sum(par.ab[1:3]) * par.ab[4]
      r2 <- expt[5:7] - digamma(par.ab[1:3]) - log(par.ab[4])
      # print(c(r1,r2)) ###
      return(c(r1,r2))
    }
    
    param.l <- log(param)
#expt %>% print
#param.l %>% print

    result <- try(multiroot(opt.vec, start=param.l[1:4])$root, silent=TRUE)
    if (class(result)=="try-error") {
      initial = rep(1,4)
      result <- multiroot(opt.vec, start = initial[1:4], rtol=1e-20)$root
    }
    param[1:4] <- exp(result)
    param[5]   <- param[4] * delta                                       # b2
    pureCor <- stat.BvZINB4(param = param, measure = "pureCor")
    if (debug) {
      lik <- lik.BvZINB4(xvec, yvec, param = param)  #debugging
      if (lik < lik.old) warnings("likelihood decreased!")
    }
    cor.trace[iter,] <<- c(iter,pureCor)
    if (showPlot & (iter %% 20 == 0)) {
      span <- min(max(iter-200+1,1),101):iter
      span2 <- max(iter-100+1,1):iter
      yspan <- c(min(0.2, min(cor.trace[span,2]-0.05)),max (max(cor.trace[span,2])+0.05,0.4))
      yspan2 <- c(min(max(cor.trace[span2,2]) - 0.001, min(cor.trace[span2,2]-0.001)),max (max(cor.trace[span2,2])+0.001,0.4))
      plot(cor.trace[span,"iter"], cor.trace[span,"pureCor"], xlab="iteration", ylab="pureCorrelation", pch=".", col="blue", ylim = yspan)
      plot(cor.trace[span2,"iter"], cor.trace[span2,"pureCor"], xlab="iteration", ylab="pureCorrelation", pch=20, col="red")
    }
    
    # boosting
    if (boosting) {
      if (iter == 6 + boost*5) {  # Creating an empty matrix
        param.boost <- matrix(NA, nrow = 5, ncol = 9)
      }
      if (iter >= 6 + boost*5 & iter <= 10 + boost*5 ) {  # Storing last ten params
        param.boost[iter - (5 + boost*5),] <- param
      }
      if (iter == 10 + boost*5) {
        param.boost <- booster(param.boost, xvec, yvec, n.cand = min(max(5, index * 2),20))
        tmp.bbbb <<-param.boost
        # print(dim(param.boost)); print(length(param.boost))
        if (showFlag) {print(param.boost)}
        
        if (is.null (dim(param.boost))) {
          param <- param.boost[1:9]
        } else {
          index <- which.max(param.boost[,10])
          param <- param.boost[index,1:9]
          if (showFlag) {print(paste0("Jump to the ",index, "th parameter"))}
        } 
        boost <- boost + 1
      }
    }
    
    #print (expt) #####
    if (showFlag) {cat("iter ", iter, "parm:", round(param,4), if (debug) {c("D.lik=", round(lik - lik.old, 2))},
                       "lik=", expt[1] * n, "p.Cor=", pureCor, "\n")} #lik: lik of previous iteration
    if (maxiter <= iter) {
      lik <- lik.BvZINB4(xvec, yvec, param = param)
      result <- c(param, lik, iter, pureCor)
      names(result) <- c("a0", "a1", "a2", "b1", "b2", paste0("p",1:4), "lik","iter", "pureCor")
      return(result)
      }
    if (max(abs(param - param.old)) <= tol) {
      lik <- lik.BvZINB4(xvec, yvec, param = param)
      result <- c(param, lik, iter, pureCor)
      names(result) <- c("a0", "a1", "a2", "b1", "b2", paste0("p",1:4), "lik","iter", "pureCor")
      return(result)
    }
    if (cor.conv & abs(pureCor - pureCor.old) <= tol) {  # if pureCor is converged, then done!
      lik <- lik.BvZINB4(xvec, yvec, param = param)
      result <- c(param, lik, iter, pureCor)
      names(result) <- c("a0", "a1", "a2", "b1", "b2", paste0("p",1:4), "lik","iter", "pureCor")
      return(result)
    }
  }
  #result <- data.frame(a0 = param[1], a1 = param[2], a2 = param[3], b1 = param[4], b2 = param[5], pi = param[6])
  #return(result)
}
ML.BvZINB4.2b <- function(xvec, yvec, ...) {
  result <- try(ML.BvZINB4.2(xvec,yvec,initial = c(as.numeric(ML.BvNB3(xvec,yvec)), .94,.02,.02,.02), ...))
  if (class(result)=="try-error") {
    result <- c(rep(NA,5+4), NA, 0)
  }
  return(result)
}
if (FALSE) {
  ML.BvZINB4.2(extractor(11), extractor(16),showFlag=TRUE, boosting=FALSE, debug=TRUE) #>-804.44 >10mins
  ML.BvZINB4.2(extractor(1), extractor(5),showFlag=TRUE, boosting=FALSE, debug=TRUE) #>-804.44 >10mins
  ML.BvZINB4.2(extractor(1), extractor(5),showFlag=TRUE)
  ML.BvZINB4.2b(extractor(11), extractor(16),showFlag=TRUE)
}


if (FALSE) {
  # making an empty shell
  MLE.BvZINB4 <- as.data.frame(matrix(NA,1,11))
  names(MLE.BvZINB4) <- c("a0", "a1", "a2", "b1", "b2", paste0("p",1:4), "lik","iter")
  MLE.BvZINB4 <- cbind(MLE.Geneset1$BP[,1:3], MLE.BvZINB4)
  
  # fill in MLE's
  tt(1)
  for (i in 1:dim(MLE.BvZINB4)[1]) {
    if (is.na(MLE.BvZINB4[i,14])) {
      if (a1 < as.numeric(MLE.BvZINB4[i,1])) {
        saveRDS(MLE.BvZINB4, "result-BvZINB4.rds")
      }
      a1 <- as.numeric(MLE.BvZINB4[i,1]); a2 <- as.numeric(MLE.BvZINB4[i,2])
      MLE.BvZINB4[i,4:14] <- ML.BvZINB4.2b(extractor(a1), extractor(a2),maxiter=200)
      print(MLE.BvZINB4[i,])
    }
  }
  saveRDS(MLE.BvZINB4, "result-BvZINB4.rds")
  tt(2)
}




#3. Simulation
rBZINB4 <- function(n, a0, a1, a2, b1, b2, p1, p2, p3, p4, param = NULL) {
  if (!is.null(param)) {
    a0 = param[1]; a1 = param[2]; a2 = param[3]
    b1 = param[4]; b2 = param[5]
    pp = param[6:9]
  } else { pp = c(p1,p2,p3,p4)}
  E <- rmultinom(n, 1, pp)
  E1 <- E[1,] + E[2,]
  E2 <- E[1,] + E[3,]
  x0 <- rgamma(n, a0, 1/b1)
  x1 <- rgamma(n, a1, 1/b1)
  x2 <- rgamma(n, a2, 1/b1)
  x <- rpois(n, x0 + x1)
  y <- rpois(n, (x0 + x2)*b2/b1)
  x <- x*E1
  y <- y*E2
  return(data.frame(x=x, y=y))
}
# rBZINB4.vec <- Vectorize(rBZINB4)
if (FALSE) {
  # param for pair 1 and 2
  param <- c(4.375187e-04,  1.012747e-02,  1.821521e-03,  6.016255e+01,  3.122548e+01,  9.486775e-01,  1.893068e-02,  1.847954e-02,  1.391224e-02)
  set.seed(1)
  tmp <- rBZINB4 (800, param=param)
  table(tmp$x, tmp$y)
  table(extractor(1),extractor(2))
  
  tt(1)
  ML.BvZINB4.2b(tmp$x, tmp$y, maxiter=500, showFlag=TRUE)
  tt(2)
  # 500 iterations 1.09 mins
  # 5.159358e-04  1.092834e-02  3.847828e-03  5.408111e+01  9.785946e+00  9.539740e-01  1.986086e-02 
  # 1.855661e-02  7.608509e-03 -3.202578e+02  5.000000e+02 
  
  
  
  # real data: 500 iterations 1.16 mins
  tt(1)
  ML.BvZINB4.2b(extractor(1), extractor(2), maxiter=500, showFlag=TRUE)
  tt(2)
  
  
  param <- c(1,1,1,1,1,  .95,  .02,  .02,  .01)
  set.seed(2)
  tmp <- rBZINB4 (800, param=param)
  table(tmp)
  tt(1)
  ML.BvZINB4.2b(tmp$x, tmp$y, maxiter=500, showFlag=TRUE)
  tt(2)  # 1.86 mins
  
}
