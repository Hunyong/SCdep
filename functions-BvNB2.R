### 1. Density, likelihood, gradient
dBvNB2 <- function(x, y, a0, a1, a2, b, log=FALSE) {
  p = (2*b + 1) /(b + 1); adj = 0
  l1 <- function(k, m) exp(lgamma(a1 + k) - lgamma(k+1) - lgamma(a1) + lgamma(x + y + a0 -m -k) - lgamma(x -k +1) - lgamma(a0 + y - m) + k *log(p))
  l2 <- function(m) sum(sapply(0:x, l1, m = m))
  #print(sapply(0:x, l1, m = 0)); print(sum(sapply(0:x, l1, m = 0))) ##
  l3 <- function(m) log(l2(m)) + lgamma(m +a2) - lgamma(m+1) - lgamma(a2) + lgamma(y +a0 -m) - lgamma(y -m +1) - lgamma(a0) + m *log(p)
  l4 <- sum(sapply(0:y, function(m) exp(l3(m)))) #%>%print
  if (is.infinite(l4)) {adj = 200; l4 <- sum(sapply(0:y, function(m) exp(l3(m) - adj)))} #%>%print
  l4 <- log(l4) - (+x+y+a0)*log(1 + 2*b) + (x+y) * log(b) -(a1 + a2)*log(1+b) + adj
  return(ifelse(log, l4, exp(l4)))
}
if (FALSE) {
  dBvNB2(1,1,1,1,1,1)
  tmp <- sapply(0:50, function(r) sapply(0:50, function(s) dBvNB2(s,r,1,1,1,1)))
  sum(tmp) #1
}
dBvNB2.vec <- Vectorize(dBvNB2)
lik.BvNB2 <- function(x, y, param) {
  sum(log(dBvNB2.vec(x, y, param[1], param[2], param[3], param[4])))
}

dBvNB2.gr <- function(x, y, a0, a1, a2, b) {
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
  l1.sum <- sum(l1.mat) 
  
  gr.b.mat <- sapply(0:x, function(k) sapply(0:y, gr.b, k=k))
  gr.b.mat <- l1.mat * gr.b.mat
  gr.b.sum <- sum(gr.b.mat)/sum(l1.mat)
  
  gr.a0.mat <- sapply(0:x, function(k) sapply(0:y, gr.a0, k=k))
  gr.a0.mat <- l1.mat * gr.a0.mat
  gr.a0.sum <- sum(gr.a0.mat)/sum(l1.mat)
  
  gr.a1.mat <- matrix(sapply(0:x, gr.a1),x+1, y+1) #%>% print
  gr.a1.mat <- l1.mat * t(gr.a1.mat)
  gr.a1.sum <- sum(gr.a1.mat)/sum(l1.mat)

  gr.a2.mat <- matrix(sapply(0:y, gr.a2), y+1, x+1) #%>% print
  gr.a2.mat <- l1.mat * gr.a2.mat
  gr.a2.sum <- sum(gr.a2.mat)/sum(l1.mat)
  
  result <- list(logdensity=log(l1.sum), gradient = c(gr.a0.sum, gr.a1.sum, gr.a2.sum, gr.b.sum))
  return(result)
}
dBvNB2.gr.vec <- Vectorize(dBvNB2.gr)

tmp <- dBvNB2.gr.vec(c(1,1,1),c(0,1,2),1,1,1,1)[2,]
sapply(tmp, cbind)



### 2. MLE
ML.BvNB2 <- function (xvec, yvec, abstol=1e-8, method="BFGS", showFlag=FALSE) {
  xy.reduced <- as.data.frame(table(xvec,yvec))
  names(xy.reduced) <- c("x", "y","freq")
  xy.reduced <- xy.reduced[xy.reduced$freq != 0,]
  xy.reduced$x <- as.numeric(as.character(xy.reduced$x))
  xy.reduced$y <- as.numeric(as.character(xy.reduced$y))
  xy.reduced$freq <- as.numeric(as.character(xy.reduced$freq))
  if (max(xvec)==0 & max(yvec)==0) {return(rep(1e-10,4))}
  #print(xy.reduced)
  
  fn.1 = function (param) {
    val <- dBvNB2.gr.vec( x = xy.reduced$x, y = xy.reduced$y,  a0 = param[1], a1 = param[2], a2 = param[3], b = param[4])
    lik <- sapply(val[1,],cbind) %*% xy.reduced$freq
    return(lik)
  }
  gr.1 = function (param) {
    val <- dBvNB2.gr.vec( x = xy.reduced$x, y = xy.reduced$y,  a0 = param[1], a1 = param[2], a2 = param[3], b = param[4])
    lik <- sapply(val[1,],cbind) %*% xy.reduced$freq
    gr <- sapply(val[2,],cbind) %*% xy.reduced$freq
    # print(c(param,gr)) ###
    return(gr)
  }
  #print(fn.1(c(1,1,1,1)))
  #print(gr.1(c(1,1,1,1)))
  
  #log-scaled params: param.l
  fn.log = function (param.l) { fn.1 (exp(param.l))}
  gr.log = function (param.l) { 
    if (showFlag) {print(exp(param.l))}
    (as.vector(gr.1 (exp(param.l))) * exp(param.l)) 
    }
  
  # initial guess
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
    
  # result <- multiroot(f = gr.log, start = c(1,1,1,1), positive=TRUE)$root
  # result <- exp(optim(par = log(initial), fn = fn.log, control=list(fnscale=-1), method = "Nelder-Mead")$par)
  result <- try(exp(optim(par = log(initial), fn = fn.log, gr = gr.log, control=list(fnscale=-1, abstol=abstol), method = method)$par), silent=TRUE)
  if (class(result)=="try-error") {
    initial = rep(1,4)
    result <- exp(optim(par = log(initial), fn = fn.log, gr = gr.log, control=list(fnscale=-1, abstol=abstol), method = method)$par)
  }
  result <- data.frame(a0 = result[1], a1 = result[2], a2 = result[3], b = result[4])
  return(result)
}
# simple tests
if (FALSE) {
  ML.BvNB2(c(10,1,1),c(10,1,2))
  a <- ML.BvNB2(c(10,1,1),c(10,1,2))
  lik.BvNB2(c(10,1,1),c(10,1,2), a)  #1.1896784112 0.0001129979 0.0003656534 3.5011913755
  lik.BvNB2(c(10,1,1),c(10,1,2), c(1.18968, 0.1, 0.0003656534, 3.5011913755)) # BFGS w/o gradient better           #-13.082
  lik.BvNB2(c(10,1,1),c(10,1,2), c(1.189528, 1.097285e-14, 4.434082e-08, 3.503274)) #Nelda-Mead w/ gradient best!! #-12.924 computation time is also quick
  lik.BvNB2(c(10,1,1),c(10,1,2), c(1.1529372, 0.4683784, 0.4427861, 3.5977963)) # multirood not good               #-13.93
  tt(1)
  a <- ML.BvNB2(extractor(1), extractor(38), abstol=1e-15)
  tt(2)
  
  ML.BvNB2(extractor(2), extractor(5), method="BFGS", abstol=1e-30)
  ML.BvNB2(extractor(2), extractor(5), method="Nelder-Mead", abstol=1e-30)
  sum(sapply(dBvNB2.gr.vec(extractor(2), extractor(5), 1.572104e-04, 1.985879e-03, 1.680592e-02, 3.230458e+01)[1,],cbind))
  sum(sapply(dBvNB2.gr.vec(extractor(2), extractor(5), 1.559623e-04, 1.989401e-03, 1.681576e-02, 3.229095e+01)[1,],cbind))
  
  a <- ML.BvNB2(extractor(1), extractor(3), method="Nelder-Mead", abstol=1e-30)
  apply(sapply(dBvNB2.gr.vec(extractor(1), extractor(3), 1.932920e-07, 8.667443e-03, 5.954793e-02 ,1.258406e+02)[2,],cbind), 1, sum)
  sum(sapply(dBvNB2.gr.vec(extractor(1), extractor(3), 1.932920e-07, 8.667443e-03, 5.954793e-02 ,1.258406e+02)[1,],cbind))
  
  tt(1)
  ML.BvNB2(extractor(1), extractor(3), method="BFGS")
  tt(2)
  tt(1)
  ML.BvNB2(extractor(1), extractor(38), method="BFGS", showFlag=TRUE)
  tt(2) #39 secs
  lik.BvNB2(extractor(1), extractor(38), c(5.788959e-03, 3.624115e-03, 7.932051e-02, 1.002944e+02))
  
  }
# some deviance tests
if (FALSE) {
  tt(1)
  a <- ML.BvNB2(extractor(11), extractor(16))
  tt(2) #15.3 sec!
  lik.BvNB2(extractor(11), extractor(16),a) #lik = -806
  
  tt(1)
  a <- ML.BvNB2(extractor(11), extractor(16))
  tt(2) #15.3 sec!
  lik.BvNB2(extractor(11), extractor(16),a) #lik = -806
  
  tt(1)
  a <- ML.BvNB2(extractor(30), extractor(50))
  tt(2) #1.99 sec!
  lik.BvNB2(extractor(30), extractor(50),a) #lik = -123
  
  tt(1)
  a <- ML.BvNB2(extractor(8), extractor(11))
  tt(2) # sec!
  
  ML.BvNB2(extractor(1), extractor(38))
  lik.BvNB2(extractor(8), extractor(11),a) #lik = -123
}

tt(1)
a <- pairwise.MLE(data=data[iset.Mm.c2[[1]],], ML.fun = ML.BvNB2, showFlag=TRUE)
tt(2)
ML.BvNB2(extractor(1), extractor(12), showFlag=TRUE)
ML.BvNB2(extractor(51), extractor(12), showFlag=TRUE)

tt(1)
MLE.Geneset1$BvNB2 <- pairwise.MLE(data[iset.Mm.c2[[1]],], ML.fun = ML.BvNB2, showFlag=TRUE)  ## 1.8hrs
tt(2)

ML.BvNB2(0,155)



### 3. Deviance
dev.BvNB2 <- function(xvec, yvec, param = NULL, a0 = NULL, a1 = NULL, a2= NULL, b = NULL) {
  # If params = NULL, apply ML.BvNB2. if else, apply those params
  if (is.null (param)) { 
    if (is.null (a0) | is.null (a1) | is.null (a2) | is.null (b)) {
      param = ML.BvNB2 (xvec = xvec, yvec = yvec)
    }
    else { param = c(a0, a1, a2, b)}
  }
  # Log-likelihood of the BvNB2 model
  lik.model <- lik.BvNB2 (x = xvec, y = yvec, param = param)
  
  # Saturated model BZIP params
  ML.BvNB2.vec <- Vectorize(ML.BvNB2)
  param.sat <- t(ML.BvNB2.vec(xvec, yvec)) #%>%print
  param.sat <- do.call(rbind,param.sat)   #"matrix of lists" into "a matrix"
  #tmp.a <<- param.sat ####
  lik.sat   <- sum(dBvNB2.vec(x= xvec, y = yvec, a0 = param.sat[,1], a1 = param.sat[,2], a2 = param.sat[,3], b = param.sat[,4], log = TRUE))
  
  return(data.frame(model.likelihood = lik.model, satrtd.likelihood = lik.sat, deviance = 2*(lik.sat - lik.model)))
}
dev.BvNB2(extractor(1), extractor(2), param=as.vector(a[,2]))  #459  vs BZIP.B:3150
dev.BvNB2(extractor(1), extractor(5), param=as.vector(a[,5]))  #861  vs BZIP.B:2877
dev.BvNB2(extractor(1), extractor(24), param=as.vector(a[,24]))  #395  vs BZIP.B:5145
MLE.Geneset2$BZIP.B$dev

MLE.Geneset1$BvNB2
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
