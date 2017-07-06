# brief test of bias of empirical MI under Bivariate ZIP
# params
mu0 = 0.11; mu1 = 1e-8; mu2 = 1e-8
MI.param = MI.BP (mu0, mu1, mu2) %>% print   # theoretical MI = 0.3569

rBP.MI.sim <- function(n.sim = 1000, n = 100, mu0, mu1, mu2, ML = FALSE, measure = mi.empirical.2) {
  result <- result2 <- list()
  for (i in 1:n.sim) {
    print(c(i," of ",n.sim))
    common.tmp <- rpois(n,mu0)
    x.tmp <- rpois(n,mu1) #%>% print
    y.tmp <- rpois(n,mu2) #%>% print
    x.tmp <- common.tmp + x.tmp    #<<- for debug!!!
    y.tmp <- common.tmp + y.tmp    #<<- for debug!!!
    result[[i]] <- measure(x.tmp, y.tmp) # %>% print
    # print(table(x.tmp,y.tmp)) #debug
    
    if (ML) {
      result2[[i]] <- MI.ML.BP (xvec = x.tmp, yvec = y.tmp) %>% print
    }
  }
  if (ML) return(list(result, result2)) else return(result)
}

set.seed(100); a <- Sys.time()
EMI.bias <- rBP.MI.sim(n.sim = 100, n =100,  mu0=mu0, mu1=mu1, mu2=mu2, ML = TRUE)
Sys.time()-a
png("plot2-8-1.png", width = 720, height = 360)  # histogram of EMI and PMI each
par(mfrow=c(1,2))
hist(unlist(EMI.bias[[1]]), xlab = "MI (truth = blue dotted line)", main= "EMI"); abline(v=c(0.3460, 0.3569), lty=c(1,2), col=c("red","blue"))
hist(unlist(EMI.bias[[2]]), xlab = "MI (truth = blue dotted line)", main= "parametric MI"); abline(v=c(0.3527, 0.3569), lty=c(1,2), col=c("red","blue"))
dev.off()

# statistics for EMI and PMI
mean(unlist(EMI.bias[[1]]))  # mean(EMI) = 0.346    but truth is 0.3569
sd(unlist(EMI.bias[[2]]))    # sd (EMI)  = 0.07
mean(unlist(EMI.bias[[2]]))  # mean(parametric MI) = 0.353    very close to the truth 0.258
sd(unlist(EMI.bias[[2]]))    # sd (parametric MI)  = 0.069
png("plot2-8-1B.png", width = 720, height = 360) # plot comparing EMI and PMI
plot(EMI.bias[[2]],EMI.bias[[1]], xlab="Parametric MI(BP)", ylab = "Empirical MI", xlim = c(0,1.2), ylim = c(0, 0.6))
abline(coef=c(0,1), col="red")
dev.off()


## General settings (BP)
grid.sim <- 10^-c(2:0)  #0.01, 0.1, 1.0
EMI.bias2 <- list(EMI.mean = array(, dim=c(3,3,3)), EMI.sd = array(, dim=c(3,3,3)), PMI.mean = array(, dim=c(3,3,3)), PMI.sd= array(, dim=c(3,3,3)))
EMI.bias2.raw <- list()
tmp.2 <- tmp.3 <- list()
set.seed(200); tt(1)
for (i in 1:length(grid.sim)) { # mu0 loop
  for (j in 1:length(grid.sim)) { # mu1 loop
    for (k in 1:length(grid.sim)) { # mu2 loop
      print(c(paste0("mu",0:2),c(i,j,k)))
      if (j <= k) {
        tmp.1 <- rBP.MI.sim(n.sim = 30, n =100,  mu0 = grid.sim[i], mu1 = grid.sim[j], mu2 = grid.sim[k], ML = TRUE, measure = mi.empirical.2)
        EMI.bias2$EMI.mean[k,j,i] <- mean(unlist(tmp.1[[1]])) # EMI
        EMI.bias2$EMI.sd[k,j,i] <- sd(unlist(tmp.1[[1]]))
        EMI.bias2$PMI.mean[k,j,i] <- mean(unlist(tmp.1[[2]])) # parametric
        EMI.bias2$PMI.sd[k,j,i] <- sd(unlist(tmp.1[[2]]))
        tmp.2[[k]] <- (tmp.1)
      } else {
        EMI.bias2$EMI.mean[k,j,i] <- EMI.bias2$EMI.mean[j,k,i] # storing symmetric components without redundant calculation
        EMI.bias2$EMI.sd[k,j,i] <- EMI.bias2$EMI.sd[j,k,i]
        EMI.bias2$PMI.mean[k,j,i] <- EMI.bias2$PMI.mean[j,k,i]
        EMI.bias2$PMI.sd[k,j,i] <- EMI.bias2$PMI.sd[j,k,i]
      } 
    }
    tmp.3[[j]] <- tmp.2
  }
  EMI.bias2.raw[[i]] <- tmp.2
}
dimnames(EMI.bias2$EMI.mean) <- list(paste0("mu1=",grid.sim),paste0("mu2=",grid.sim),paste0("mu0=",grid.sim))
dimnames(EMI.bias2$EMI.sd) <- dimnames(EMI.bias2$PMI.mean) <- dimnames(EMI.bias2$PMI.sd) <- dimnames(EMI.bias2$EMI.mean)
tt(2)
rm(tmp.1,tmp.2,tmp.3)

# true MI
EMI.bias2.TRUTH <- array(,c(3,3,3))
for (i in 1:3) {
  EMI.bias2.TRUTH[,,i] <- array(unlist(lapply(grid.sim, function(s) {sapply(grid.sim, function(x) {MI.BP(grid.sim[i], s, x)})})),c(3,3))
}
# upper triangular index
tmp.a <- array(,c(3,3,3)); for (i in 1:3) (tmp.a[,,i] <- upper.tri(tmp.a[,,1],diag=TRUE))
tmp.b <- array(,c(3,3,3)); for (i in 1:3) {for(j in 1:3) {for(k in 1:3) {tmp.b[k,j,i] <- paste0(10^(i-3),"-",10^(j-3),"-",10^(k-3))}}}

EMI.bias2$summary <- data.frame(EMI = as.numeric(EMI.bias2$EMI.mean), PMI =as.numeric(EMI.bias2$PMI.mean), MI.true = as.numeric(EMI.bias2.TRUTH))
EMI.bias2$summary <- EMI.bias2$summary[tmp.a,]
EMI.bias2$summary$mu012 <- tmp.b[tmp.a]
EMI.bias2$summary <- EMI.bias2$summary[,c("mu012","MI.true","PMI","EMI")]
EMI.bias2$summary$PMI.sd <- EMI.bias2$PMI.sd[tmp.a]
EMI.bias2$summary$EMI.sd <- EMI.bias2$EMI.sd[tmp.a]

cbind(EMI.bias2$summary$mu012,round(EMI.bias2$summary[,2:4],4),round(EMI.bias2$summary[,5:6],3))
plot2.8.2 <- ggplot(EMI.bias2$summary, aes(1:18)) + 
  geom_line(aes(y = MI.true, colour = "true MI")) + 
  geom_line(aes(y = EMI, colour = "EMI")) +
  geom_line(aes(y = PMI, colour = "PMI"))
png2(plot2.8.2)



################ Bias of EMI for binned data EMI.bias3

## General settings (BP)
grid.sim <- 10^-c(2:0)  #0.01, 0.1, 1.0
EMI.bias3 <- list(EMI.mean = array(, dim=c(3,3,3)), EMI.sd = array(, dim=c(3,3,3)), EMI.St.mean = array(, dim=c(3,3,3)), EMI.St.sd = array(, dim=c(3,3,3)), EMI.Sc.mean = array(, dim=c(3,3,3)), EMI.Sc.sd = array(, dim=c(3,3,3)))
EMI.bias3.raw <- list()
tmp.2 <- tmp.3 <- list()
set.seed(200); tt(1) # time check
for (i in 1:length(grid.sim)) { # mu0 loop
  for (j in 1:length(grid.sim)) { # mu1 loop
    for (k in 1:length(grid.sim)) { # mu2 loop
      print(c(paste0("mu",0:2),c(i,j,k)))
      if (j <= k) {
        # mi.empirical.3 = Sturge, ML=FALSE = no PMI, mi.empirical.3.SC = Scott
        tmp.1 <- rBP.MI.sim(n.sim = 30, n =1000,  mu0 = grid.sim[i], mu1 = grid.sim[j], mu2 = grid.sim[k], ML = FALSE, measure = mi.empirical.2) 
        tmp.1.St <- rBP.MI.sim(n.sim = 30, n =1000,  mu0 = grid.sim[i], mu1 = grid.sim[j], mu2 = grid.sim[k], ML = FALSE, measure = mi.empirical.3) 
        tmp.1.Sc <- rBP.MI.sim(n.sim = 30, n =1000,  mu0 = grid.sim[i], mu1 = grid.sim[j], mu2 = grid.sim[k], ML = FALSE, measure = mi.empirical.3.SC) 
        EMI.bias3$EMI.mean[k,j,i] <- mean(unlist(tmp.1)) # unbinned
        EMI.bias3$EMI.sd[k,j,i] <- sd(unlist(tmp.1))
        EMI.bias3$EMI.St.mean[k,j,i] <- mean(unlist(tmp.1.St)) # Sturge
        EMI.bias3$EMI.St.sd[k,j,i] <- sd(unlist(tmp.1.St))
        EMI.bias3$EMI.Sc.mean[k,j,i] <- mean(unlist(tmp.1.Sc)) # Scott
        EMI.bias3$EMI.Sc.sd[k,j,i] <- sd(unlist(tmp.1.Sc))
        tmp.2[[k]] <- list(original = unlist(tmp.1), Sturge = unlist(tmp.1.St), Scott = unlist(tmp.1.Sc))
      } else {
        EMI.bias3$EMI.mean[k,j,i] <- EMI.bias3$EMI.mean[j,k,i] # storing symmetric components without redundant calculation
        EMI.bias3$EMI.sd[k,j,i] <- EMI.bias3$EMI.sd[j,k,i]
        EMI.bias3$EMI.St.mean[k,j,i] <- EMI.bias3$EMI.St.mean[j,k,i] # storing symmetric components without redundant calculation
        EMI.bias3$EMI.St.sd[k,j,i] <- EMI.bias3$EMI.St.sd[j,k,i]
        EMI.bias3$EMI.Sc.mean[k,j,i] <- EMI.bias3$EMI.Sc.mean[j,k,i]
        EMI.bias3$EMI.Sc.sd[k,j,i] <- EMI.bias3$EMI.Sc.sd[j,k,i]
      } 
    }
    tmp.3[[j]] <- tmp.2
  }
  EMI.bias3.raw[[i]] <- tmp.2
}
dimnames(EMI.bias3$EMI.mean) <- dimnames(EMI.bias3$EMI.sd) <- list(paste0("mu1=",grid.sim),paste0("mu2=",grid.sim),paste0("mu0=",grid.sim))
dimnames(EMI.bias3$EMI.St.mean) <-dimnames(EMI.bias3$EMI.St.sd) <- dimnames(EMI.bias3$EMI.Sc.mean) <- dimnames(EMI.bias3$EMI.Sc.sd) <- dimnames(EMI.bias3$EMI.mean)
tt(2)
rm(tmp.1, tmp.1.St, tmp.1.Sc, tmp.2, tmp.3)

#EMI.bias2.TRUTH # true MI

EMI.bias3$summary <- data.frame(EMI = as.numeric(EMI.bias3$EMI.mean), EMI.St = as.numeric(EMI.bias3$EMI.St.mean), EMI.Sc =as.numeric(EMI.bias3$EMI.Sc.mean), MI.true = as.numeric(EMI.bias2.TRUTH))
EMI.bias3$summary <- EMI.bias3$summary[tmp.a,]
EMI.bias3$summary$mu012 <- tmp.b[tmp.a]
EMI.bias3$summary <- EMI.bias3$summary[,c("mu012","MI.true","EMI", "EMI.St", "EMI.Sc")]
EMI.bias3$summary$EMI.sd <- EMI.bias3$EMI.sd[tmp.a]
EMI.bias3$summary$EMI.St.sd <- EMI.bias3$EMI.St.sd[tmp.a]
EMI.bias3$summary$EMI.Sc.sd <- EMI.bias3$EMI.Sc.sd[tmp.a]
cbind(EMI.bias3$summary$mu012,round(EMI.bias3$summary[,2:5],4),round(EMI.bias3$summary[,6:8],3))

#EMI.bias3$summary$PMI <- EMI.bias2$summary$PMI
plot2.8.3 <- ggplot(EMI.bias3$summary, aes(1:18)) + 
  geom_line(aes(y = MI.true, colour = "true MI")) + 
  geom_line(aes(y = EMI, colour = "EMI(original)")) +
  geom_line(aes(y = EMI.St, colour = "EMI(Sturge)")) +
  geom_line(aes(y = EMI.Sc, colour = "EMI(Scott)"))
png2(plot2.8.3)
plot2.8.3B <- ggplot(EMI.bias3$summary, aes(1:18)) + 
  geom_line(aes(y = EMI, colour = "EMI(original)")) +
  geom_line(aes(y = EMI.St, colour = "EMI(Sturge)"))
png2(plot2.8.3B)




################ Bias of EMI for length
mi.empirical.2.len <- function(x, y, length=list(500, 600, 700, 800, 900, 1000, 2500, 5000, 7500, 10000)) {
  result <- sapply (length, function(i) {mi.empirical.2(x[1:i], y[1:i])})
  names(result) <- unlist(length)
  return(result)
}

set.seed(1)
tmp.1 <- rBP.MI.sim(n.sim = 100, n = 10000,  mu0 = .01,  mu1 = .01, mu2 = .01, ML = FALSE, measure = mi.empirical.2.len)
tmp.2 <- apply(do.call(rbind,tmp.1),2,mean)
png("plot2-8-4.png"); plot(tmp.2, xlab="sample size"); dev.off()



set.seed(2)
tmp.1 <- rBP.MI.sim(n.sim = 100, n = 10000,  mu0 = .1,  mu1 = .1, mu2 = .1, ML = FALSE, measure = mi.empirical.2.len)
apply(do.call(rbind,tmp.1),2,mean)
plot(apply(do.call(rbind,tmp.1),2,mean), xlab = "sample size")
  