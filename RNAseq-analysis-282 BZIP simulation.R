# Saved results
# readRDS("sim_result_BZIP_170615/sim_result_BZIP(nsim 30).rds")
# readRDS("sim_result_BZIP_170615/sim_result_BP(nsim 30).rds")
# readRDS("sim_result_BZIP_170615/sim_result_EMI(nsim 30).rds")
# readRDS("sim_result_BZIP_170615/param_BZIP(nsim 30).rds")
library(reshape2)

param.BZIP <- expand.grid(c(.7,.8,.9), 10^-(2:0), 10^-(2:0), 10^-(2:0))
names(param.BZIP) <- c("pi", "m0", "m1", "m2")
param.BZIP <- param.BZIP[param.BZIP$m1 >= param.BZIP$m2,]
dim(param.BZIP)
n.sim = 100
sim.result.BZIP <- sim.result.BP <- sim.result.EMI <- list()

tt(1) #time check
set.seed(100)
tmp.2BZIP <- data.frame(pp = NA, mu0 = NA, mu1 = NA, mu2 = NA)
tmp.2BP <- data.frame(mu0 = NA, mu1 = NA, mu2 = NA)
tmp.2EMI <- data.frame(EMI = NA)

for (i in 1:dim(param.BZIP)) {
  tmp.1 <- rBZIP.vec(rep(100,n.sim), param.BZIP$pi[i], param.BZIP$m0[i], param.BZIP$m1[i], param.BZIP$m2[i])
  for (j in 1:n.sim) {
    print(c(j, "out of ",n.sim, "//"))
    print(c(i,"th param :",as.numeric(param.BZIP[i,])))
    tmp.2BZIP[j,] <- ML.BZIP(tmp.1[,j]$x,tmp.1[,j]$y, showFlag=FALSE)
    tmp.2BP[j,] <- ML.BP(tmp.1[,j]$x,tmp.1[,j]$y)
    tmp.2EMI[j,] <- mi.empirical.2(tmp.1[,j]$x,tmp.1[,j]$y)
  }
  sim.result.BZIP[[i]] <- tmp.2BZIP
  sim.result.BP[[i]] <- tmp.2BP
  sim.result.EMI[[i]] <- tmp.2EMI
}
tt(2) #21 hrs

# Adding PMI(BZIP) result to param.BZIP
tmp.3 <- data.frame(t(sapply(1:54,function(s) apply(sim.result.BZIP[[s]],2,mean,na.rm=T))))
names(tmp.3)[1:4] <- c("pp.est","m0.est","m1.est","m2.est")
tmp.3$match <- abs(tmp.3$pp.est / param.BZIP$pi - 1) < 0.2
param.BZIP <- cbind(param.BZIP,tmp.3[,1:5])
param.BZIP[param.BZIP$match==F,]

# Inspecting discrepancy
set.seed(1)
A <- rBZIP.vec(100, .7, .1, .1, .1)
tmp <- ML.BZIP(A[,1]$x,A[,1]$y) %>% print    # Estimate is far away from truth
lik.BZIP(A[,1]$x,A[,1]$y,tmp)           
lik.BZIP(A[,1]$x,A[,1]$y,c(.7,.1,.1,.1))     # Likelihood is smaller for truth.
rm(A)

# sim.result.BZIP.backup <- sim.result.BZIP  # back up
# sim.result.BZIP <- sim.result.BZIP.backup  # restore
## PMI estimates(BZIP)
for (i in 1:dim(param.BZIP)[1]) {
  print(i); tt(1)
  sim.result.BZIP[[i]] <- sim.result.BZIP[[i]][complete.cases(sim.result.BZIP[[i]]),] # remove NA rows
  b <- entropy.BZIP.vec(pp = sim.result.BZIP[[i]]$pp, 
                        m0 = sim.result.BZIP[[i]]$mu0, 
                        m1 = sim.result.BZIP[[i]]$mu1, 
                        m2 = sim.result.BZIP[[i]]$mu2)
  b <- matrix(as.numeric(t(b)),ncol=4)
  sim.result.BZIP[[i]] <- cbind(sim.result.BZIP[[i]], b)
  names(sim.result.BZIP[[i]])[5:8] <- c("H1", "H2", "H12", "MI")
  print(sim.result.BZIP[[i]])
  tt(2) %>% print
}

## PMI estimates(BP)
for (i in 1:dim(param.BZIP)[1]) {
  print(i); tt(1)
  # sim.result[[i]] <- sim.result[[i]][complete.cases(sim.result[[i]]),] # remove NA rows
  b <- entropy.BP.vec(m0 = sim.result.BP[[i]]$mu0, 
                      m1 = sim.result.BP[[i]]$mu1, 
                      m2 = sim.result.BP[[i]]$mu2)
  b <- matrix(as.numeric(t(b)),ncol=4)
  sim.result.BP[[i]] <- cbind(sim.result.BP[[i]], b)
  names(sim.result.BP[[i]])[4:7] <- c("H1", "H2", "H12", "MI")
  print(sim.result.BP[[i]])
  tt(2) %>% print  #4.5 min * 54 = 4.5 hr
}


## True MI
b <- entropy.BZIP.vec(pp = param.BZIP$pi, m0 = param.BZIP$m0, m1 = param.BZIP$m1, m2 = param.BZIP$m2)
b <- matrix(as.numeric(t(b)),ncol=4)
param.BZIP <- cbind(param.BZIP, b)
names(param.BZIP)[10:13] <- c("H1","H2","H12","MI")
param.BZIP$MI.sim.BZIP.mean <- sapply(1:54,function(x) mean(sim.result.BZIP[[x]]$MI))
param.BZIP$MI.sim.BP.mean <- sapply(1:54,function(x) mean(sim.result.BP[[x]]$MI))
param.BZIP$MI.sim.EMI.mean <- sapply(1:54,function(x) mean(sim.result.EMI[[x]]$EMI))
#param.BZIP[,-c(4:7)]
# MSE : PMI(BZIP) performs best on average

param.BZIP$MI.sim.BZIP.MSE <- sapply(1:54,function(s) MSE(sim.result.BZIP[[s]]$MI,param.BZIP[s,"MI"]))
param.BZIP$MI.sim.BP.MSE <- sapply(1:54,function(s) MSE(sim.result.BP[[s]]$MI,param.BZIP[s,"MI"]))
param.BZIP$MI.sim.EMI.MSE <- sapply(1:54,function(s) MSE(sim.result.EMI[[s]]$EMI,param.BZIP[s,"MI"]))
param.BZIP$MI.sim.BZIP.SD <- sapply(1:54,function(s) sd(sim.result.BZIP[[s]]$MI))
param.BZIP$MI.sim.BP.SD <- sapply(1:54,function(s) sd(sim.result.BP[[s]]$MI))
param.BZIP$MI.sim.EMI.SD <- sapply(1:54,function(s) sd(sim.result.EMI[[s]]$EMI))
#param.BZIP$param.no <- rep(1:54)
param.BZIP$param.no <- paste0("(", format(param.BZIP$m0, nsmall = 2), ", ", 
                              format(param.BZIP$m1, nsmall = 2), ", ", 
                              format(param.BZIP$m2, nsmall = 2), ")")
param.BZIP$param.no <- paste0("10^(", log10(param.BZIP$m0), ", ", 
                              log10(param.BZIP$m1), ", ", 
                              log10(param.BZIP$m2), ")")

# Bias plot : EMI is over estimating in general.
param.BZIP.long <- melt(param.BZIP[,c(1:4,13:16,23)],id.vars=names(param.BZIP)[c(1:4,23)])
plot2.8.2B.1 <- ggplot(param.BZIP.long, aes(x=param.no, y=value, fill = variable)) + 
  geom_bar(position = "dodge", stat="identity", width=1) +
  facet_grid(pi~., labeller = as_labeller(c("0.7" = "pi=0.7", "0.8" = "pi=0.8", "0.9" = "pi=0.9"))) +
  ylab("Mean of MI estiamtors") +
  xlab("parameter settings") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual("Estimators",values=c("#CC79A7", "#0072B2", "#009E73", "#D55E00")) +
  scale_x_discrete(limits=param.BZIP$param.no)
png2(plot2.8.2B.1)


#SD
param.BZIP.long <- melt(param.BZIP[,c(1:4,20:23)],id.vars=names(param.BZIP)[c(1:4,23)])
plot2.8.2B.2 <- ggplot(param.BZIP.long, aes(x=param.no, y=value, fill = variable)) + 
  geom_bar(position = "dodge", stat="identity", width=0.5) +
  facet_grid(pi~., labeller = as_labeller(c("0.7" = "pi=0.7", "0.8" = "pi=0.8", "0.9" = "pi=0.9"))) +
  ylab("SD (std. dev) of MI estiamtors") +
  xlab("parameter settings") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual("Estimators",values=c("#0072B2", "#009E73", "#D55E00"))+
  scale_x_discrete(limits=param.BZIP$param.no)
png2(plot2.8.2B.2)

xvec <- extractor(1); yvec <- extractor(2)
table(xvec,yvec)
dev.BZIP.B(xvec, yvec)#dev = 841    vs 861 (reference chi2(800-6, 95%))
dev.BZIP(xvec, yvec)  #dev = 1040   vs 863 (reference chi2(800-4, 95%))
dev.BZIP(xvec, yvec, c(0.9525, 7.699164e-09, 14.52631, 0.6842103))
dev.BP(xvec, yvec)    #dev = 4201   vs 864
qchisq(.95,796) #863
qchisq(.95,797) #864

xvec <- extractor(3); yvec <- extractor(4)   # This has some problem with large number (338,4) that makes the likelihood -Inf screws everything
dev.BZIP(xvec, yvec)  #dev = 1,694  vs 863 (reference chi2(800-4, 95%))
dev.BP(xvec, yvec)    #dev = 17,815 vs 864
dev.BZIP.B(xvec, yvec) 

xvec <- extractor(5); yvec <- extractor(6)
dev.BZIP.B(xvec, yvec)#dev = 1386   vs 861 (reference chi2(800-6, 95%))
dev.BZIP(xvec, yvec)  #dev = 2030   vs 863 (reference chi2(800-4, 95%))
dev.BP(xvec, yvec)    #dev = 5074   vs 864

xvec <- extractor(8); yvec <- extractor(11)
dev.BZIP.B(xvec, yvec)#dev = 4279   vs 861 (reference chi2(800-6, 95%))
dev.BZIP(xvec, yvec)  #dev = 7088   vs 863 (reference chi2(800-4, 95%))
dev.BP(xvec, yvec)    #dev = 18925   vs 864

xvec <- extractor(30); yvec <- extractor(50)
dev.BZIP.B(xvec, yvec)#dev = 239   vs 861 (reference chi2(800-6, 95%))
dev.BZIP(xvec, yvec)  #dev = 250   vs 863 (reference chi2(800-4, 95%))
dev.BP(xvec, yvec)    #dev = 486   vs 864


# BZIP is sig'ly better than BP
# BZIP is not sig'ly better than saturated model
# Try BvZINB!!


# BZIP simulated data example
dev.BZIP(tmp.1[,1]$x, tmp.1[,1]$y) # 111
dev.BP(tmp.1[,1]$x, tmp.1[,1]$y)   # 186
qchisq(.95,96) #119
qchisq(.95,97) #120


# Saving data files
# saveRDS(sim.result.BZIP,"sim_result_BZIP(nsim 30).rds")
# saveRDS(sim.result.BP,"sim_result_BP(nsim 30).rds")
# saveRDS(param.BZIP,"param_BZIP(nsim 30).rds")
# saveRDS(sim.result.EMI,"sim_result_EMI(nsim 30).rds")
# saveRDS(sim.result.BZIP,"sim_result_BZIP(nsim100).rds")
# saveRDS(sim.result.BP,"sim_result_BP(nsim100).rds")
# saveRDS(param.BZIP,"param_BZIP(nsim100).rds")
# saveRDS(sim.result.EMI,"sim_result_EMI(nsim100).rds")


a<- readRDS("sim_result_BZIP_170615/param_BZIP(nsim 30).rds")
a[a$match==FALSE,]

# Check 8th setting: SD(BZIP) > SD(BP)
sim.result.BZIP[[8]][1:5,]
sim.result.BP[[8]][1:5,]
round(param.BZIP[,c(1:8,13,14)],3)


# Relation between Pi estimates and MI estimates
plot2.8.2B.5 <- ggplot(param.BZIP, aes(x=pi-pp.est, y=MI.sim.BZIP.mean/MI)) + 
  geom_point() + ylab("ratio of MI to MI.hat") + xlab("pi - pi.hat")
plot2.8.2B.6 <- ggplot(param.BZIP, aes(x=pi-pp.est, y=MI.sim.BZIP.mean-MI)) + 
  geom_point() + ylab("MI.hat - MI") + xlab("pi - pi.hat")
png2(plot2.8.2B.5); png2(plot2.8.2B.6)

#High ratio >2 points: 28th and 29th: MI=0.006, 0.005 while PMI(BZIP) = 0.014, 0.013
