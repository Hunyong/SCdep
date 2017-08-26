#### Solved! ############################
# Expt function included error (density for nonzero cells did not consider (1-pi)!)



#### Simulation for EM diagnostic
### 0. Symptom
ML.BvZINB2(c(rep(0,9),1),c(rep(0,9),1), showFlag=TRUE)
# At thrid~40th iterations likelihood falls down
# Take 5th iteration
c(0.126980,  0.000010, 0.000010,  2.041530,  0.640850) #5th
c(0.131700,  0.000010, 0.000010,  1.974500,  0.618190) #6th

### 1. M-step check
a1 <- apply(dBvZINB2.Expt.vec(c(rep(0,9),1),c(rep(0,9),1),0.126980,  0.000010, 0.000010,  2.041530,  0.640850),1,mean)
# lik.log = function (param.l) { lik.BvZINB2 (x = c(rep(0,9),1), y = c(rep(0,9),1), exp(param.l))}
# exp(optim(par=log(c(0.126980,  0.000010, 0.000010,  2.041530,  0.640850)), fn = lik.log, control=list(fnscale=-1))$par)

opt.vec <- function(par.ab, expt) {
  par.ab <- exp(par.ab)
  r1 <- sum(expt[2:4]) - sum(par.ab[1:3]) * par.ab[4]
  r2 <- expt[5:7] - digamma(par.ab[1:3]) - log(par.ab[4])
  # print(c(r1,r2)) ###
  return(c(r1,r2))
}
exp(multiroot(opt.vec, start=log(c(0.126980,  0.000010, 0.000010,  2.041530)), expt=a1)$root)
# produces (1.317013e-01 9.999976e-06 9.999976e-06 1.974506e+00)
lik.BvZINB2 (x = c(rep(0,9),1), y = c(rep(0,9),1), c(1.317013e-01, 9.999976e-06, 9.999976e-06, 1.974506e+00,  0.618190))
lik.BvZINB2 (x = c(rep(0,9),1), y = c(rep(0,9),1), c(0.131700,  0.000010, 0.000010,  1.974500,  0.618190))


opt.vec2 <- function(par.ab, expt) {
  par.ab <- exp(par.ab)
  b <- par.ab[4]
  l = sum(par.ab[1:3] * expt[5:7])
  l = l - sum(par.ab[1:3])*log(b)
  l = l - sum(log(gamma(par.ab[1:3])))
  l = l - sum(expt[2:4]*c(1+b*2, 1+ b, 1+b)/b)
  return(-l)
}
exp(optim(par= log(c(0.126980,  0.000010, 0.000010,  2.041530)), opt.vec2, expt=a1)$par)
# produces (1.316441e-01 9.999524e-06 9.996755e-06 1.976954e+00)
lik.BvZINB2 (x = c(rep(0,9),1), y = c(rep(0,9),1), c(1.316441e-01, 9.999524e-06, 9.996755e-06, 1.976954e+00,  0.618190))


[1]  5.000000  0.126980  0.000010  0.000010  2.041530  0.640850 -4.562934
[1]  6.000000  0.131700  0.000010  0.000010  1.974500  0.618190 -4.599228
[1]  7.000000  0.136140  0.000010  0.000010  1.886370  0.599870 -4.629382



### 2. E-step check: OKAY!!!!
# After >5000*5000 Every error was < 2e-4 (mostly ~= e-5 ~ e-6)
#Simulation
n.sim <- 1e+5; n.sim2 <- 1e+5
# ML.BvZINB2(c(rep(0,9),1),c(rep(0,9),1), initial=c(1,1,1,1,.5), showFlag=TRUE)
tmp.1 <- data.frame(n0=1:n.sim2, R0=NA, R1=NA, R2=NA, lR0=NA, lR1 = NA, lR2=NA, Z=NA, V=NA)
tmp.2 <- data.frame(n0=1:n.sim2, R0=NA, R1=NA, R2=NA, lR0=NA, lR1 = NA, lR2=NA, Z=NA, V=NA)
to.be <- t(dBvZINB2.Expt.vec(c(0,1),c(0,1),1,1,1,1,.5))
#      logdensity     xi0.E     xi1.E     xi2.E  log.xi0.E  log.xi1.E  log.xi2.E       z.E      v.E
#[1,] -0.6131045 0.9487179 0.9615385 0.9615385 -0.6617243 -0.6305347 -0.6305347 0.9230769 1.846154
#[2,] -2.7011298 0.6551724 0.7586207 0.7586207 -0.8482417 -0.7531215 -0.7531215 0.0000000 1.000000
set.seed(1)
for (i in 1:1e+5)
  {
  par <- c(1,1,1,1,.5)
  rmat <- matrix(rgamma(n.sim*3, par[1:3], 1/par[4]),n.sim,3)
  rmat2 <- rmat
  rmat2[,3] <- rmat2[,1] + rmat2[,3]
  rmat2[,2] <- rmat2[,1] + rmat2[,2]
  rmat2 <- rmat2[,2:3]
  uv <- matrix(rpois(n.sim*2, rmat2),n.sim,2)
  rmat2[1:10,];uv[1:10,]
  z <- rbinom(n.sim, 1, .5)
  z <- matrix(z,n.sim,2)
  
  xy <- uv * (1-z)
  z[1:10,]; xy[1:10,]
  
  # current solution
  tmp.1[i,2:4] <- apply(rmat[xy[,1] == 0 & xy[,2] == 0,],2,mean) #E[R|(0,0)]
  tmp.2[i,2:4] <- apply(rmat[xy[,1] == 1 & xy[,2] == 1,],2,mean) #E[R|(1,1)]
  tmp.1[i,5:7] <- apply(log(rmat[xy[,1] == 0 & xy[,2] == 0,]),2,mean) #E[log R|(0,0)]
  tmp.2[i,5:7] <- apply(log(rmat[xy[,1] == 1 & xy[,2] == 1,]),2,mean) #E[log R|(1,1)]
  tmp.1[i,8] <- apply(z[xy[,1] == 0 & xy[,2] == 0,],2,mean)[1] #E[Z|(0,0)]
  tmp.2[i,8] <- apply(z[xy[,1] == 1 & xy[,2] == 1,],2,mean)[1] #E[Z|(1,1)]
  tmp.1[i,9] <- apply(uv[xy[,1] == 0 & xy[,2] == 0,],2,mean)[2] #E[V|(0,0)]
  tmp.2[i,9] <- apply(uv[xy[,1] == 1 & xy[,2] == 1,],2,mean)[2] #E[V|(1,1)]
  #print(tmp.1[i,])
  #print(tmp.2[i,])
  tmp.3 <- rbind(c(i,apply(tmp.1[,-1],2,mean,na.rm=TRUE)),
                 c(i,apply(tmp.2[,-1],2,mean,na.rm=TRUE)))
  tmp.3[,-1] <- tmp.3[,-1] - to.be[,-1]
  print(round(tmp.3,7))
}
apply(tmp.1,2,mean,na.rm=TRUE)
apply(tmp.2,2,mean,na.rm=TRUE)


rm(a1,a2,lik.log, n.sim,n.sim2, par, rmat, rmat2, uv, z, xy, tmp.1,tmp.2, opt.vec, to.be)
