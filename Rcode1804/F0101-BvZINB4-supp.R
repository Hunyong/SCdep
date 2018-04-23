#### robustness check of MLE for MI estimation
# stat.BvZINB4 corrected for error: (pureCor should include sqrt(b1*b2/(b1+1)/(b2+1)))
stat.BvZINB4 <- function(a0, a1, a2, b1, b2, p1, p2, p3, p4, param = NULL, measure = "both") {
  if (!is.null(param)) {
    a0 = param[1]; a1 = param[2]; a2 = param[3]
    b1 = param[4]; b2 = param[5]
    pp = param[6:9]
  } else { pp = c(p1,p2,p3,p4)}
  p1 <- pp[1]; p2 <- pp[2]; p3 <- pp[3]; p4 <- pp[4]
  if (measure=="both") {
    cov.x.y <- p1 * b1 * b2 *(a0 + (a0 + a1) * (a0 + a2)) - 
      (p1 + p2) * (p1 + p3) * (a0 + a1) * (a0 + a2) * b1 * b2
    cov.x <- (a0 + a1) * b1 * (p1 + p2) * (b1 + 1 + (a0 + a1) * b1 * (p3 + p4))
    cov.y <- (a0 + a2) * b2 * (p1 + p3) * (b2 + 1 + (a0 + a2) * b2 * (p2 + p4))
    cor <- cov.x.y / sqrt(cov.x * cov.y)    
    MI <- entropy.generic(a0=a0, a1=a1, a2=a2, b1=b1, b2=b2, p1=p1, p2=p2, p3=p3, p4=p4, densityfn=dBvZINB4)
    result <- as.numeric(c(cov.x.y, cov.x, cov.y, cor, MI))
    names(result) <- c("cov.x.y", "cov.x", "cov.y", "cor", "H1", "H2", "H12", "MI")
  }
  if (measure=="cor") {
    cov.x.y <- p1 * b1 * b2 *(a0 + (a0 + a1) * (a0 + a2)) - 
      (p1 + p2) * (p1 + p3) * (a0 + a1) * (a0 + a2) * b1 * b2
    cov.x <- (a0 + a1) * b1 * (p1 + p2) * (b1 + 1 + (a0 + a1) * b1 * (p3 + p4))
    cov.y <- (a0 + a2) * b2 * (p1 + p3) * (b2 + 1 + (a0 + a2) * b2 * (p2 + p4))
    cor <- cov.x.y / sqrt(cov.x * cov.y)    
    result <- as.numeric(c(cov.x.y, cov.x, cov.y, cor))
    names(result) <- c("cov.x.y", "cov.x", "cov.y", "cor")
  }
  if (measure=="pureCor") {
    pureCor <- a0 / sqrt((a0 + a1) * (a0 + a2)) * sqrt(b1*b2 / (b1+1) / (b2+1))
    result <- pureCor
    names(result) <- c("pureCor")
  }
  if (measure=="MI") {
    MI <- entropy.generic(a0=a0, a1=a1, a2=a2, b1=b1, b2=b2, p1=p1, p2=p2, p3=p3, p4=p4, densityfn=dBvZINB4)
    result <- MI
  }
  return(result)
}

if (FALSE) {
  param <- c(1,1,1,1,1,.94,.02,.02,.02)
  stat.BvZINB4(1,1,1,1,2,.25,.25,.25,.25, measure="cor")
  cor = 0.09128709
  
  set.seed(1)
  tmp <- rBvZINB4(800, 1,1,1,1,2,.25,.25,.25,.25)
  ML.BvZINB4.2(tmp[,"x"], tmp[,"y"], maxiter = 10000, showFlag = TRUE, showPlot = TRUE)
  # 0.06 (10000 iterations) vs truth=0.09129
  
  a <- list(tmp, param=c(1,1,1,1,2,.25,.25,.25,.25), data=tmp)
  saveRDS(a, "trace-sim1.rds")
  
  
  set.seed(2)
  tmp <- rBvZINB4(100, 1,1,1,1,2,.25,.25,.25,.25)
  ML.BvZINB4.2(tmp[,"x"], tmp[,"y"], maxiter = 10000, showFlag = TRUE, showPlot = TRUE)
  a <- list(tmp, param=c(1,1,1,1,2,.25,.25,.25,.25), data=tmp)
  saveRDS(a, "trace-sim2.rds")
  # -0.04240468 (2676 iterations) vs truth=0.09129
  
  
  set.seed(3)
  tmp <- rBvZINB4(100, 1,1,1,1,2,.25,.25,.25,.25)
  ML.BvZINB4.2(tmp[,"x"], tmp[,"y"], maxiter = 10000, showFlag = TRUE, showPlot = TRUE)
  a <- list(tmp, param=c(1,1,1,1,2,.25,.25,.25,.25), data=tmp)
  saveRDS(a, "trace-sim3.rds")
  # 0.1058056 (10000 iterations) vs truth=0.09129
  
  set.seed(4)
  tmp <- rBvZINB4(100, 1,1,1,1,2,.25,.25,.25,.25)
  ML.BvZINB4.2(tmp[,"x"], tmp[,"y"], maxiter = 10000, showFlag = TRUE, showPlot = TRUE)
  a <- list(tmp, param=c(1,1,1,1,2,.25,.25,.25,.25), data=tmp)
  saveRDS(a, "trace-sim3.rds")
  # 0.1351091 (10000 iterations) vs truth=0.09129
  
  
  param <- c(4.375187e-04,  1.012747e-02,  1.821521e-03,  6.016255e+01,  3.122548e+01,  9.486775e-01,  1.893068e-02,  1.847954e-02,  1.391224e-02)
  stat.BvZINB4(4.375187e-04,  1.012747e-02,  1.821521e-03,  6.016255e+01,  3.122548e+01,  9.486775e-01,  1.893068e-02,  1.847954e-02,  1.391224e-02, measure="cor")
  cor = 0.08578824
  
  set.seed(11)
  tmp <- rBZINB4(800, param=param)
  ML.BvZINB4.2(tmp[,"x"], tmp[,"y"], maxiter = 10000, showFlag = TRUE, showPlot = TRUE)
  # 0.02318066 (10000 iterations) vs truth=0.085788
  a <- list(tmp, param=param, data=tmp)
  saveRDS(a, "trace-sim2-1.rds")
  
  
  set.seed(12)
  tmp <- rBZINB4(100, param=param)
  ML.BvZINB4.2(tmp[,"x"], tmp[,"y"], maxiter = 10000, showFlag = TRUE, showPlot = TRUE)
  # 0.02318066 (10000 iterations) vs truth=0.085788
  a <- list(tmp, param=param, data=tmp)
  saveRDS(a, "trace-sim2-2.rds")
  
  set.seed(13)
  tmp <- rBZINB4(100, param=param)
  ML.BvZINB4.2(tmp[,"x"], tmp[,"y"], maxiter = 10000, showFlag = TRUE, showPlot = TRUE)
  # -1.141378e-02 (10000 iterations) vs truth=0.085788
  
  set.seed(17)  #seed 14-16 error!
  tmp <- rBZINB4(100, param=param)
  ML.BvZINB4.2(tmp[,"x"], tmp[,"y"], maxiter = 10000, showFlag = TRUE, showPlot = TRUE)
  #  -8.111036e-03 (10000 iterations) vs truth=0.085788
  
  
  set.seed(18)
  tmp <- rBZINB4(100, param=param)
  ML.BvZINB4.2(tmp[,"x"], tmp[,"y"], maxiter = 10000, showFlag = TRUE, showPlot = TRUE)
  #  -8.111036e-03 (10000 iterations) vs truth=0.085788
  
  
  set.seed(3)
  tmp <- rBZINB4(800, param=param)
  ML.BvZINB4.2(tmp[,"x"], tmp[,"y"], maxiter = 10000, showFlag = TRUE, showPlot = TRUE)
  table(tmp)
  ggplot(tmp, aes(x,y)) +
    geom_jitter()
  plot(tmp[,"x"], tmp[,"y"])
  
  
  #### estimation of correlation for a simulated data
  n.sim = 100
  
  ### setting 1: cor = 0.10, % zero pairs = 0.9
  param <- c(.04,.1,.1,1,1, .4,.2,.2,.2)
  stat.BvZINB4(param=param, measure="cor") # cor=0.097
  set.seed(3)
  tmp <- rBZINB4(800, param=param)
  table(tmp$x, tmp$y)[1,1]/800
  
  sim.1 <- list(param = param, 
                cor = stat.BvZINB4(param=param, measure="cor")["cor"],
                sim.stat = data.frame(matrix(NA, n.sim, length(ML.BvZINB4.2(1,1, maxiter = 2))+2)))
  names(sim.1$sim.stat) <- c("no", "n.zero", names(ML.BvZINB4.2(1,1, maxiter = 2)))
  
  set.seed(21)
  for (i in 1:n.sim) {
    tmp <- rBZINB4(200, param=param)
    tmp.1 <- ML.BvZINB4.2(tmp[,"x"], tmp[,"y"], maxiter = 1000, showFlag = TRUE, showPlot = TRUE)
    tmp.2 <- c(no = i, n.zero = table(tmp$x, tmp$y)[1,1]/200, tmp.1)
    sim.1$sim.stat[i,] <- tmp.2
  }
  saveRDS(sim.1, "sim-1-170924.rds")
  
  hist(sim.1$sim.stat$cor)
  mean(sim.1$sim.stat$cor);sd(sim.1$sim.stat$cor)
  sim.1$cor
  apply(sim.1$sim.stat, 2, mean)
  c(sim.1$param, sim.1$cor)
  
  
  ### setting 2: cor = 0.20, % zero pairs = 0.85
  param <- c(.15,.1,.1,1,1, .4,.2,.2,.2)
  stat.BvZINB4(param=param, measure="cor") # cor=0.097
  
  sim.2 <- list(param = param, 
                cor = stat.BvZINB4(param=param, measure="cor")["cor"],
                sim.stat = data.frame(matrix(NA, n.sim, length(ML.BvZINB4.2(1,1, maxiter = 2))+2)))
  names(sim.2$sim.stat) <- c("no", "n.zero", names(ML.BvZINB4.2(1,1, maxiter = 2)))
  
  tt(1)
  set.seed(22)
  for (i in 1:n.sim) {
    tmp <- rBZINB4(200, param=param)
    tmp.1 <- ML.BvZINB4.2(tmp[,"x"], tmp[,"y"], maxiter = 1000, showPlot = TRUE)
    tmp.2 <- c(no = i, n.zero = table(tmp$x, tmp$y)[1,1]/200, tmp.1)
    sim.2$sim.stat[i,] <- tmp.2
    print(i); print(tmp.2)
  }
  tt(2)
  saveRDS(sim.2, "sim-2-170924.rds")
  
  hist(sim.2$sim.stat$cor)
  mean(sim.2$sim.stat$cor);sd(sim.2$sim.stat$cor)
  sim.2$cor
  round(apply(sim.2$sim.stat, 2, mean, trim=0.1),3)
  c(sim.2$param, sim.2$cor)
  
  
  
  
  ### setting 3: cor = 0.21, % zero pairs = 0.78
  param <- c(.25,.15,.15,1,1, .4,.2,.2,.2)
  stat.BvZINB4(param=param, measure="cor") # cor=0.205
  
  sim.3 <- list(param = param, 
                cor = stat.BvZINB4(param=param, measure="cor")["cor"],
                sim.stat = data.frame(matrix(NA, n.sim, length(ML.BvZINB4.2(1,1, maxiter = 2))+2)))
  names(sim.3$sim.stat) <- c("no", "n.zero", names(ML.BvZINB4.2(1,1, maxiter = 2)))
  
  tt(1)
  set.seed(23)
  for (i in 1:n.sim) {
    tmp <- rBZINB4(200, param=param)
    tmp.1 <- ML.BvZINB4.2(tmp[,"x"], tmp[,"y"], maxiter = 1000, showPlot = TRUE)
    tmp.2 <- c(no = i, n.zero = table(tmp$x, tmp$y)[1,1]/200, tmp.1)
    sim.3$sim.stat[i,] <- tmp.2
    print(i); print(tmp.2)
  }
  tt(2)
  saveRDS(sim.3, "sim-3-170924.rds")
  
  hist(sim.3$sim.stat$cor)
  mean(sim.3$sim.stat$cor);sd(sim.3$sim.stat$cor)
  sim.3$cor
  round(apply(sim.3$sim.stat, 2, mean, trim=.1),3)
  c(sim.3$param, sim.3$cor)
  
  
  
  ### Simulation with n.sim=1000
  ### setting 1: cor = 0.10, % zero pairs = 0.9
  a <- readRDS("Rcode/sim-1-170925.rds")
  hist(a$sim.stat$cor)
  mean(a$sim.stat$cor);sd(a$sim.stat$cor)
  a$cor
  round(apply(a$sim.stat, 2, mean, trim=.1),3)
  round(apply(a$sim.stat, 2, mean),3)
  c(a$param, a$cor)
  
  
}
