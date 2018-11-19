# BvZINB4.2s: BvZINB4.2 + Std error estimation + deviance
# BvZINB4(BvZINB3 + varying zero inflation parameters) + .2(boosting)
library(mle.tools) # for std.err estimates

ML.BvZINB4.3 <- function (xvec, yvec, initial = NULL, tol=1e-8, maxiter=200, showFlag=FALSE, showPlot=FALSE, cor.conv = FALSE) {
  parm = ML.BvZINB4.2(xvec, yvec, initial = initial, tol=tol, maxiter=maxiter, showFlag=showFlag, showPlot=showPlot, cor.conv = cor.conv)
  var.r = var.rho(param[1:9], xvec, yvec)  # variance of pureCor
  dev = dev.BvZINB4(xvec, yvec, model.lik = parm[10])
  
  parm = c(parm, sePureCor = sqrt(var.r), dev = dev[3], df = dev[4])
  
    dev.BvZINB4 <- function(xvec, yvec, model.lik = NULL,
                            param = NULL, 
                            a0 = NULL, a1 = NULL, a2= NULL, b1 = NULL, b2 = NULL,
                            p1 = NULL, p2 = NULL, p3 = NULL, p4 = NULL) 