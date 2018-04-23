#### Numerical derivation of Hessian for var(correlation)

library(numDeriv) # for hessian
#library(nlme)     # for hessian (fdHess)
lik.BvZINB4.reorder = function(param, xvec, yvec) lik.BvZINB4(x=xvec, y=yvec, param = param)
hess.BvZINB4 <- function (param, xvec, yvec) {numDeriv::hessian(lik.BvZINB4.reorder, x = param, xvec=xvec, yvec=yvec,
                                                                method.args=list(eps=min(param[1:5])*.3, d=0.1, zero.tol=sqrt(.Machine$double.eps/7e-7), r=4, v=2, show.details=FALSE))}

if (FALSE) {fdHess(a[1:9], lik.BvZINB4.reorder, 1:5, 2:6)} # nlme package

pureCor <- function(param, xvec, yvec) {
  
  ## 0. renaming
  a0 <- param[1]
  a1 <- param[2]
  a2 <- param[3]
  b1 <- param[4]
  b2 <- param[5]
  
  ## 1. pureCor
  rho = a0 / sqrt((a0 + a1) * (a0 + a2)) * sqrt(b1*b2 / (b1+1) / (b2+1))
  
  ## 2. variance of parameteres alphas and betas
  hess <- hess.BvZINB4(param = param,  xvec=xvec, yvec=yvec)
  var.ab <- solve(-hess[1:5,1:5])  # take inverse to get the variance
  if (any(diag(hess)[1:5] > 0)) {
    warning("Numerical calculation of Hessian is not negative-definite.")
    return(c(pureCor = rho, se.rho = NA, p.val = NA))
  }
# print(hess)  
  ## 3. Delta method - get gradient of correlation function
  grad.a0 = (a0 * a1 + a0 * a2 + 2* a1 * a2)/(2*a0 *(a0 + a1) *(a0 + a2)) * rho
  grad.a1 = - rho / 2 / (a0 + a1)
  grad.a2 = - rho / 2 / (a0 + a2)
  grad.b1 = - rho / 2 / b1 / (b1 + 1)
  grad.b2 = - rho / 2 / b2 / (b2 + 1)
  
  grad = c(grad.a0, grad.a1, grad.a2, grad.b1, grad.b2)

  ## 4. Delta method - result
  var.rho = t(grad) %*% var.ab %*% grad
  se.rho = as.numeric(sqrt(var.rho))
  p = 2*(1-pnorm(rho/se.rho))
  
  return(c(pureCor = rho, se.rho = se.rho, p.val = p))
}

if (FALSE) {
  # numerically instable probably for small sample
  a <- ML.BvZINB4.2(1, 1)
  pureCor( a[1:9], 1,1)
  
  a <- ML.BvZINB4.2(1:5, 2:6)
  pureCor( a[1:9], 1:5,2:6)
  
  # But works well for large sample
  a <- ML.BvZINB4.2(extractor(1), extractor(11), maxiter=20)
  tt(1); pureCor( a[1:9], extractor(1), extractor(11)); tt(2)
  
  a <- ML.BvZINB4.2(extractor(3), extractor(7), maxiter=20)
  tt(1); pureCor( a[1:9], extractor(3), extractor(7)); tt(2)
  
  a <- ML.BvZINB4.2(extractor(15), extractor(20), maxiter=20)
  tt(1); pureCor( a[1:9], extractor(15), extractor(20)); tt(2)
  
  pureCor( c(1,1,1,1,1,.25,.25,.25,.25), 1:5,2:6)
}

# argument added and order changed (2018-04-15) !!!
# if model likelihood is given, no need to calculate it.
dev.BvZINB4 <- function(xvec, yvec, model.lik = NULL,
                        param = NULL, 
                        a0 = NULL, a1 = NULL, a2= NULL, b1 = NULL, b2 = NULL,
                        p1 = NULL, p2 = NULL, p3 = NULL, p4 = NULL) {
  
  if (is.null(model.lik)) {
    # If params = NULL, apply ML.BvNB3. if else, apply those params
    if (is.null (param)) { 
      if (is.null (a0) | is.null (a1) | is.null (a2) | is.null (b1) | is.null (b2)) {
        param = ML.BvZINB4.2 (xvec = xvec, yvec = yvec)
      }
      else { param = c(a0, a1, a2, b1, b2, p1, p2, p3, p4)}
    }
    # Log-likelihood of the BvNB2 model
    model.lik <- lik.BvZINB4 (x = xvec, y = yvec, param = param)
  }

  # Reduced calculation
  xy.reduced <- as.data.frame(table(xvec,yvec))
  names(xy.reduced) <- c("x", "y","freq")
  xy.reduced <- xy.reduced[xy.reduced$freq != 0,]
  xy.reduced$x <- as.numeric(as.character(xy.reduced$x))
  xy.reduced$y <- as.numeric(as.character(xy.reduced$y))
  xy.reduced$freq <- as.numeric(as.character(xy.reduced$freq))
  # Saturated model BZIP params

  param.sat <- sapply(1:(dim(xy.reduced)[1]), function(i) {
    x=xy.reduced$x[i]; y=xy.reduced$y[i]
# print(c(x=xy.reduced$x[i], y=xy.reduced$y[i]))
    if (x * y == 0 & x + y > 30) {  # For a (zero, large nonzero) pair, ML.BvZINB4.2 fails. replaced by ML.BvNB3
      c(as.numeric(ML.BvNB3(x, y))[1:5], 0,x>0,y>0,0) # %>% print
    } else if (x * y > 0) {
      c(as.numeric(ML.BvNB3(x, y))[1:5], 1,0,0,0) # %>% print
    } else {
      as.numeric(ML.BvZINB4.2(xy.reduced$x[i], xy.reduced$y[i], tol=1e-6, cor.conv=TRUE, maxiter = 2000))[1:9]  # %>% print
    }
  })
  param.sat <- t(param.sat)
# print(param.sat)  
  sat.lik   <- sum(dBvZINB4.vec(x= xy.reduced$x, y = xy.reduced$y, 
                                a0 = param.sat[,1], a1 = param.sat[,2], a2 = param.sat[,3], 
                                b1 = param.sat[,4], b2 = param.sat[,5],  
                                p1 = param.sat[,6], p2 = param.sat[,7], p3 = param.sat[,8], p4 = param.sat[,9], 
                                log = TRUE) * xy.reduced$freq)
  
  return(data.frame(model.likelihood = model.lik, satrtd.likelihood = sat.lik, deviance = 2*(sat.lik - model.lik), df = length(xvec) - 8))
}

if (FALSE) {
  dev.BvZINB4(extractor(15), extractor(20), model.lik = -160.85)
  ML.BvZINB4.2(extractor(15), extractor(20), tol=1e-6, cor.conv=TRUE, maxiter = 2000)  
}

