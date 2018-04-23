library(pracma)  # numerical integral
library(distrEx)   # Gauss-Legendre integration (GLIntegrate)
### Model fit: NB is way better than Poisson
library(MASS) # for glm.NB

if (FALSE) {
  a.nb <- sapply(1:50, function(s) {a <- try(glm.nb(I(extractor(s))~1)$deviance,silent=TRUE)
  return(if (class(a)=="try-error") NA else a)})
  a.poi <- sapply(1:50, function(s) glm(I(extractor(s))~1, family=poisson)$deviance)
  round(cbind(a.nb,a.poi),1)
  mean(a.nb,na.rm=TRUE);mean(a.poi)  
}



### Base function
bfunction.base <- function (w, x, y, a, b, b1, b2) { # using beta function the integration becomes very stable for small values of parameters
  return((1 + w * b1)^ -{a + b + x} * (1 + w * b2)^ -{a + b + y} * dbeta(w, a + x + y , b) * b1 ^ (a+b+x) * b2 ^ (a+b+y) * beta(a + x + y , b))
}
bfunction <- function (x, y, a, b, b1, b2, integrand = bfunction.base,  GL = FALSE, subdivisions =  2000) {
  if (!GL) {
    result <- try(integrate(integrand, lower = 0, upper = 1, x = x, y = y, a = a, b = b, b1 = b1, b2 = b2, subdivisions = subdivisions)$value)
  } else {
    fn1 <- function(w) {integrand (w, x, y, a, b, b1, b2)}
    result <- try(distrEx::GLIntegrate(f = fn1, lower = 0, upper = 1, order = subdivisions))
  }
  result <- if (class(result) == "try-error") {print(c(a,b,b1,b2)); return(0)} else result
  return(result)
}
# old functions with numerical problems
if (FALSE) {
  bfunction.base.old <- function (w, x, y, a, b, b1, b2) {
    return((1+1/w/b1)^{-x-a-b} * (1+1/w/b2)^{-y-a-b} * w^{-a-2*b-1} * (1-w)^{b-1})
  }
  bfunction.base.old2 <- function (w, x, y, a, b, b1, b2) {
    return((1 + w * b1)^ -{a + b + x} * (1 + w * b2)^ -{a + b + y} * w^{a + x + y - 1} * (1-w)^{b-1} * b1 ^ (a+b+x) * b2 ^ (a+b+y))
  }
  bfunction.base.log <- function (w, x, y, a, b, b1, b2) { # using beta function the integration becomes very stable for small values of parameters
  result <- -{a + b + x} * log (1 + w * b1) - {a + b + y} * log (1 + w * b2) + log (dbeta(w, a + x + y , b)) + (a+b+x) * log (b1) + (a+b+y) * log (b2) + lbeta(a + x + y , b)
  return(exp(result))
  
  # bfunction.GL uses Gauss-Legendre quadrature: alread incorporated in bfunction
  bfunction.GL <- function (x, y, a, b, b1, b2, integrand = bfunction.base, order = 1000) {
    fn1 <- function(w) {integrand (w, x, y, a, b, b1, b2)}
    result <- try(GLIntegrate(f = fn1, lower = 0, upper = 1, order = order))
    result <- if (class(result) == "try-error") 0 else result
    return(result)
  }
  
  
  }
  
  bfunction.new <- function (x, y, a, b, b1, b2, integrand = bfunction.base, reltol = 1e-8) {
    fn1 <- function(w) {integrand (w, x, y, a, b, b1, b2)}
    result <- try(integral(fn1, xmin = 0, xmax = 1, method="Kronrod", reltol = reltol))
    result <- if (class(result) == "try-error") 0 else result
  return(result)
}
}

bfunction.gr.a <- function (w, x, y, a, b, b1, b2) {
  bfunction.base(w, x, y, a, b, b1, b2) * (log(1+1/w/b1) + log(1+1/w/b2) + log(w))
}
bfunction.gr.b <- function (w, x, y, a, b, b1, b2, d=1) {
  bfunction.base(w, x, y, a, b, b1, b2) * (log(1+1/w/b1) + log(1+1/w/b2) + 2*log(w) - log(1-w))
}
bfunction.gr.aa <- function (w, x, y, a, b, b1, b2) {
  bfunction.base(w, x, y, a, b, b1, b2) * (log(1+1/w/b1) + log(1+1/w/b2) + log(w))^2
}
bfunction.gr.bb <- function (w, x, y, a, b, b1, b2) {
  bfunction.base(w, x, y, a, b, b1, b2) * (log(1+1/w/b1) + log(1+1/w/b2) + 2*log(w) - log(1-w))^2
}
bfunction.gr.ab <- function (w, x, y, a, b, b1, b2) {
  bfunction.base(w, x, y, a, b, b1, b2) * (log(1+1/w/b1) + log(1+1/w/b2) + log(w)) * (log(1+1/w/b1) + log(1+1/w/b2) + 2*log(w) - log(1-w))
}

if (FALSE) {
  bfunction(1,0,1,1,1,1)
  bfunction(0,1,1,1,1,1)
  bfunction(0,0,.01,1,2,1)
  integrate(bfunction.base.test, lower = 0, upper = 1, subdivisions = 200000, abs.tol = 1e-20)$value
  bfunction.new(0,0,.01,1,2,1)  
}



dBvNB1 <- function (x, y, a, b, b1, b2, log = FALSE, GL=FALSE, ...) {
  bfn <- bfunction(a = a, b = b, b1 = b1, b2 = b2, x = x, y = y, GL = GL, ...)
  result <- lgamma(x + a + b) + lgamma(y + a + b) - lgamma(a) - lgamma(b) - lgamma(a+b) - lgamma(x + 1) - lgamma(y + 1) - (a+b)*{log(b1) + log(b2)} + log(bfn) 
  if (is.infinite(result)) {result <- -Inf}  # for some large x, density = Inf (numerical imprecision)
  if (!log) {result <- exp(result)}
  return(result)
}
dBvNB1.vec <- Vectorize(dBvNB1, vectorize.args = c("x", "y"))
if (FALSE) {
  dBvNB1(0,1,1,1,1,1)
  dBvNB1(1,1,100,1,1,1)
  dBvNB1(0,0,.01,1,2,1)  
}



# Checking sum to 1 (Numerical imprecision for large x and y compared to small mean parameters)
if (FALSE) {
  a <- sapply(0:100, function(r) sapply(0:100, function(s) dBvNB1(s,r,100,10,2,1)))
  sum(a) - 1
  sum(a[1:50,1:50])
  sum(a[-1])
  a[80:100,80:100]
  a[1:10,1:10]
  persp(0:100,0:100, a)
}
# Convexity Check
if (FALSE) {
  # Checking convexity numerically (wrt a, b)
  a <- sapply(1:100, function(r) sapply(1:100, function(s) dBvNB1(1,1,s/10,r/10,1,1)))
  persp((1:100)/10,(1:100)/10, a)
  
  # Checking convexity numerically (wrt a, b)
  a <- sapply(1:100, function(r) sapply(1:100, function(s) dBvNB1(10,10,s/3,r/3,1,1)))
  persp((1:100)/10,(1:100)/10, a)
  
  # Checking convexity numerically (wrt betas)
  a <- sapply(1:100, function(r) sapply(1:100, function(s) dBvNB1(1,1,1,1,s/10,r/10)))
  persp((1:100)/10,(1:100)/10, a)
}

lik.BvNB1 <- function(x, y, param, a = NULL, b  = NULL, b1  = NULL, b2  = NULL, GL = FALSE){
  if (is.null(a)|is.null(b)|is.null(b1)|is.null(b2)) {
    a = param[1]; b = param[2]; b1 = param[3]; b2 = param[4]
  } 
  sum(dBvNB1.vec(x = x, y = y, a = a, b = b, b1 = b1, b2 = b2, GL=TRUE, log=TRUE))
}

gr.BvNB1.a <- function (x, y, a, b, b1, b2, GL=FALSE) { #d li / d a
  result <- digamma(a + b + x) + digamma(a + b + y) - digamma (a) - digamma (a + b) - log (b1) - log(b2) 
  result <- result - bfunction(x, y, a, b, b1, b2, bfunction.gr.a, GL=GL) / bfunction(x, y, a, b, b1, b2, bfunction.base, GL=GL)
  return(result)
}
gr.BvNB1.b <- function (x, y, a, b, b1, b2, GL=FALSE) { #d li / d a
  result <- digamma(a + b + x) + digamma(a + b + y) - digamma (b) - digamma (a + b) - log (b1) - log(b2) 
  result <- result - bfunction(x, y, a, b, b1, b2, bfunction.gr.b, GL=GL) / bfunction(x, y, a, b, b1, b2, bfunction.base, GL=GL)
  return(result)
}
gr.BvNB1.b1 <- function (x, y, a, b, b1, b2, GL=FALSE) { #d li / d a
  result <- (a + b + x) /b1^2 * bfunction(x, y-1, a+1, b, b1, b2, bfunction.base, GL=GL) / bfunction(x, y, a, b, b1, b2, bfunction.base, GL=GL)
  result <- result - (a + b) / b1
  return(result)
}
gr.BvNB1.b2 <- function (x, y, a, b, b1, b2, GL=FALSE) { #d li / d a
  result <- (a + b + y) /b2^2 * bfunction(x-1, y, a+1, b, b1, b2, bfunction.base, GL=GL) / bfunction(x, y, a, b, b1, b2, bfunction.base, GL=GL)
  result <- result - (a + b) / b2
  return(result)
}
gr.BvNB1.a.vec <- Vectorize(gr.BvNB1.a, vectorize.args = c("x", "y"))
gr.BvNB1.b.vec <- Vectorize(gr.BvNB1.b, vectorize.args = c("x", "y"))
gr.BvNB1.b1.vec <- Vectorize(gr.BvNB1.b1, vectorize.args = c("x", "y"))
gr.BvNB1.b2.vec <- Vectorize(gr.BvNB1.b2, vectorize.args = c("x", "y"))
gr.BvNB1.n <- function (x, y, a, b, b1, b2, GL=FALSE) {
  return(c(sum(gr.BvNB1.a.vec (x, y, a, b, b1, b2, GL=GL)),
           sum(gr.BvNB1.b.vec (x, y, a, b, b1, b2, GL=GL)),
           sum(gr.BvNB1.b1.vec (x, y, a, b, b1, b2, GL=GL)),
           sum(gr.BvNB1.b2.vec (x, y, a, b, b1, b2, GL=GL))))
}

gr.BvNB1.aa <- function (x, y, a, b, b1, b2, GL=FALSE) { #d^2 li / d a^2
  result <- trigamma(a + b + x) + trigamma(a + b + y) - trigamma (a) - trigamma (a + b)
  result <- result - bfunction(x, y, a, b, b1, b2, bfunction.gr.aa, GL=GL) / bfunction(x, y, a, b, b1, b2, bfunction.base, GL=GL)
  result <- result + (bfunction(x, y, a, b, b1, b2, bfunction.gr.a, GL=GL) / bfunction(x, y, a, b, b1, b2, bfunction.base, GL=GL))^2
  return(result)
}
gr.BvNB1.bb <- function (x, y, a, b, b1, b2, GL=FALSE) { #d^2 li / d a^2
  result <- trigamma(a + b + x) + trigamma(a + b + y) - trigamma (b) - trigamma (a + b)
  result <- result - bfunction(x, y, a, b, b1, b2, bfunction.gr.bb, GL=GL) / bfunction(x, y, a, b, b1, b2, bfunction.base, GL=GL)
  result <- result + (bfunction(x, y, a, b, b1, b2, bfunction.gr.b, GL=GL) / bfunction(x, y, a, b, b1, b2, bfunction.base, GL=GL))^2
  return(result)
}
gr.BvNB1.ab <- function (x, y, a, b, b1, b2, GL=FALSE) { #d^2 li / d a^2
  result <- trigamma(a + b + x) + trigamma(a + b + y) - trigamma (a + b)
  result <- result - bfunction(x, y, a, b, b1, b2, bfunction.gr.ab, GL=GL) / bfunction(x, y, a, b, b1, b2, bfunction.base, GL=GL)
  result <- result + bfunction(x, y, a, b, b1, b2, bfunction.gr.a, GL=GL) * bfunction(x, y, a, b, b1, b2, bfunction.gr.b, GL=GL) / bfunction(x, y, a, b, b1, b2, bfunction.base, GL=GL)^2
  return(result)
}
gr.BvNB1.b1b1 <- function (x, y, a, b, b1, b2, GL=FALSE) { #d^2 li / d beta1^2
  bf0 <- bfunction(x, y, a, b, b1, b2, bfunction.base, GL=GL)
  bf1 <- bfunction(x, y-1, a+1, b, b1, b2, bfunction.base, GL=GL)
  bf2 <- bfunction(x, y-2, a+2, b, b1, b2, bfunction.base, GL=GL)
  result <- -2*(a + b + x) /b1^3 * bf1 / bf0
  result <- result + (a + b + x)/b1^2 * (bf2/bf0 - bf1^2/bf0^2)
  result <- result - (a + b) / b1^2
  return(result)
}

hessian.BvNB1 <- function (x, y, a, b, b1, b2, GL=FALSE) { #d^2 li / d beta1^2
  bf <- bfunction(x, y, a, b, b1, b2, bfunction.base, GL=GL)
  bf.y1 <- bfunction(x, y-1, a+1, b, b1, b2, bfunction.base, GL=GL)
  bf.y2 <- bfunction(x, y-2, a+2, b, b1, b2, bfunction.base, GL=GL)
  bf.x1 <- bfunction(x-1, y, a+1, b, b1, b2, bfunction.base, GL=GL)
  bf.x2 <- bfunction(x-2, y, a+2, b, b1, b2, bfunction.base, GL=GL)
  bf.xy <- bfunction(x-1, y-1, a+2, b, b1, b2, bfunction.base, GL=GL)
  bf.aa <- bfunction(x, y, a, b, b1, b2, bfunction.gr.aa, GL=GL)
  bf.a  <- bfunction(x, y, a, b, b1, b2, bfunction.gr.a, GL=GL)
  bf.bb <- bfunction(x, y, a, b, b1, b2, bfunction.gr.bb, GL=GL)
  bf.b  <- bfunction(x, y, a, b, b1, b2, bfunction.gr.b, GL=GL)
  bf.ab <- bfunction(x, y, a, b, b1, b2, bfunction.gr.ab, GL=GL)
  bf.a.x1 <- bfunction(x-1, y, a+1, b, b1, b2, bfunction.gr.a, GL=GL)
  bf.a.y1 <- bfunction(x, y-1, a+1, b, b1, b2, bfunction.gr.a, GL=GL)
  bf.b.x1 <- bfunction(x-1, y, a+1, b, b1, b2, bfunction.gr.b, GL=GL)
  bf.b.y1 <- bfunction(x, y-1, a+1, b, b1, b2, bfunction.gr.b, GL=GL)
  tri.gamma.ab <- trigamma(a + b + x) + trigamma(a + b + y) - trigamma (a + b)
  
  # hessians for a and b
  aa <- tri.gamma.ab - trigamma (a) - bf.aa / bf + (bf.a / bf)^2
  bb <- tri.gamma.ab - trigamma (b) - bf.bb / bf + (bf.b / bf)^2
  ab <- tri.gamma.ab - bf.ab / bf + (bf.a * bf.b / bf^2)

  # hessians for beta's
  b1b1 <- -2*(a + b + x) /b1^3 * bf.y1 / bf
  b1b1 <- b1b1 + (a + b + x)/b1^2 * (bf.y2/bf - bf.y1^2/bf^2)
  b1b1 <- b1b1 + (a + b) / b1^2

  b2b2 <- -2*(a + b + y) /b2^3 * bf.x1 / bf
  b2b2 <- b2b2 + (a + b + y)/b2^2 * (bf.x2/bf - bf.x1^2/bf^2)
  b2b2 <- b2b2 + (a + b) / b2^2
  
  b1b2 <- (a + b + x) * (a + b + y) /b1^2 /b2^2 * (bf.xy/bf - bf.x1 *bf.y1 /bf^2)
  
  # hessians for beta - ab
  a.b1 <- - 1/b1 + 1/b1^2 *bf.y1/bf - (a+b+x)/b1^2 * (bf.a.y1/bf - bf.y1*bf.a/bf^2)
  a.b2 <- - 1/b2 + 1/b2^2 *bf.x1/bf - (a+b+y)/b2^2 * (bf.a.x1/bf - bf.x1*bf.a/bf^2)
  b.b1 <- - 1/b1 + 1/b1^2 *bf.y1/bf - (a+b+x)/b1^2 * (bf.b.y1/bf - bf.y1*bf.b/bf^2)
  b.b2 <- - 1/b2 + 1/b2^2 *bf.x1/bf - (a+b+y)/b2^2 * (bf.b.x1/bf - bf.x1*bf.b/bf^2)
  
  result <- matrix(c(aa, ab, a.b1, a.b2, NA, bb, b.b1, b.b2, NA, NA, b1b1, b1b2, NA, NA, NA, b2b2),4,4)
  result[upper.tri(result)] = t(result)[upper.tri(result)]
  return(result)
}
hessian.BvNB1.vec <- Vectorize(hessian.BvNB1, vectorize.args = c("x", "y"))

ML.BvNB1 <- function (x, y, GL=FALSE, tol=1e-6) {
  xy.reduced <- as.data.frame(table(x,y))
  names(xy.reduced) <- c("x", "y","freq")
  xy.reduced <- xy.reduced[xy.reduced$freq != 0,]
  xy.reduced$x <- as.numeric(as.character(xy.reduced$x))
  xy.reduced$y <- as.numeric(as.character(xy.reduced$y))
  xy.reduced$freq <- as.numeric(as.character(xy.reduced$freq))
  # print(xy.reduced)
  
  fn.1 = function (param) {
    lik <- dBvNB1( x = xy.reduced$x, y = xy.reduced$y,  a = param[1], b = param[2], b1 = param[3], b2 = param[4], GL=GL, log=TRUE)
    result <- lik %*% xy.reduced$freq
    #print(param)
    #print(result)
    return(result)
  }
  gr.1 = function (param) {
    a <- gr.BvNB1.a.vec ( x = xy.reduced$x, y = xy.reduced$y, a = param[1], b = param[2], b1 = param[3], b2 = param[4], GL=GL)
    b <- gr.BvNB1.b.vec ( x = xy.reduced$x, y = xy.reduced$y, a = param[1], b = param[2], b1 = param[3], b2 = param[4], GL=GL)
    b1 <- gr.BvNB1.b1.vec ( x = xy.reduced$x, y = xy.reduced$y, a = param[1], b = param[2], b1 = param[3], b2 = param[4], GL=GL)
    b2 <- gr.BvNB1.b2.vec ( x = xy.reduced$x, y = xy.reduced$y, a = param[1], b = param[2], b1 = param[3], b2 = param[4], GL=GL)
    result <- t(cbind(a,b,b1,b2)) %*% xy.reduced$freq
    #print(param)
    #print(result)
    return(result)
  }
  hs.1 = function (param) {
    hessian <- hessian.BvNB1.vec(x = xy.reduced$x, y = xy.reduced$y,  a = param[1], b = param[2], b1 = param[3], b2 = param[4], GL=GL)
    hessian <- hessian %*% xy.reduced$freq
    hessian <- matrix(hessian, 4, 4)
    return(hessian)
  }
  
  #log-scaled params: param.l
  fn.log = function (param.l) { fn.1 (exp(param.l)) }
  gr.log = function (param.l) { gr.1 (exp(param.l)) * exp(param.l) }
  hs.log = function (param.l) { 
    hs.log <- diag(exp(param.l)) %*% hs.1 (exp(param.l)) %*% diag(exp(param.l)) 
    hs.log <- diag(gr.1 (exp(param.l)) * exp(param.l)) + hs.log 
    return(hs.log)
    }
  
  #print(hs.1(c(1,1,1,1)))
  #print(hs.log(c(1,1,1,1)))
  
  # Newton-Raphson with exponentiated parameters
  if (TRUE) {
    iter = 0
    param.l = rep(0,4)
    repeat {
      param.old <- param.l # saving old parameters
      iter = iter + 1
      fn <- fn.log(param.l) # %>% print
      gr <- gr.log(param.l) %>% print
      hs <- hs.log(param.l) %>% print
      param.l <- as.vector(param.l - solve(hs) %*% gr)
      print(param.l)
      print(c("iter:",iter, round(exp(param.l),4), "lik: ",round(fn,2)))
      
      if (sum(abs(param.l - param.old)) <= tol) {
        return(exp(param.l))
        break
      } else if (iter >= 200) {
        return(NA)
        break
      }
    }
  }
  
  # require(rootSolve)
  # result <- multiroot(f = gr.log, start = c(1,1,1,1), positive=TRUE)
  # result <- optim(par = c(1,1,1,1), fn = fn.1, gr = gr.1, lower=rep(1e-20,4), control=list(fnscale=-1), method = "L-BFGS-B")$par
  # result <- optim(par = c(1,1,1,1), fn = fn.log, control=list(fnscale=-1), method = "BFGS")
  return(result)
}

if (FALSE) {
  tt(1)
  ML.BvNB1(extractor(11)[1:30], extractor(16)[1:30], GL=TRUE)
  tt(2) #14 mins: 136 iterations
  # 0.144747 0.231625 1.560235 3.167277 #lik=-51.19
  
  tt(1)
  ML.BvNB1(extractor(11), extractor(16), GL=TRUE)
  tt(2) # 3.4 hours after stopping with 200 iterations
  # "0.1265"   "0.2143"   "20.9506"  "0.9819"   "lik: "    "-2153.51" -> Not optimized
  
  rm(GL, xy.reduced, hs.1)
  tt(1)
  lik.BvNB1(extractor(1)[1:30], extractor(2)[1:30], param = c(.01,3.45,1.25,1.009), GL=TRUE)
  tt(2)
  gr.BvNB1.n(extractor(1), extractor(2), 1,1,1,1)
  lik.BvNB1(extractor(1), extractor(2), param = c(1,1,1,1))
  gr.BvNB1.n(extractor(1), extractor(2), 1,1,1,1)
  bfunction(1,1,2,1,1,1)
  
}
