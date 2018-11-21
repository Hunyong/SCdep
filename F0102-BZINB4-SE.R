# BvZINB4: getting standard errors of the estimator.
rm(list=ls())
# x = 3; y = 5; a0 = 1; a1 = 2; a2 = 3; b1 = 1; b2 = 2; p1 = .3; p2 = .4; p3 = .1; p4 = .2

info.i.BvZINB4 <- function(x, y, a0, a1, a2, b1, b2, p1, p2, p3, p4, param = NULL) {
  if (!is.null(param)) {
    a0 = param[1]; a1 = param[2]; a2 = param[3]
    b1 = param[4]; b2 = param[5]
    p1 = param[6]; p2 = param[7]; p3 = param[8]; p4 = param[9]
  }
  ### basic materials
    t1 = (b1 + b2 + 1) /(b1 + 1); t2 = (b1 + b2 + 1) /(b2 + 1)
    adj = 0
    l1 <- function(k, m) exp(lgamma(a1 + k) - lgamma(k+1) - lgamma(a1) + lgamma(x + y + a0 -m -k) - lgamma(x -k +1) - lgamma(a0 + y - m) 
                             + lgamma(m + a2) - lgamma(m+1) - lgamma(a2) + lgamma(y +a0 -m) - lgamma(y -m +1) - lgamma(a0))
    l1.B <- - (+x+y+a0)*log(1 + b1 + b2) + x * log(b1) + y * log(b2) - a1 * log(1 + b1) - a2 * log(1 + b2)
    l1.B <- exp(l1.B)  #pi1 not multiplied here. (!!! F0101-BZINB file has pi1 multiplied !!!)
    
    l1.C <- function(k, m) exp(k *log(t1) + m *log(t2))
    
    
    l.A.mat <- sapply(0:x, function(k) sapply(0:y, l1, k=k))  # %>% print
    l.C.mat <- sapply(0:x, function(k) sapply(0:y, l1.C, k=k))  # %>% print
    if (is.infinite(sum( l.A.mat))) {
      adj = 200
      l.A.mat <- sapply(0:x, function(k) sapply(0:y, function(m) {l1(k =k, m = m) *exp(-adj)}))
    } #%>%print
    #adjustment is cancelled out for each Expectation, so can be ignored. But for the final likelihood it should be adjusted at the end.
    
    log.R0.E1 <- function(k, m) {digamma(x - k + y - m + a0)}
    log.R1.E1 <- function(k) {digamma(k + a1)}
    log.R2.E1 <- function(m) {digamma(m + a2)}
    R0.E1 <- function(k, m) {x - k + y - m + a0}
    R1.E1 <- function(k) {k + a1}
    R2.E1 <- function(m) {m + a2}
    
    log.R0.mat <- sapply(0:x, function(k) sapply(0:y, log.R0.E1, k=k))
    log.R1.mat <- t(matrix(sapply(0:x, log.R1.E1), x+1, y+1))
    log.R2.mat <- matrix(sapply(0:y, log.R2.E1), y+1, x+1)
    
    R0.mat <- sapply (0:x, function(k) sapply(0:y, R0.E1, k=k))/(b1 + b2 + 1)
    R1.mat <- t(matrix(sapply(0:x, R1.E1), x+1, y+1))/(b1 + 1)
    R2.mat <- matrix(sapply(0:y, R2.E1), y+1, x+1)/(b2 + 1)
    
    
    digamma.a <- digamma(c(a0, a1, a2))
    log.b <- log(c(b1 + b2 + 1, b1 + 1, b2 + 1))
    dga.logb = digamma.a + log.b
    
    
    l.AC.mat = l.A.mat * l.C.mat
    sum.AC = sum(l.AC.mat)   # sum.AC * l1.B = sum(sum((H0))
    
    dNB = function(x, a, b) exp(lgamma(x + a) - lgamma(x + 1) - lgamma(a) + x * log(b) - (x + a) * log(b + 1))
    dNB.x = dNB(x, a0 + a1, b1)
    dNB.y = dNB(y, a0 + a2, b2)
    dBNB.xy = l1.B * sum(l.AC.mat)
    zx = ifelse(x == 0, 1, 0)
    zy = ifelse(y == 0, 1, 0)
    dBvZINB.xy = 
      dBNB.xy * p1 + dNB.x * zy * p2 + 
      dNB.y * zx * p3 + zx * zy * p4
    # cat("d.xy=", dBvZINB.xy, "d.x = ", dNB.x, "d.y =", dNB.y)    
    
  ### score function
    D1 = l1.B * (sum(l.AC.mat * log.R0.mat) - sum.AC * dga.logb[1])
    D2 = l1.B * (sum(l.AC.mat * log.R1.mat) - sum.AC * dga.logb[2])
    D3 = l1.B * (sum(l.AC.mat * log.R2.mat) - sum.AC * dga.logb[3])
    D4 = l1.B * (sum(l.AC.mat * (- R0.mat - R1.mat)) + sum.AC * (x/b1))
    D5 = l1.B * (sum(l.AC.mat * (- R0.mat - R2.mat)) + sum.AC * (y/b2))
    
    D6 = dNB.x * (digamma(a0 + a1 + x) - digamma(a0 + a1) - log(b1 + 1))
    D7 = dNB.x * (x/b1 - (a0 + a1 + x)/(b1 + 1))
    D8 = dNB.y * (digamma(a0 + a2 + y) - digamma(a0 + a2) - log(b2 + 1))
    D9 = dNB.y * (y/b2 - (a0 + a2 + y)/(b2 + 1))
    
    score = c(p1 * D1 + p2 * zy * D6 + p3 * zx * D8,
              p1 * D2 + p2 * zy * D6,
              p1 * D3 + p3 * zx * D8,
              p1 * D4 + p2 * zy * D7,
              p1 * D5 + p3 * zx * D9,
              dBNB.xy - zx * zy,
              zy * dNB.x - zx * zy,
              zx * dNB.y - zx * zy)
    
    score = score / dBvZINB.xy
    
  ### Materials for information
    tri12 <- function(k, m) {trigamma(x - k + y - m + a0)}
    tri12.mat <- sapply(0:x, function(k) sapply(0:y, tri12, k=k))
    tri1 <- function(k) {trigamma(k + a1)}
    tri1.mat <- t(matrix(sapply(0:x, tri1), x+1, y+1))
    tri2 <- function(m) {trigamma(m + a2)}
    tri2.mat <- matrix(sapply(0:y, tri2), y+1, x+1)
    tra = trigamma(c(a0, a1, a2))
    
  ### information
    A1 = (log.R0.mat - dga.logb[1])^2 + (tri12.mat - tra[1])
    A1 = l1.B * sum(l.AC.mat * A1)
    
    A2 = (log.R0.mat - dga.logb[1]) * (log.R1.mat - dga.logb[2])
    A2 = l1.B * sum(l.AC.mat * A2)
    
    A3 = (log.R0.mat - dga.logb[1]) * (log.R2.mat - dga.logb[3])
    A3 = l1.B * sum(l.AC.mat * A3)
    
    A4 = (log.R0.mat - dga.logb[1]) * (- R0.mat - R1.mat + x/b1) - 1/(b1 + b2 + 1)
    A4 = l1.B * sum(l.AC.mat * A4)
    
    A5 = (log.R0.mat - dga.logb[1]) * (- R0.mat - R2.mat + y/b2) - 1/(b1 + b2 + 1)
    A5 = l1.B * sum(l.AC.mat * A5)
    
    A6 = (log.R1.mat - dga.logb[2])^2 + (tri1.mat - tra[2])
    A6 = l1.B * sum(l.AC.mat * A6)
  
    A7 = (log.R1.mat - dga.logb[2]) * (log.R2.mat - dga.logb[3])
    A7 = l1.B * sum(l.AC.mat * A7)
    
    A8 = (log.R1.mat - dga.logb[2]) * (- R0.mat - R1.mat + x/b1) - 1/(b1 + 1)
    A8 = l1.B * sum(l.AC.mat * A8)
    
    A9 = (log.R1.mat - dga.logb[2]) * (- R0.mat - R2.mat + y/b2)
    A9 = l1.B * sum(l.AC.mat * A9)
    
    A10 = (log.R2.mat - dga.logb[3])^2 + (tri2.mat - tra[3])
    A10 = l1.B * sum(l.AC.mat * A10)
    
    A11 = (log.R2.mat - dga.logb[3]) * (- R0.mat - R1.mat + x/b1)
    A11 = l1.B * sum(l.AC.mat * A11)
    
    A12 = (log.R2.mat - dga.logb[3]) * (- R0.mat - R2.mat + y/b2) - 1/(b2 + 1)
    A12 = l1.B * sum(l.AC.mat * A12)
    
    A13 = (- R0.mat - R1.mat + x/b1)^2 - (- R0.mat/(b1 + b2 + 1) - R1.mat/(b1 + 1) + x/b1^2)
    A13 = l1.B * sum(l.AC.mat * A13)
    
    A14 = (- R0.mat - R1.mat + x/b1) * (- R0.mat - R2.mat + y/b2) - (- R0.mat/(b1 + b2 + 1))
    A14 = l1.B * sum(l.AC.mat * A14)
    
    A15 = (- R0.mat - R2.mat + y/b2)^2 - (- R0.mat/(b1 + b2 + 1) - R2.mat/(b2 + 1) + y/b2^2)
    A15 = l1.B * sum(l.AC.mat * A15)
    
    B1 = (digamma(a0 + a1 + x) - digamma(a0 + a1) - log(b1 + 1))^2 + 
      (trigamma(a0 + a1 + x) - trigamma(a0 + a1))
    B1 = dNB.x * B1  
    
    B2 = (digamma(a0 + a1 + x) - digamma(a0 + a1) - log(b1 + 1)) *
      (x/b1 - (a0 + a1 + x)/(b1 + 1)) - 1/(b1 + 1)
    B2 = dNB.x * B2
    
    B3 = (x/b1 - (a0 + a1 + x)/(b1 + 1))^2 - (x/b1^2 - (a0 + a1 + x)/(b1 + 1)^2)
    B3 = dNB.x * B3
    
    C1 = (digamma(a0 + a2 + y) - digamma(a0 + a2) - log(b2 + 1))^2 + 
      (trigamma(a0 + a2 + y) - trigamma(a0 + a2))
    C1 = dNB.y * C1  
    
    C2 = (digamma(a0 + a2 + y) - digamma(a0 + a2) - log(b2 + 1)) *
      (y/b2 - (a0 + a2 + y)/(b2 + 1)) - 1/(b2 + 1)
    C2 = dNB.y * C2
    
    C3 = (y/b2 - (a0 + a2 + y)/(b2 + 1))^2 - (y/b2^2 - (a0 + a2 + y)/(b2 + 1)^2)
    C3 = dNB.y * C3
  
  ### result  
    info = matrix(NA, 8, 8)
    info[lower.tri(info, diag = TRUE)] =
      c(p1 * A1 + p2 * zy * B1 + p3 * zx * C1,  # 1st column
        p1 * A2 + p2 * zy * B1,
        p1 * A3 + p3 * zx * C1,
        p1 * A4 + p2 * zy * B2,
        p1 * A5 + p3 * zx * C2,
        D1, D6 * zy, D8 * zx,
        p1 * A6 + p2 * zy * B1,            # 2nd column
        p1 * A7,
        p1 * A8 + p2 * zy * B2,
        p1 * A9,
        D2, D6 * zy, 0,
        p1 * A10 + p3 * zx * C1,           # 3rd column
        p1 * A11,
        p1 * A12 + p3 * zx * C2,
        D3, 0, D8 * zx,
        p1 * A13 + p2 * zy * B3,           # 4th column
        p1 * A14,
        D4, D7 * zy, 0,
        p1 * A15 + p3 * zx * C3,           # 5th column
        D5, 0, D9 * zx,
        0, 0, 0, 0, 0, 0)             # 6th - 8th columns
    
    # fill in off-diagonals
    info[upper.tri(info)] = t(info)[upper.tri(t(info))]
    
    info = info / dBvZINB.xy  - score %*% t(score)
    info = - info
      
    return(list (score = score, info = info))
}
info.i.BvZINB4.vec <- Vectorize(info.i.BvZINB4, vectorize.args = c("x", "y"))

if (FALSE) {
  info.i.BvZINB4.vec(x = 3, y = 5, a0 = 1, a1 = 2, a2 = 3, b1 = 1, b2 = 2, p1 = .3, p2 = .4, p3 = .1, p4 = .2)[[1]]
  info.i.BvZINB4.vec(x = 3, y = 5, a0 = 1, a1 = 2, a2 = 3, b1 = 1, b2 = 2, p1 = .3, p2 = .4, p3 = .1, p4 = .2)[[2]] %>% diag
  info.i.BvZINB4.vec(x = 4, y = 10, a0 = 1, a1 = 2, a2 = 3, b1 = 1, b2 = 2, p1 = .3, p2 = .4, p3 = .1, p4 = .2)[[2]] %>% diag
  ML.BvZINB4.2b(x = c(3,4,5), y = c(5,10,2))
  info.i.BvZINB4.vec(x = c(3,4,5), y = c(5,10,2), a0 = 1.27, a1 = 1.407e+4, a2 = 4.2423, b1 = 2.84e-04, 
               b2 = 1.028, p1 = 1, p2 = 0, p3 = 0, p4 = 0)[1,]
}

info.BvZINB4 <- function(x, y, ...) {
  result = info.i.BvZINB4.vec(x, y, ...)
  n = length(result[1,])
  score = Reduce("+", result[1,])
  info = Reduce("+", result[2,]) / n
print(diag(info))
  vcov = solve(info)
print(diag(vcov))
  info
  #result
}
if (FALSE) {
  info.BvZINB4(x = c(3,4,7,0), y = c(4,10,3,0), a0 = 1, a1 = 2, a2 = 3, b1 = 1, b2 = 2, p1 = .3, p2 = .4, p3 = .1, p4 = .2)
  ML.BvZINB4.2b(c(3,4,7,0), c(4,10,3,0), boosting=FALSE)
  ML.BvZINB4.2b(c(3,4,7), c(4,10,3), maxiter=1000, boosting=FALSE)
  info.BvZINB4(x = c(3,4,7), y = c(4,10,3), a0 = 0.556, a1 = 4622, a2 = 8.32, b1 = 1e-3, b2 = 0.6384, p1 = 1, p2 = 0, p3 = 0, p4 = 0)[[1]]
  info.BvZINB4(x = c(3,4,7,0), y = c(4,10,3,0), a0 = 1, a1 = 7.4286, a2 = 9.405e-04, b1 = 2.403e-02, b2 = 3.034e-02, p1 = 0.75003, p2 = 0, p3 = 0, p4 = 2.49968e-01)[[2]] %>% solve
}


