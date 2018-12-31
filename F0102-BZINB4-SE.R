### Final SE as of Dec 20, 2018

score.i <- function(x, y, a0, a1, a2, b1, b2, p1, p2, p3, p4) {
  expt <- dBvZINB4.Expt(x, y, a0, a1, a2, b1, b2, p1, p2, p3, p4, debug = FALSE)
  score <- 
    c(a0 = expt["log.R0.E"] - log(b1) - digamma(a0),
      a1 = expt["log.R1.E"] - log(b1) - digamma(a1),
      a2 = expt["log.R2.E"] - log(b1) - digamma(a2),
      b1 = - (expt["v.E"] + a0 + a1 + a2)/b1 + (expt["R0.E"] + expt["R2.E"]) * (b2 + 1)/b1^2 + 
        expt["R1.E"]/b1^2,
      b2 = expt["v.E"]/b2 - (expt["R0.E"] + expt["R2.E"]) /b1,
      p1 = expt["E1.E"]/p1 - expt["E4.E"]/p4,
      p2 = expt["E2.E"]/p2 - expt["E4.E"]/p4,
      p3 = expt["E3.E"]/p3 - expt["E4.E"]/p4)
  names(score) <- c(paste0("a", 0:2), paste0("b", 1:2), paste0("p", 1:3))
  score
}
score <- Vectorize(score.i, vectorize.args = c("x", "y"))

info <- function(x, y, a0, a1, a2, b1, b2, p1, p2, p3, p4, exact = FALSE) {
  s <- score(x, y, a0, a1, a2, b1, b2, p1, p2, p3, p4)
  s.sum <- apply(s, 1, sum)
s.tmp <<- s
  result <- s %*% t(s)
  if (exact) {result <- result - s.sum %*% t(s.sum)/length(x)}
  result
}
info(x = 3, y=4, a0 = 1, a1 = 2, a2 = 3, b1 = 1, b2 = 2, p1 = .3, p2 = .4, p3 = .1, p4 = .2)

BZINB4.se <- function(x, y, a0, a1, a2, b1, b2, p1, p2, p3, p4, ...) {
  if (any(!is.finite(c(a0,a1,a2,b1,b2,p1,p2,p3,p4)))) {return(rep(NA, 10))}
  info.mat <- info(x, y, a0, a1, a2, b1, b2, p1, p2, p3, p4, ...)
  cov.mat <- solve(info.mat)
  
  # variance of rho hat
  rho <- a0/sqrt((a0 + a1) * (a0 + a2)) *sqrt(b1 *b2 /(b1 + 1) /(b2 + 1))
  d.g <- rho * c(1/a0 - 1/{2*(a0 + a1)} - 1/{2*(a0 + a2)}, - 1/{2*(a0 + a1)}, - 1/{2*(a0 + a2)},
           1/{2 *b1 *(b1 + 1)}, 1/{2 *b2 *(b2 + 1)})
  var.rho <- t(d.g) %*% cov.mat[1:5, 1:5] %*% d.g
  
  # variance of p4 hat
  d.p4 <- -c(1, 1, 1)
  var.p4 <- t(d.p4) %*% cov.mat[6:8, 6:8] %*% d.p4
  
  se <- sqrt(c(diag(cov.mat), p4 = var.p4, rho=var.rho))
  se
}

# !!!this is complete info, for obs'd info mis|obs info should be adjusted.!!!
# info.c <- function(x, y, a0, a1, a2, b1, b2, p1, p2, p3, p4) {
#   expt <- dBvZINB4.Expt(x, y, a0, a1, a2, b1, b2, p1, p2, p3, p4, debug = FALSE)
#   a012 = -digamma(c(a0, a1, a2))
#   b1b1 = (expt["v.E"] + a0 + a1 + a2)/b1 - 2* (expt["R0.E"] + expt["R2.E"]) * (b2 + 1)/b1^3 - 
#     2 * expt["R1.E"]/b1^3
#   b2b2 = - expt["v.E"] /b1^2
#   b1b2 = (expt["R0.E"] + expt["R2.E"]) /b1^2
#   pp = - expt[c("E1.E", "E2.E", "E3.E")]/c(p1, p2, p3)^2
#   
#   hess.c = diag(c(a012, b1b1, b2b2, pp))
#   colnames(hess.c) <- rownames(hess.c) <- c(paste0("a", 0:2), paste0("b", 1:2), paste0("p", 1:3))
#   hess.c[c("a0", "a1", "a2"), "b1"] <- -1/b1
#   hess.c["b1", c("a0", "a1", "a2")] <- -1/b1
#   hess.c["b1", "b2"] <- hess.c["b2", "b1"] <- b1b2
#   
#   -hess.c
#   
# }
# info.c(x = 3, y=4, a0 = 1, a1 = 2, a2 = 3, b1 = 1, b2 = 2, p1 = .3, p2 = .4, p3 = .1, p4 = .2)

if (FALSE) {
  ### test!!
  set.seed(10)
  xy <- rBvZINB4(1000, 1,1,1,1,1,.5,.2,.2,.1)
  x = xy[,1] %>% as.numeric
  y = xy[,1] %>% as.numeric
  info(x, y, 1,1,1,1,1,.5,.2,.2,.1) %>% round(0)
  info.BvZINB4(x, y, 1,1,1,1,1,.5,.2,.2,.1) %>% round(0)
  info(x, y, 3.38,10e-06,10e-06,0.7583,0.7583,.5555,3.3e-133,2.25e-133,0.4445) %>% round(0)
  info.BvZINB4(x, y,  3.38,10e-06,10e-06,0.7583,0.7583,.5555,3.3e-133,2.25e-133,0.4445) %>% round(0)
  
}
