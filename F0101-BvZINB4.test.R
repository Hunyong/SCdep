library(tidyverse)

#### basis functions ####
  H0 <- function(k, m, a0, a1, a2, b1, b2, x, y) {
    gamma(a1 + k) * gamma(k+1)^-1 * gamma(a1)^-1 * gamma(x + y + a0 -m -k) * 
      gamma(x - k +1)^-1 * gamma(a0 + y - m)^-1 * 
      b1^x * b2^y * (b1 + b2 + 1)^(k + m - x - y -a0) *
      (b1 + 1)^-(k + a1) * (b2 + 1)^-(m + a2)}
  H1 <- function(k, a0, a1, a2, b1, b2, x, y) { # being lazy (y, a2, b2 not needed)
    gamma(a0 + x - k) * gamma(a0)^-1 * gamma(x - k +1)^-1 *
      gamma(a1 + k) * gamma(k + 1)^-1 * gamma(a1)^-1 * 
      b1^x * (b1 + 1)^-(x + a0 + a1)}
  H2 <- function(m, a0, a1, a2, b1, b2, x, y) {
    gamma(a0 + y - m) * gamma(a0)^-1 * gamma(y - m +1)^-1 *
      gamma(a2 + m) * gamma(m + 1)^-1 * gamma(a2)^-1 * 
      b2^x * (b2 + 1)^-(y + a0 + a2)}

#### functions for ElogR0 ####
  logR0.H0 <- function(k, m, a0, a1, a2, b1, b2, p1, p2, p3, p4, x, y) {
    H0(k, m, a0, a1, a2, b1, b2, x, y) *
      (digamma(x + y - k - m + a0) + log(b1/(1 + b1 + b2))) * p1
  }
   
  logR0.H1 <- function(k, a0, a1, a2, b1, b2, p1, p2, p3, p4, x, y) {
    if (y == 0) {
      H1(k, a0, a1, a2, b1, b2, x, y) *
        (digamma(x - k + a0) + log(b1/(1 + b1))) * p2 
    } else {0}
  }
  
  logR0.H2 <- function(m, a0, a1, a2, b1, b2, p1, p2, p3, p4, x, y) {
    if (x == 0) {
      H1(m, a0, a1, a2, b1, b2, x, y) *
        (digamma(y - m + a0) + log(b1/(1 + b2))) * p3 
    } else {0}
  }

#### functions for ElogR1 ####
  logR1.H0 <- function(k, m, a0, a1, a2, b1, b2, p1, p2, p3, p4, x, y) {
    H0(k, m, a0, a1, a2, b1, b2, x, y) *
      (digamma(k + a1) + log(b1/(1 + b1))) * p1
  }
  
  logR1.H1 <- function(k, a0, a1, a2, b1, b2, p1, p2, p3, p4, x, y) {
    if (y == 0) {
      H1(k, a0, a1, a2, b1, b2, x, y) *
        (digamma(k + a1) + log(b1/(1 + b1))) * p2 
    } else {0}
  }
  
  logR1.H2 <- function(m, a0, a1, a2, b1, b2, p1, p2, p3, p4, x, y) {
    if (x == 0) {
      H1(m, a0, a1, a2, b1, b2, x, y) *
        (digamma(a1) + log(b1)) * p3 
    } else {0}
  }

  
#### functions for ElogR2 ####
  logR2.H0 <- function(k, m, a0, a1, a2, b1, b2, p1, p2, p3, p4, x, y) {
    H0(k, m, a0, a1, a2, b1, b2, x, y) *
      (digamma(m + a2) + log(b1/(1 + b2))) * p1
  }
  
  logR2.H1 <- function(k, a0, a1, a2, b1, b2, p1, p2, p3, p4, x, y) {
    if (y == 0) {
      H1(k, a0, a1, a2, b1, b2, x, y) *
        (digamma(a2) + log(b1)) * p2 
    } else {0}
  }
  
  logR2.H2 <- function(m, a0, a1, a2, b1, b2, p1, p2, p3, p4, x, y) {
    if (x == 0) {
      H1(m, a0, a1, a2, b1, b2, x, y) *
        (digamma(m + a2) + log(b1/(1 + b2))) * p3 
    } else {0}
  }
  

dBvZINB4.Expt.test <- function(x, y, a0, a1, a2, b1, b2, p1, p2, p3, p4) {
  dBNB = dBvZINB4(x, y, a0, a1, a2, b1, b2, 1, 0, 0, 0, log=FALSE)
  dNB.x = dnbinom(x, a0 + a1, b1/(1+b1))
  dNB.y = dnbinom(y, a0 + a2, b2/(1+b2))
  dBZINB = dBvZINB4(x, y, a0, a1, a2, b1, b2, p1, p2, p3, p4, log=FALSE)
  dBZINB.0 = dBvZINB4(x, y, a0 + 1, a1, a2, b1, b2, p1, p2, p3, p4, log=FALSE)
  dBZINB.1 = dBvZINB4(x, y, a0, a1 + 1, a2, b1, b2, p1, p2, p3, p4, log=FALSE)
  dBZINB.2 = dBvZINB4(x, y, a0, a1, a2 + 1, b1, b2, p1, p2, p3, p4, log=FALSE)
  
  
  #### ER0 - ER2 ####
  ER0 = a0 * b1 * dBZINB.0 / dBZINB
  ER1 = a1 * b1 * dBZINB.1 / dBZINB
  ER2 = a2 * b1 * dBZINB.2 / dBZINB
  
  
  #### ER0 - ER2 ####
  EE1 = dBNB / dBZINB * p1
  EE2 = if (y == 0) {dNB.x / dBZINB * p2} else {0}
  EE3 = if (x == 0) {dNB.y / dBZINB * p3} else {0}
  EE4 = if (x + y == 0) {p4 / dBZINB} else {0}
  EX2 = ifelse(x + y > 0, y, 0) + (a0 + a2) * b1 /dBZINB * (ifelse(y == 0, dNB.x * p2, 0) + ifelse(x + y == 0, p4, 0))
  
  
  #### ElogR0 ####
  ElogR0.H0 = 
    sapply(0:y, function(m) {
      sapply(0:x, function(k) {
        logR0.H0(k, m, a0, a1, a2, b1, b2, p1, p2, p3, p4, x, y)
        }) %>% sum
    }) %>% sum
  
  ElogR0.H1 = 
    sapply(0:x, function(k) {
      logR0.H1(k, a0, a1, a2, b1, b2, p1, p2, p3, p4, x, y)
    }) %>% sum
  
  ElogR0.H2 = 
    sapply(0:y, function(m) {
      logR0.H2(m, a0, a1, a2, b1, b2, p1, p2, p3, p4, x, y)
    }) %>% sum
  
  ElogR0 = ElogR0.H0 + ElogR0.H1 + ElogR0.H2 + if (x + y == 0) {p4 *(digamma(a0) + log(b1))}
  

  #### ElogR1 ####
  ElogR1.H0 = 
    sapply(0:y, function(m) {
      sapply(0:x, function(k) {
        logR1.H0(k, m, a0, a1, a2, b1, b2, p1, p2, p3, p4, x, y)
      }) %>% sum
    }) %>% sum
  
  ElogR1.H1 = 
    sapply(0:x, function(k) {
      logR1.H1(k, a0, a1, a2, b1, b2, p1, p2, p3, p4, x, y)
    }) %>% sum
  
  ElogR1.H2 = 
    sapply(0:y, function(m) {
      logR1.H2(m, a0, a1, a2, b1, b2, p1, p2, p3, p4, x, y)
    }) %>% sum
  
  ElogR1 = ElogR1.H0 + ElogR1.H1 + ElogR1.H2 + if (x + y == 0) {p4 *(digamma(a1) + log(b1))}
  
  
  #### ElogR2 ####
  ElogR2.H0 = 
    sapply(0:y, function(m) {
      sapply(0:x, function(k) {
        logR2.H0(k, m, a0, a1, a2, b1, b2, p1, p2, p3, p4, x, y)
      }) %>% sum
    }) %>% sum
  
  ElogR2.H1 = 
    sapply(0:x, function(k) {
      logR2.H1(k, a0, a1, a2, b1, b2, p1, p2, p3, p4, x, y)
    }) %>% sum
  
  ElogR2.H2 = 
    sapply(0:y, function(m) {
      logR2.H2(m, a0, a1, a2, b1, b2, p1, p2, p3, p4, x, y)
    }) %>% sum
  
  ElogR2 = ElogR2.H0 + ElogR2.H1 + ElogR2.H2 + if (x + y == 0) {p4 *(digamma(a2) + log(b1))}
  if (length(ElogR0) == 0) {ElogR0 = NA}
  if (length(ElogR1) == 0) {ElogR1 = NA}
  if (length(ElogR2) == 0) {ElogR2 = NA}
  
  return(c(logdensity = log(dBZINB),
           ER0 = ER0, ER1 = ER1, ER2 = ER2, ElogR0 = ElogR0, ElogR1 = ElogR1, ElogR2 = ElogR2,
           EE1 = EE1, EE2 = EE2, EE3 = EE3, EE4 = EE4, EX2 = EX2))
}
dBvZINB4.Expt.test.vec <- Vectorize(dBvZINB4.Expt.test)
