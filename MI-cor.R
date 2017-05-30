# relation between MI

set.seed(100)
a1 <- 1:1000; a2 <- a1 + rnorm(1000)*2
cor(a1,a2); mi.empirical.2(a1,a2) # 1, 2.8

a2 <- 1000 - a1 + rnorm(1000)*2
cor(a1,a2); mi.empirical.2(a1,a2) # -1, 2.85

a2 <- 1+rnorm(1000)
cor(a1,a2); mi.empirical.2(a1,a2) # 0, 0.2

a2 <- a1; mi.empirical.2(a1,a2)   # 1, 3
a2 <- c(1:500, 500:1)
cor(a1,a2); mi.empirical.2(a1,a2) # 0, 2.3

a1 <- c(rep(0,10),1,1,1,10,10)
a2 <- c(rep(0,10),1,1,1,10,10)
cor(a1,a2); mi.empirical.2(a1,a2) # 1, 0.9

a1 <- c(rep(0,10),1,1,1,10,10)
a2 <- c(10,10,1,1,1,rep(0,10))
cor(a1,a2); mi.empirical.2(a1,a2) # -0.2, 0.2

a1 <- c(rep(0,5),1,1,1,10,10)
a2 <- c(10,10,1,1,1,rep(0,5))
cor(a1,a2); mi.empirical.2(a1,a2) # -0.35, 0.7

a1 <- c(rep(0,5),1,1,1,10,100)
a2 <- c(100,10,1,1,1,rep(0,5))
cor(a1,a2); mi.empirical.2(a1,a2) # -0.14, 0.1

a1 <- c(0,1,0,1)
a2 <- c(1,0,1,0)
cor(a1,a2); mi.empirical.2(a1,a2) # -1, 0.7

a1 <- as.numeric(data[iset.Mm.c2[[1]][1],-1])
a2 <- as.numeric(data[iset.Mm.c2[[1]][2],-1])
bin=20;discretize2d(a1,a2,numBins1=bin,numBins2=bin)
table(a1,a2)
mi.empirical.2(a1,a2)   #0.01  double zero count: 762/800   # sum(a1==0 & a2==0)
mi.empirical.2(a1,a1)   #0.28  zero count: 767/800
mi.empirical.2(a2,a2)   #0.05  zero count: 794/800
mi.empirical(discretize2d(a1,a1,160,160))  #0.28

a3 <- a1[a1!=0| a2!=0]
a4 <- a2[a1!=0| a2!=0]
mi.empirical.2(a3,a4)  #0.51 double zero: 0/38
mi.empirical.2(a3,a3)  #2.6  zero count: 5/38
mi.empirical.2(a4,a4)  #0.6  zero count: 32/38

# empirical MI does not take into the order (and quantitative) information. It treats the counts as just categories.
# contingecy tables have the same MI
# # 0 1 2     # 4 3 7
# 0 10 . .    0 10 . .
# 3 1 1 2     1 1 1 2
# 5 . 5 1     4 . 5 1
a1 <- c(rep(0,10), 0,1,2,2,rep(1,5),2)
a2 <- c(rep(0,10), rep(3,4),rep(5,6))
mi.empirical.2(a1,a2)
a1 <- c(rep(4,10), 4,3,7,7,rep(3,5),7)
a2 <- c(rep(0,10), rep(1,4),rep(4,6))
mi.empirical.2(a1,a2)

a1 <- c(0,0,1,2,4)
a2 <- c(1,1,2,3,5)
mi.empirical.2(a1,a2)


