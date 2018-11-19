source("Rcode1804/F0101-BvZINB4.test.R")
data <-read.delim("sox9ko_paneth_800_cell_counts.txt" , sep="\t")

tt(1)
a1 <- dBvZINB4.Expt.test.vec(data[1,-1], data[2,-1], 1, 1, 1, 1, 1, 0.6, 0.2, 0.1, 0.1)
tt(2)  #2.5 sec

tt(1)
a2 <- dBvZINB4.Expt.vec(data[1,-1], data[2,-1], 1, 1, 1, 1, 1, 0.6, 0.2, 0.1, 0.1)
tt(2)  #0.4 sec

cbind(a1[,1:5],a2[,1:5])
a = cbind(a1, a2)
a[, 1:800 * 2 - 1] = a1
a[, 1:800 * 2] = a2
colnames(a) = paste0(rep(1:800, each=2), c(".1", ".2"))

a[,1:10] # ElogR0~ ElogR2 differ, EX2 differs
a[,21:30]
which(abs(a[1,]  + 1.49165) > 0.1)
a[,which(abs(a[1,]  + 1.49165) > 0.1)]

