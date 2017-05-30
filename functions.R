### making functions

## various types of ggplots
corMIplot <- function(data, group = as.formula(". ~ I(non0.min>.1) + I(non0.min>.05)"), type = "joint") {
  if (type == "joint") {x=data$MI; y= data$cor; title = "MI x Correlation by zero counts" }
  else if (type == "cor") {x=1; y= data$cor; title = "Correlation by zero counts"}
  else if (type == "MI") {x=1; y= data$MI; title = "MI by zero counts"}
  
  ggplot(data, aes(x=x, y=y, label=pair)) +
    facet_grid(group, labeller=labeller(.rows = label_both, .cols = label_value)) +
    geom_jitter(aes(color=1-non0.min), width = (type != "joint")*.4 ) +
    ggtitle(title)
}

## 1.1 pairwise association macro
pairwise.dep <- function(data, measure, name, nondata.col = 1, p.nonzero = TRUE, rm.id.vec = TRUE, log.data = FALSE, add.log = 1, p.value = FALSE, ...) {
  # data: a dataframe of which pairs are to be analyzed
  # measure can be either a function or a list of functions
  # nondata.col: the column number to be ignored for measurement eg. first col (gene names)
  # name: list of names of the association measures (need to find some code to extract from the list of functions)
  # log, add.log: log tansformation, addition before log transformation
  # p.nonzero: show proportion of non-zero values
  # p.value: Calculate the p-values by permutation test
  data <- as.matrix(data[,-nondata.col])
  if (log.data) {data <- log(data + add.log)}
  dim.p <- dim(data)[1]
  comb <- expand.grid(1:dim.p,1:dim.p)
  if (rm.id.vec)  {rw <- which(comb[,1] > comb[,2])} else {rw <- which(comb[,1] >= comb[,2])}
  comb <- data.frame(comb[rw,c(2,1)])
  # [lower.tri(matrix(1:dim.p^2,dim.p),TRUE),][,c(2,1)]
  # comb <- data.frame(t(combn(1:dim.p,2)))
  comb$pair <- apply(comb,1,function(x) paste(x,collapse="-"))
  
  # function of a measurement
  fun.measure <- function(measure0, comb0 = comb, data0 = data, p.value0 = p.value, ...) {
    a <- sapply(1:dim(comb0)[1], function(i) {
      x1 <- data0[comb0[i,1],]
      x2 <- data0[comb0[i,2],]
      point.est <- try(measure0(x1, x2, ...),silent=TRUE)
      point.est <- ifelse(class(point.est)=="try-error",NA, point.est)
      if (p.value0 == TRUE) {
        if (is.na(point.est) == FALSE) {
          p.val <- p.permutation(x1, x2, measure = measure0, point.est = point.est, ...)$p.val
        } else {p.val <- NA}
        return(data.frame(point.est=point.est, p.val=p.val))
      } else {return(data.frame(point.est=point.est, p.val=NA))}
      }
    )
   return(as.data.frame(t(a)))
  }
  if (class(measure) != "list") {measure <- list(measure)}
  ##sapply(measure, fun.measure)
  
  storage <- list(); j = 0
  for (x in measure) {
    j = j + 1
    storage[[j]] <- fun.measure(x)
    names(storage)[j] <- name[j]
  }
  for (k in 1:length(measure)) {comb <- cbind(comb, unlist(storage[[k]]$point.est))}
  for (k in 1:length(measure)) {comb <- cbind(comb, unlist(storage[[k]]$p.val))}
  names(comb) <- c(1,2,"pair", name, paste0(name,".p.val"))
  
  if (p.nonzero == TRUE) {
    n <- dim(data)[2]
    comb$non0.1 <- sapply(1:dim(comb)[1], function(i) {sum(data[comb[i,1],] != 0) / n})
    comb$non0.2 <- sapply(1:dim(comb)[1], function(i) {sum(data[comb[i,2],] != 0) / n})
    comb$non0.min <- pmin(comb$non0.1,comb$non0.2)
  }
  return(comb)
}
  # example
  pairwise.dep(matrix(1:100,5), mi.empirical.2, name = "mi", p.value=TRUE,n.permutation=20)
  pairwise.dep(matrix(1:100,5), cor, name = "cor")
  pairwise.dep(matrix(0,5,20), mi.empirical.2, name = "mi")  #all zero data -> all NA
  pairwise.dep(matrix(1:100,5), list(cor, mi.empirical.2), name = c("cor","MI"), log = TRUE)
  pairwise.dep(matrix(1:100,5), list(cor, mi.empirical.2), name = c("cor","MI"), p.value=TRUE, n.permutation=20)

# 1.1 mi.empirical.2: helper fn to make mi.empirical's arguments two vectors (not single matrix) to be applicable to pairwise.dep fn
# mi.empirical.2 <- function(x, y, bin=20, ...) {mi.empirical(discretize2d(x,y,numBins1=bin,numBins2=bin),...)}
mi.empirical.2 <- function(x, y, ...) {mi.empirical(table(x,y),...)}

# Aggregating
nBins <- function(x, type="Sturge") { # number of bins
  rng <- range(x)
  IQR <- IQR(x)
  SD <- sd(x)
  n <- length(x)
  if (type == "FD") {Dx <- rng*n^(1/3)/IQR/2}
  else if (type == "Scott") {Dx <- rng*n^(1/3)/SD/3.5}
  else if (type == "Sturge") {Dx <- 1 + log(n, base = 2)}
  return(max(2,round(Dx)))
}
mi.empirical.3 <- function(x, y, ...) {mi.empirical(discretize2d(x, y, numBins1 = nBins(x), numBins2 = nBins(y)),...)}
mi.empirical.3.SC <- function(x, y, ...) {mi.empirical(discretize2d(x, y, numBins1 = nBins(x, type="Scott"), numBins2 = nBins(y, type="Scott")),...)}

# log-scale aggregating
discretize2d.log <- function(x, y, log = TRUE, base=2, BinSize = 1, numBins1=nBins(x), numBins2=nBins(y), ...) {
  if (!log) {
    return (discretize2d(x,y,...))
  } else {
    x.log <- log(x+1, base)
    y.log <- log(y+1, base)
    r1 <- range(x.log)
    r2 <- range(y.log)
    if (is.null(BinSize)) { #different bin sizes (determined by numBin1,2)
      b1 <- seq(from = r1[1], to = r1[2], length.out = numBins1 +1)
      b2 <- seq(from = r2[1], to = r2[2], length.out = numBins2 +1)
    } else {
      b1 <- seq(from = r1[1], to = r1[2] + BinSize, by = BinSize)
      b2 <- seq(from = r2[1], to = r2[2] + BinSize, by = BinSize)
    }
    y2d = table(cut(x, breaks = base^b1-1, include.lowest = TRUE), 
                cut(y, breaks = base^b2-1, include.lowest = TRUE))
    return(y2d)
  }
}

mi.empirical.log <- function(x, y, ...) {mi.empirical(discretize2d.log(x,y,log=TRUE, BinSize = NULL, numBins1 = nBins(x), numBins2 = nBins(y)),...)}
mi.empirical.log2 <- function(x, y, ...) {mi.empirical(discretize2d.log(x,y,log=TRUE, BinSize = 1),...)}

  # example
  x <- as.numeric(data[1,-1]); y <- as.numeric(data[3,-1])
  discretize2d(x, y, numBins1 = nBins(x), numBins2 = nBins(y)) #plot5 (11 bins)
  discretize2d.log(x, y, log = TRUE, BinSize = NULL, numBins1 = nBins(x), numBins2 = nBins(y)) #plot5.2-1 (11 bins each, exponentiated)
  discretize2d.log(x, y, log = TRUE, BinSize = 1) #plot5.2-1 (equal size bins(2^k))


## p-value calculation (2-tail): for correlation, just standard p-value
p.permutation <- function(x, y, measure, point.est = NULL ,n.sample=1000, ...) {
  if (identical(measure, cor)) {     # if measure==cor, then no need for perm test but standard test
    a <- cor.test(x,y)
    return(data.frame(point.est = a$est, p.val = a$p.value))
  } else {
    if (is.null(point.est)) {point.est   <- measure(x,y,...)}
    replicates <- sapply(1:n.sample, function(a) {
      k <- sample(x); # print(table(k,y))
      l <- measure(k,y,...)
      # print(l)
      return(l)
      # measure(sample(x),y,...)
    })
    p.val <- mean(abs(replicates) > abs(point.est)) + mean(abs(replicates) == abs(point.est))/2
    return(data.frame(point.est=point.est, p.val=p.val))
  }
}
  # example
  p.permutation(1:5, 7:3, mi.empirical.2, n.sample=20) #p-val=0.5 always
  p.permutation(c(rep(0,10),1,1,2,2,3,3,3,3), c(1,7,3,rep(0,10),1,3,4,5,6), mi.empirical.2, n.sample=20) # p-val=0.25
  p.permutation(c(rep(0,10),1,1,2,2,3,3,3,3), c(1,7,3,rep(0,10),1,3,4,5,6), mi.empirical.2, n.sample=20) # p-val=0.25
  p.permutation(as.numeric(data[iset.Mm.c2[[1]][1],-1]),as.numeric(data[iset.Mm.c2[[1]][6],-1]), mi.empirical.2)


## matrix conversion: triangular elements to symmetric matrix
tri2mat <- function (trivec, diag = FALSE) {
  #diag: including the diagonal elements?
  d <- length(trivec)
  n <- round(sqrt(8*d+1)/2-1/2,0)
  if (n^2+n != 2*d) stop("the length is not n(n+1)/2.")
  if (!diag) {n <- n+1}
  X <- diag(NA, n)
  X[lower.tri(X, diag=diag)] <- trivec
  X <- X + t(X) - diag(diag(X))
  return(X)
}
# example
tri2mat(1:15)


# counts to repeated vectors
countsToCases <- function(x, countcol = "Freq") {
  # Get the row indices to pull from x
  idx <- rep.int(seq_len(nrow(x)), x[[countcol]])
  
  # Drop count column
  x[[countcol]] <- NULL
  
  # Get the rows from x
  x[idx, ]
}

# nonredundant result to duplicate result
 # (restoring the redundant pairs 1-2 => 1-2 & 2-1)
 # input = pairwise.dep function's output
redundant <- function(x) {
  tmp2 <- x
  tmp2$"1" <- x$"2"; tmp2$"2" <- x$"1"
  tmp2$non0.1 <- x$non0.2
  tmp2$non0.2 <- x$non0.1
  return(x <- rbind(x, tmp2))
}

# saving graphics
png2 <- function(image, filename = NULL, width = 720, height = 360,  units = "px"){
  if (is.null(filename)) {
    filename <- deparse(substitute(image))
    filename <- gsub("\\.","-",filename)
    filename <- paste0(filename,".png")
  }
  png(filename, width = width, height = height, units = "px")
  print(image)
  dev.off()
}

  # example
  png2(plot2.3.1)


## for 2.5 NMI
NMI <- function(MI) sqrt(1-exp(-2*MI))

NMI25 <- function(x, y, ...) {
  MI <- mi.empirical(table(x,y),...)
  HX <- entropy(table(x)); HY <- entropy(table(y))
  HXY <- HX + HY - MI
  # if (HXY < 0) {HXY <- 0}
  NMI2 <- 2* MI / (HX + HY)  ; if (!is.finite(NMI2)) NMI2 <- 0
  NMI3 <- MI / HXY        ; if (!is.finite(NMI3)) NMI3 <- 0
  NMI4 <- MI / min(HX, HY); if (!is.finite(NMI4)) NMI4 <- 0
  NMI5 <- MI / sqrt(HX*HY); if (!is.finite(NMI5)) NMI5 <- 0
  return(data.frame(NMI2 = NMI2, NMI3 = NMI3, NMI4 = NMI4, NMI5 = NMI5))
}
NMI2 <- function(x, y, ...) NMI25(x,y,...)$NMI2
NMI3 <- function(x, y, ...) NMI25(x,y,...)$NMI3
NMI4 <- function(x, y, ...) NMI25(x,y,...)$NMI4
NMI5 <- function(x, y, ...) NMI25(x,y,...)$NMI5
# computationally redundant, but done for the sake of convenience

  mi.empirical.2(as.numeric(data[iset.Mm.c2[[1]][7],-1]),as.numeric(data[iset.Mm.c2[[1]][9],-1])) #MI
  entropy(table(as.numeric(data[iset.Mm.c2[[1]][7],-1]),as.numeric(data[iset.Mm.c2[[1]][9],-1]))) #HXY
  entropy(table(as.numeric(data[iset.Mm.c2[[1]][9],-1]))) #HX
  NMI5(as.numeric(data[iset.Mm.c2[[1]][7],-1]),as.numeric(data[iset.Mm.c2[[1]][9],-1])) #NMI

## for 2.6 quant
polar <- function(cor, MI) {
  theta <- atan2(abs(cor),MI)
  norm <- sqrt(cor^2 + MI^2)
  return(data.frame(theta = theta, norm = norm))
}
polar2 <- function(data, a = "cor",  b = "MI", return.vec = FALSE) {
  q <- polar(data[,a], data[,b])
  m <- weighted.mean(q$theta, q$norm, na.rm = TRUE)
  v <- weighted.mean(q$theta^2, q$norm, na.rm = TRUE) - m^2
  rst <- data.frame(mean=m, sd = sqrt(v))
  if (return.vec) {rst <- list(stat = rst,vec = q)}
  return(rst)
}
  # example
  polar2(Mm.c2.micor.tmp[[1]])

polar3 <- function(data, a = "cor",  b = "MI") {
  q <- polar(data[,a], data[,b])
  data$theta <- q$theta
  data$norm <- q$norm
  return(data)
}

extractor <- function(gene.no , geneset.no = 1, dat=data, nondata.col = 1, genesetlist = iset.Mm.c2) {
  extracted <- dat[genesetlist[[geneset.no]][gene.no],-nondata.col]
  extracted <- as.numeric(extracted)
  return(extracted)
}
extractor(1)
