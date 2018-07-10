# 0.1 settings
  # setwd("/Users/hycho/Documents/1STAT_Research/201704 SingleCellDep/Rcode1804")
  rm(list=ls())
  iset.Mm.c2 <- readRDS("iset.Mm.c2.rds")
  data <-read.delim("sox9ko_paneth_800_cell_counts.txt" , sep="\t")
  
  args = commandArgs(trailingOnly=TRUE)  # passed from script
  # test if there is the argument: if not, return an error
  if (length(args)<2) stop("Two arguments must be supplied (the random set replicate, cutoff)", call.=FALSE)
  copy = args[1]
  cutoff = args[2]
  
# 0.2 library
  source("F0001-base.R")
  source("F0101-BvZINB4.R")
  source("F0101-BvZINB4-supp.R")
  source("F0101-BvZINB4-supp2.R")
  source("F0101-BvNB3.R")
  source("F0101-BvNB2.R")
  source("F0101-BvNB.R")
  source("F0101-BP.R")
  library(rootSolve)
  library(dplyr);library(magrittr)
  library(ggplot2); library(gridExtra)
  
  ML.BvZINB4.2.time = function (xvec, yvec, initial = NULL, tol=1e-8, maxiter=200, timelimit = 3600,
                                showFlag=FALSE, showPlot=FALSE, cor.conv = FALSE) {
    tryCatch({R.utils::withTimeout(ML.BvZINB4.2(xvec=xvec, yvec=yvec, initial = initial, tol=tol, maxiter=maxiter, 
                                      showFlag=showFlag, showPlot=showPlot, cor.conv = cor.conv), 
                         timeout=timelimit)}, TimeoutException = function(ex) {
                             cat("Timeout. Skipping.\n")
                             return(rep(NaN, 12)) # a0 a1 a2 b1 b2 p2 p3 p4 lik iter pureCor
                         })
    }
  
  
## 1. sampling & screening
  
  # sampling real gene sets, and calculate screened sizes
  set.seed(100)
  realGeneSet.index = sample(length(iset.Mm.c2), 100)  # randomly picked 100 real genesets no.
  realGeneSet = iset.Mm.c2[realGeneSet.index]          # list of 100 real genesets
  geneset.size.10 = sapply(realGeneSet, function(x) screened.length(data=data, geneset=x, cut = .10))
  geneset.size.05 = sapply(realGeneSet, function(x) screened.length(data=data, geneset=x, cut = .05))
  geneset.size.01 = sapply(realGeneSet, function(x) screened.length(data=data, geneset=x, cut = .01))
  geneset.size.02 = sapply(realGeneSet, function(x) screened.length(data=data, geneset=x, cut = .02))
  
  # screened data (random sampling frame)
  frame.01 = screen.zero (data, 1:dim(data)[1], cut = .01) %>% as.double
  frame.02 = screen.zero (data, 1:dim(data)[1], cut = .02) %>% as.double
  frame.05 = screen.zero (data, 1:dim(data)[1], cut = .05) %>% as.double
  frame.10 = screen.zero (data, 1:dim(data)[1], cut = .10) %>% as.double
  
  # generating random geneset with matched size
  set.seed(copy)  # for each simulation replicate, the random sets are different.
  randGeneSet.01 = lapply(geneset.size.01, function(x) sample(frame.01, x))
  randGeneSet.02 = lapply(geneset.size.02, function(x) sample(frame.02, x))
  
  set.seed(copy)  # for each simulation replicate, the random sets are different.
  randGeneSet.05 = lapply(geneset.size.05, function(x) sample(frame.05, x))
  randGeneSet.10 = lapply(geneset.size.10, function(x) sample(frame.10, x))
  
    # check lengths
    identical(sapply(randGeneSet.01, length), geneset.size.01)
    identical(sapply(randGeneSet.05, length), geneset.size.05)
    identical(sapply(randGeneSet.10, length), geneset.size.10)
    
## 2. fitting MLE
  
  
  for (i in 1:60) {
    
    if (geneset.size.05[i] < 2 | geneset.size.05[i] > 70) next
    
    data.tmp <- data[randGeneSet.05[[i]],]
    
    randGeneStat <- pairwise.MLE(data = data.tmp, ML.fun = ML.BvZINB4.2, nondata.col = 1, showFlag=FALSE, tol=1e-5, cor.conv=TRUE, maxiter = 1000)
    randGeneStat <- list(cutoff = ifelse(cutoff == 1, "nonzero > 5%","nonzero > 10%"),
                         genes = randGeneStat[[i]], 
                         table = randGeneStat)
    
    saveRDS(randGeneStat, paste0("output/Set.rand",copy,".",i, ifelse(cutoff==1, ".cut05", ".cut10"),".list.rds"))
  }
  

## 3. p-value distribution.
  randset.10.cor <- data.frame(matrix(NA,30,10)) #rows: sets, cols: simulation copies
  names(randset.cor) <- paste0("copy",2:11)
  randset.10.cor <- randset.cor
  moreThanTen = TRUE #count only pairs with nonzero more than 10%.
  
  for (copy in 9) {
  print(c(copy=copy))
    
    if (moreThanTen) {
      pdf(paste0("output/cor.plots-random",copy,".10.pdf"))
    } else {
      pdf(paste0("output/cor.plots-random",copy,".pdf"))
    }
    
    for (i in 1:30) {
      rand.set = readRDS(paste0("output/RandSet",copy,".index.rds"))
      if (length(screen.zero(data=data, geneset = rand.set[[i]], cut=.05, output="which"))<2) next
      print(paste("random ",i, "th geneset (copy #: ", copy,")"))
      a <- readRDS(paste0("output/RandSet", copy, ".", i,".list.rds"))
      # update names
      names(a$table)[4:15] <- c("a0","a1","a2", "b1", "b2", "p1", "p2", "p3", "p4", "lik","iter", "pureCor") 
      
      #augmenting - for symmetry
      if (moreThanTen) {
        a$table <- a$table[a$table$non0.min>=.1,]
      }
      if (dim(a$table)[1] == 0) next
      a2 <- a$table; a2[,c("1","2","non0.1","non0.2")] <- a2[,c("2","1","non0.2","non0.1")]
      a2 <- rbind(a$table, a2)
      
      
      
      cor.max = max(a2$pureCor)
      cor.mean <- round(mean(a2$pureCor, na.rm=TRUE),3)
      if (moreThanTen) {
        randset.10.cor[i, copy-1] <- cor.mean
      } else {
        randset.cor[i, copy-1] <- cor.mean
      }
      
      a2 %>% ggplot(aes(non0.1, non0.2, col=pureCor)) + geom_jitter() +
        scale_colour_gradient2() + ggtitle(paste0(i, "th random geneset: ", length(a$genes), " genes (nonzero > ", ifelse(moreThanTen, 10, 5), "%) (Copy# = ", copy, ")")) -> p1
      a$table %>% ggplot(aes(pureCor)) + geom_histogram() + 
        geom_density(aes(fill="red"), alpha=.4) + xlim(c(0, max(cor.max, 0.4))) +
        ggtitle(paste0(i, "th random geneset (Copy# = ", copy, ")")) + guides(fill=FALSE) +
        annotate("text", x= 0.3, y = 5, label = paste("Mean cor = ", cor.mean), col="blue") -> p2
      
      
      
      #if (k %% 10 == 0) {print(paste("printing", k, "th out of", length(gene.id)))}
      
      grid.arrange( p1, p2, ncol=1)
      
    } 
    
    dev.off()
  }  
  