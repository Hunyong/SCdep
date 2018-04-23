## 0. library  
  # 0.1 data
  rm(list=ls())
  setwd("/Users/hycho/Documents/1STAT_Research/201704 SingleCellDep/Rcode1804")
  iset.Mm.c2 <- readRDS("iset.Mm.c2.rds")
  data <-read.delim("sox9ko_paneth_800_cell_counts.txt" , sep="\t")
  
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
  
## 1. Screening: Nonzero proportion > 10% ?
  if (FALSE) {
    # geneset 1
    screen.zero(data=data, geneset = iset.Mm.c2[[1]], cut=.1, output="stat")
    
    # First 100 genesets: For 10%, average 13% genes are above threshold
    1:100 %>% sapply(function(i) screen.zero(data=data, geneset = iset.Mm.c2[[i]], cut=.1, output="stat")) %>%
      apply(1, mean)
    
    # First 100 genesets: For 5%, average 29% genes are above threshold
    1:100 %>% sapply(function(i) screen.zero(data=data, geneset = iset.Mm.c2[[i]], cut=.05, output="stat")) %>%
      apply(1, mean)
  }
  
## 2. Fitting BvZINB with genes (> 5% nonzero)
  
  for (i in 1:30) { ## Fitting 1 ~ 30th genesets
    print(paste(i, "th geneset"))
    gene.tmp <- screen.zero(data=data, geneset = iset.Mm.c2[[i]], cut=.05, output="which")
    if (length(gene.tmp) <2) next
    data.tmp <- data[gene.tmp,]
    
    geneset <- pairwise.MLE(data = data.tmp, ML.fun = ML.BvZINB4.2, nondata.col = 1, showFlag=FALSE, tol=1e-5, cor.conv=TRUE, maxiter = 1000)
    geneset.list <- list(cutoff = "nonzero > 5%",
                         genes = gene.tmp, 
                         table = geneset)
    
    saveRDS(geneset.list, paste0("output/geneset.",i,".list.rds"))
  }
  tt(2)
  
  if (FALSE) {
     tt(1)
     ML.BvZINB4.2(data.tmp[1,-1] %>% as.numeric, data.tmp[2,-1]%>% as.numeric, maxiter=3000, showFlag=TRUE, tol=1e-5) -> a
     tt(2) # single pair: 38 minutes.
    # 2000th: 0.00283
    # 2500th: 0.00259
    # 2600th: 0.00256
    # 2700th: 0.00253
    # 2800th: 0.00250
    # 2900th: 0.00248
    # 3000th: 0.00246
    
     tt(1)
     ML.BvZINB4.2(data.tmp[1,-1] %>% as.numeric, data.tmp[5,-1]%>% as.numeric, maxiter=1000, showFlag=TRUE, tol=1e-5, cor.conv=TRUE) -> a
     tt(2) # single pair: 12 minutes.  (122 iter)
    # 0th: 0.56
    # 30th: 0.42
    # 100th: 0.395
    # 122th: 0.394  stop!
    # 400th: 0.3934
    # 430th: 0.3934
    
     tt(1)
     ML.BvZINB4.2(data.tmp[2,-1] %>% as.numeric, data.tmp[6,-1]%>% as.numeric, maxiter=1000, showFlag=TRUE, tol=1e-5, cor.conv=TRUE) -> a
     tt(2) # single pair: 1.3 minutes.  (352 iter)
    # 0th: 0.2
    # 50th: 0.106
    # 100th: 0.060
    # 200th: 0.043
    # 300th: 0.031
    
    
     tt(1)
     ML.BvZINB4.2(data.tmp[1,-1] %>% as.numeric, data.tmp[21,-1]%>% as.numeric, maxiter=1000, showFlag=TRUE, tol=1e-5, cor.conv=TRUE) -> a
     tt(2) # single pair: 3 minutes.  (92 iter)
    # 0th: 0.35
    # 50th: 0.202
    # 92th: 0.1819
    
     tt(1)
     ML.BvZINB4.2(data.tmp[3,-1] %>% as.numeric, data.tmp[16,-1]%>% as.numeric, maxiter=1000, showFlag=TRUE, tol=1e-5, cor.conv=TRUE) -> a
     tt(2) # single pair: 18sec.  (64 iter)
    # 0th: 0.09
    # 64th: 0.0823
  }

  
  