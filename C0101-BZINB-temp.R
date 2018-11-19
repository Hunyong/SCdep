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
library(ggplot2); library(gridExtra)


for (i in 1:30) { ## Fitting 1 ~ 30th genesets
  print(paste(i, "th geneset"))
  gene.tmp <- screen.zero(data=data, geneset = iset.Mm.c2[[i]], cut=.05, output="which")
  if (length(gene.tmp) <2) next
  data.tmp <- data[gene.tmp,]
  # read previous output
### change number here!!!
  a <- readRDS(paste0("output/geneset.", i, ".list.rds"))
  
  if (TRUE) {
    # update names
    names(a$table)[4:15] <- c("a0","a1","a2", "b1", "b2", "p1", "p2", "p3", "p4", "lik","iter", "pureCor") 
  
    # recalculate the pureCor, and add se, pval.
    tmp.pureCor = apply(a$table[1:2,], 1, function(x) {
print(x)
      if (any(is.na(x[4:12]))) {rep(NA,3)} else {
        pureCor (as.numeric(x[4:12]), data.tmp[as.numeric(x[1]), -1], data.tmp[as.numeric(x[2]), -1])  
      }})
    a$pureCor = t(tmp.pureCor)[,1]
    a$pureCor.se = t(tmp.pureCor)[,2]
    a$pureCor.p = t(tmp.pureCor)[,3]
  } 
  print(c(i=i, line="line40"))
  # add deviance
  if (FALSE) {
    dev = apply(a$table, 1, function(x) {
      print(x)
      dev.BvZINB4(  as.numeric(data.tmp[as.numeric(x[1]), -1]), 
                    as.numeric(data.tmp[as.numeric(x[2]), -1]),
                    model.lik = as.numeric(x[13]))})
    
    dev = do.call(rbind, dev)  # convert a list to a stacked-up data.frame
    a$table = cbind(a$table, dev[,-1])
    
  
  saveRDS(a, paste0("output/Pgeneset.", i, ".list.rds"))
  }
}


### a list of objects
    geneset.cor.list <- lapply(1:30, function(i) { ## Fitting 1 ~ 30th genesets
      print(paste(i, "th geneset"))
      if (length(screen.zero(data=data, geneset = iset.Mm.c2[[i]], cut=.05, output="which"))<2) {
        NULL
      } else {
        a <- readRDS(paste0("output/geneset.",i,".list.rds"))
        
        # update names
        names(a$table)[4:15] <- c("a0","a1","a2", "b1", "b2", "p1", "p2", "p3", "p4", "lik","iter", "pureCor") 
        a    
      }
    })

    # relation btw size and mean cor
    sapply(geneset.cor.list, function(x) {
      if (is.null(x)) {c(NA,NA,NA)} else {
        c(mean.cor = x$table$pureCor %>% mean(na.rm=TRUE), 
          size=dim(x$table)[1],
          nonzero = mean(x$table$non0.1 + x$table$non0.2)/2)}
      }) %>% t %>% data.frame -> tmp

    lm(mean.cor ~ size, tmp %>% data.frame ) %>% summary
    
    lm(mean.cor ~ nonzero, tmp %>% data.frame ) %>% summary
    lm(nonzero ~ size, tmp %>% data.frame ) %>% summary
    
    tmp %>% ggplot(aes(size, nonzero)) + geom_point()
    tmp %>% ggplot(aes(size, mean.cor)) + geom_point()
    tmp %>% ggplot(aes(nonzero, mean.cor)) + geom_point()
    
### plots!!
    moreThanTen = TRUE #count only pairs with nonzero more than 10%.
    if (moreThanTen) {
      geneset.10.cor = rep(NA,30)
      pdf("output/cor.plots.10.pdf")
      } else {
      geneset.cor = rep(NA,30)
      pdf("output/cor.plots.pdf")
      }
    
    for (i in 1:30) { ## Fitting 1 ~ 30th genesets
      print(paste(i, "th geneset"))
      if (length(screen.zero(data=data, geneset = iset.Mm.c2[[i]], cut=.05, output="which"))<2) next
      a <- readRDS(paste0("output/geneset.",i,".list.rds"))
      
      # update names
      names(a$table)[4:15] <- c("a0","a1","a2", "b1", "b2", "p1", "p2", "p3", "p4", "lik","iter", "pureCor") 
      
      if (moreThanTen) {
        a$table <- a$table[a$table$non0.min>=.1,]
      }
      if (dim(a$table)[1] == 0) next
      
      #augmenting - for symmetry
      a2 <- a$table; a2[,c("1","2","non0.1","non0.2")] <- a2[,c("2","1","non0.2","non0.1")]
      
      a2 <- rbind(a$table, a2)
      cor.mean <- round(mean(a2$pureCor, na.rm=TRUE),3)
        if (moreThanTen) {
          geneset.10.cor[i] <- cor.mean
        } else {
          geneset.cor[i] <- cor.mean
        }
      
      cor.max = max(a2$pureCor)
      a2 %>% ggplot(aes(non0.1, non0.2, col=pureCor)) + geom_jitter() +
        scale_colour_gradient2() + ggtitle(paste0(i, "th geneset: ", length(a$genes), " genes (nonzero > 5%)")) -> p1
      a$table %>% ggplot(aes(pureCor)) + geom_histogram() + 
        geom_density(aes(fill="red"), alpha=.4) + xlim(c(0, max(cor.max, 0.4))) +
        ggtitle(paste0(i, "th geneset")) + guides(fill=FALSE) +
        annotate("text", x= 0.3, y = 10, label = paste("Mean cor = ", cor.mean), col="blue") -> p2
      
      
      
      #if (k %% 10 == 0) {print(paste("printing", k, "th out of", length(gene.id)))}
      
      grid.arrange( p1, p2, ncol=1)
      
    }
    dev.off()

### real vs random genesets
    # compare mean correlation of genesets and of random sets
    rbind(data.frame(cor = geneset.cor, set="geneset"),
          data.frame(cor = as.vector(as.matrix(randset.cor)), set=paste("random set",rep(1:10, each=30)))) %>% 
      ggplot(aes(set, cor, col=set)) + geom_boxplot() + geom_jitter(height=0)
    ggsave("../Documentation/plot4-1-1r.png")
    
    # sets w/ more than 10% of nonzero gene pairs
    rbind(data.frame(cor = geneset.10.cor, set="geneset (>10%)"),
          data.frame(cor = as.vector(as.matrix(randset.10.cor)), set=paste("random set",rep(1:10, each=30), "(>10%)"))) %>% 
      ggplot(aes(set, cor, col=set)) + geom_boxplot() + geom_jitter(height=0)
    ggsave("../Documentation/plot4-1-2r.png")
    
    
### relation btw setsize and avg.cor
    # real genesets
    data.frame(set.size=sapply(geneset.cor.list, function(x) ifelse(is.null(x),0,length(x$genes))),
               avg.cor = geneset.cor) %>% 
      ggplot(aes(set.size, avg.cor)) + geom_point() +
      ggtitle("set size : average correlation for genesets")
    ggsave("../Documentation/plot4-1-1-size.png")
    
    # random genesets
    data.frame(set.size=sapply(rand.set, length)[1:41],
               avg.cor = randset.cor) %>% 
      ggplot(aes(set.size, avg.cor)) + geom_point() +
    ggtitle("set size : average correlation for random sets")
    ggsave("../Documentation/plot4-1-2-size.png")
    
### relation btw nonzero proportion and avg.cor
    # real genesets
    data.frame(nonzero=sapply(geneset.cor.list, function(x) ifelse(is.null(x),NA,mean(x$table$non0.1 + x$table$non0.2)/2)),
               set.size=sapply(geneset.cor.list, function(x) ifelse(is.null(x),0,length(x$genes)))) %>% 
      ggplot(aes(nonzero, set.size)) + geom_point() +
      ggtitle("setsize : nonzero-proportion for genesets")
    ggsave("../Documentation/plot4-1-1-nonzero.png")
    
