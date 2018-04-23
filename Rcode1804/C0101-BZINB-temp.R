library(ggplot2); library(gridExtra)
for (i in 1:30) { ## Fitting 1 ~ 30th genesets
  print(paste(i, "th geneset"))
  gene.tmp <- screen.zero(data=data, geneset = iset.Mm.c2[[i]], cut=.05, output="which")
  if (length(gene.tmp) <2) next
  data.tmp <- data[gene.tmp,]
  # read previous output
### change number here!!!
  a <- readRDS("output/geneset.1.list.rds")
  
  # update names
  names(a$table)[4:15] <- c("a0","a1","a2", "b1", "b2", "p1", "p2", "p3", "p4", "lik","iter", "pureCor") 
  
  # recalculate the pureCor, and add se, pval.
  pureCor = apply(a$table[1:2,], 1, function(x) {
    pureCor (as.numeric(x[4:12]), data.tmp[as.numeric(x[1]), -1], data.tmp[as.numeric(x[2]), -1])})
  a$pureCor = t(pureCor)[,1]
  a$pureCor.se = t(pureCor)[,2]
  a$pureCor.p = t(pureCor)[,3]
  
  # add deviance
  dev = apply(a$table, 1, function(x) {
print(x)
    dev.BvZINB4(  as.numeric(data.tmp[as.numeric(x[1]), -1]), 
                  as.numeric(data.tmp[as.numeric(x[2]), -1]),
                         model.lik = as.numeric(x[13]))})
  ### DEBUG ERRORS HERE !!!!!!          
  lik = dev
}

pdf("output/cor.plots.pdf")
for (i in 1:30) { ## Fitting 1 ~ 30th genesets
  print(paste(i, "th geneset"))
  a <- readRDS(paste0("output/geneset.",i,".list.rds"))
  # update names
  names(a$table)[4:15] <- c("a0","a1","a2", "b1", "b2", "p1", "p2", "p3", "p4", "lik","iter", "pureCor") 
  
  #augmenting - for symmetry
  a2 <- a$table; a2[,c("1","2","non0.1","non0.2")] <- a2[,c("2","1","non0.2","non0.1")]
  a2 <- rbind(a$table, a2)

  a2 %>% ggplot(aes(non0.1, non0.2, col=pureCor)) + geom_jitter() +
    scale_colour_gradient2() + ggtitle(paste0(i, "th geneset: ", length(a$genes), " genes (nonzero > 5%)")) -> p1
  a$table %>% ggplot(aes(pureCor)) + geom_histogram() + 
    geom_density(aes(fill="red"), alpha=.4) +
    ggtitle(paste0(i, "th geneset")) + guides(fill=FALSE) -> p2
  
  
  
  #if (k %% 10 == 0) {print(paste("printing", k, "th out of", length(gene.id)))}
  
  grid.arrange( p1, p2, ncol=1)
  
} 

dev.off()
