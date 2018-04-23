

### a random geneset
## 1. sampling & screening
  set.seed(1)
  rand.set <- sample(1:dim(data)[1], 30)
  
  gene.tmp <- screen.zero(data=data, geneset = rand.set, cut=.05, output="which")
  if (length(gene.tmp) <2) next
  data.tmp <- data[gene.tmp,]
  
## 2. fitting (MLE)  
  geneset <- pairwise.MLE(data = data.tmp, ML.fun = ML.BvZINB4.2, nondata.col = 1, showFlag=FALSE, tol=1e-5, cor.conv=TRUE, maxiter = 1000)
  geneset.list <- list(cutoff = "nonzero > 5%",
                       genes = gene.tmp, 
                       table = geneset)

## 3. p-value distribution.
  pdf("output/cor.plots-random.pdf")
  {
    print(paste("random geneset"))
    a <- geneset.list
    i = "random"
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
  