### This is solely for C0101-BZINB-randomsets(supp).R file

### input the simulated result file "sim.result"[[firstlevel]][[secondlevel]]
  ## first level: 1 correlation (nonzero cutoff 5%), 2 correlation (nonzero cutoff 10%), 
  ##              3 nonzero info (nonzero cutoff 5%), 4 nonzero info (nonzero cutoff 10%)
  ## second level: simulation replicates 1:60

### parameters: mean, quantile(q), cutoff = .05, .1
  
  backgroundNoisePlot <- function(sim.result, cutoff = 0.05, FUN = "mean", p = 0.5, cor.lim=c(-0.01, 0.4), n.sim=10, index=NULL) {
    
    if (is.null(index)) {
      index = c(1,4) 
      if (cutoff == 0.1) {index = index + 1} else if (cutoff == 0.01) {index = index + 2} else if
      (cutoff == 0.02) {index = c(1,2)}
    }
    
    setsize = if (cutoff == 0.05) {
        geneset.size.05
      } else if (cutoff==.1) {
        geneset.size.10
      } else if (cutoff == .01) {
        geneset.size.01
      } else if (cutoff == .02) {
        geneset.size.02}
    if (FUN == "mean") {
      FUN = function(x) mean(x, na.rm=TRUE)
      fun.name = "mean"
    } else if (FUN == "quantile") {
      FUN = function(x) {
        x = x[!is.na(x)]
        quantile(x, probs = p)
      }
      fun.name = paste0(p*100, "th quantile")
    } else if (FUN == "logMean") {
      FUN = function(x) {
        x = x[!is.na(x)]
        2^mean(log2(x))
      }
      fun.name = "geographic mean"
    }
    
    ### getting %non-zero
    sapply(1:60, function(i) 
    {if (is.null(sim.result[[index[2]]][[i]])) rep(NA, n.sim + 1) else
      sapply(sim.result[[index[2]]][[i]], function(sim) {
        if (is.null(sim)) NA else {mean(as.matrix(sim))}
      })}) %>% matrix(n.sim + 1, 60, dimnames=list(c("real", rep("random",n.sim)), 1:60)) -> nonzero
  #print(nonzero)
    ### mean + set.size + %nonzero
    tmp = sapply(sim.result[[index[1]]], function(sim) {if (is.null(sim)) {rep(NA,n.sim + 1)} else {apply(sim, 2, FUN)}}) %>% as.vector
    tmp<- data.frame(cor = tmp, group = c("real", rep("random",n.sim)), set = rep(1:60, each=n.sim + 1))
    tmp$setSize <- rep(setsize[1:60], each=n.sim + 1)
    #tmp$set.and.size <- paste(tmp$setSize,tmp$set, sep = ".")

    # print(tmp %>% head)    
    tmp$nonzero <- as.vector(nonzero)
    # global avg by group
    cor.rand <- tmp%>%filter(group=="random")%>%"["("cor") %>% (function(x) {x[!is.na(x)]}) %>% mean
    cor.real <- tmp%>%filter(group=="real")%>%"["("cor") %>% (function(x) {x[!is.na(x)]}) %>% mean
    
    tmp %>% filter(group=="real") %>% "["("cor") %>% order -> set.order
    tmp$set2 <- factor(tmp$set, levels = set.order)  #set2 is the releveled factor
    
    # first 50 non-NA sets
    set.50 <<- tmp %>% filter(group=="real", !is.na(cor)) %>% "["(1:50, "set")
# tmp <<- tmp
    # arrange(group) is needed to have random set points plotted first than the real set points.
    tmp %>% arrange(group) %>% filter(set %in% set.50) %>% 
      ggplot(aes(factor(setSize), cor, col=factor(group))) + 
      geom_hline(yintercept = cor.rand, col="red") + 
      geom_hline(yintercept = cor.real, col="#33CCCC") +
      geom_boxplot(aes(x=factor(setSize), y=cor), outlier.shape = NA, data=subset(tmp, group=="random")) +
      geom_point(aes(size = factor(group))) + 
      scale_size_manual(values=c(0.3, 1.5)) + 
      ylim(cor.lim) + theme(legend.position = "bottom") +
      xlab ("Set size (# genes)") + 
      ylab (paste0("Correlation (", fun.name, ")")) +
      ggtitle(paste0("Real set V.S. background correlation \n (Summarized by ", fun.name, ", nonzero cutoff = ", cutoff, ")")) 
  }
  if (FALSE) {
    backgroundNoisePlot(sim.result=sim.result, cutoff=0.05, FUN = "quantile", p =.5)
  }
  
  ### p-values given a function (e.g. mean)
  p.functional <- function(object, FUN = mean, quantile = 0.5, real.position = 1) {
    n.sim = dim(object)[2]
    result = apply(object, 2, FUN, probs = quantile, na.rm=TRUE) %>% percent_rank %>% "[["(real.position)
    result = (1 - result*n.sim/(n.sim+1)) #flip 0-1 to 1-0, and adjust for continuity
    result
  }
  if (FALSE) {
    p.functional(sim.result[[1]][[3]])
  }
  sapply(sim.result[[1]], p.functional)
  