###################################################################################################
### Section 1. Data step   by Dr. Di Wu
### Data: data <- sox9ko.pan.800 <- "sox9ko_paneth_800_cell_counts.txt"
### Gene sets: Mm.c2 (unofficial names), Mm.c2.symb (official names), 
###            iset.Mm.c2 (locations corresponding to data)
###################################################################################################

source("https://bioconductor.org/biocLite.R")
biocLite("biomaRt")

library(biomaRt)

ensembl = useMart("ensembl")
listDatasets(ensembl) 
ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl")
filters = listFilters(ensembl)
attributes = listAttributes(ensembl) # entrezgene
attributes[1:5,]
gene.emsl.mmu=a=getBM(attributes=c("mgi_symbol","entrezgene"),    mart=ensembl   ) 
dim(gene.emsl.mmu)
table(is.na(gene.emsl.mmu[,2]))

#FALSE  TRUE 
#21366 27845

load("mouse_c2_v4.rdata")
ls()
#Mm.c2


Mm.c2.symb<-vector("list", length(Mm.c2))
for (iset in 1:length(Mm.c2)){
  Mm.c2.symb[[iset]]<-unique(gene.emsl.mmu[gene.emsl.mmu[,2]%in%Mm.c2[[iset]], 1] )}
names(Mm.c2.symb)<-names(Mm.c2)

Mm.c2.symb[[1]]

sox9ko.pan.800<-read.delim("sox9ko_paneth_800_cell_counts.txt" , sep="\t")
dim( sox9ko.pan.800)
data<-sox9ko.pan.800

iset.Mm.c2=NULL
for (k in 1:length(Mm.c2))
{iset.Mm.c2[[k]]<-which(  (data$Gene)  %in%   as.character(Mm.c2.symb[[k]][Mm.c2.symb[[k]]!=""])    )}
names(iset.Mm.c2) <- names(Mm.c2)

# "iset.Mm.c2" are index that tell which genes in the data matrix are in the ith gene set. 

###################################################################################################
### Section 2. pairwise association - MI & cor
###################################################################################################


###################################################################################################
### 2.1 MI & cor - computing associations
###################################################################################################
library(entropy)  # https://cran.r-project.org/web/packages/entropy/entropy.pdf
                  # https://en.wikipedia.org/wiki/Mutual_information
# function: pairwise.dep, mi.empirical.2

length(Mm.c2) # 4722 Gene sets

# Mm.c2.mi.tmp <- lapply(1:2, function(i) pairwise.dep(data[iset.Mm.c2[[i]],], mi.empirical.2, name = "mi"))
# Mm.c2.mi <- lapply(1:length(Mm.c2), function(i) pairwise.dep(data[iset.Mm.c2[[i]],], mi.empirical.2))
# Mm.c2.cor.tmp <- lapply(1:2, function(i) pairwise.dep(data[iset.Mm.c2[[i]],], cor, name = "cor"))
set.seed(100); a <- Sys.time() # 12:34 ~ 1.2 hours
Mm.c2.micor.tmp <- lapply(1:10, function(i) pairwise.dep(data[iset.Mm.c2[[i]],], list(cor, mi.empirical.2),name = c("cor","MI"), p.value=TRUE, n.sample=1000))
Sys.time() - a # 
# Mm.c2.micor.tmp.log <- lapply(1:2, function(i) pairwise.dep(data[iset.Mm.c2[[i]],], list(cor, mi.empirical.2),name = c("cor","MI"), log.data = TRUE))

  # cor & MI for the 1st & 6th pair: cor -0.01, MI 0.59
  Mm.c2.micor.tmp[[1]][5,]
  # nonzero cells of the 1st & 6th observations
  t(data[iset.Mm.c2[[1]][c(1,6)],data[iset.Mm.c2[[1]][1],]!=0 | data[iset.Mm.c2[[1]][6],]!=0])
  mi.empirical.2(data[iset.Mm.c2[[1]][1],-1], data[iset.Mm.c2[[1]][6],-1])
  freqs.empirical(rbind(data[iset.Mm.c2[[1]][1],-1], data[iset.Mm.c2[[1]][6],-1]))
  mi.plugin(freqs.empirical(rbind(data[iset.Mm.c2[[1]][1],-1], data[iset.Mm.c2[[1]][6],-1])))


###################################################################################################
### 2.2 MI & cor - computing associations - visualization (scatter plots, p-value plots)
###################################################################################################

library(ggplot2)
library(gridExtra)


# distribution of nonzero counts
hist(Mm.c2.micor.tmp[[1]]$non0.min)
non0.cut = c(0,0.01,0.05,0.1,1)
# Storing plots
plot2.2.1 <- lapply(1:10, function(x) ggplot(na.omit(transform(Mm.c2.micor.tmp[[x]], nonzero = cut(non0.min, non0.cut))), aes(MI, cor, label=pair)) +
                  facet_grid(. ~ nonzero, labeller=labeller(.rows = label_both, .cols = label_value)) +
                  geom_point(aes(color=1-non0.min)) + ggtitle(paste0("MI by zero counts - Geneset ",x)))
png("plot2-2-1.png", width = 720, height = 360, units = "px")
plot2.2.1[[1]]  #color = p-value(MI)
dev.off()

plot2.2.2 <- lapply(1:10, function(x) ggplot(na.omit(transform(Mm.c2.micor.tmp[[x]], nonzero = cut(non0.min, non0.cut))), aes(MI, cor, label=pair)) +
                  facet_grid(. ~ nonzero, labeller=labeller(.rows = label_both, .cols = label_value)) +
                  geom_point(aes(color=MI.p.val)) + ggtitle(paste0("MI by zero counts - Geneset ",x)))
png("plot2-2-2.png", width = 720, height = 360, units = "px")
plot2.2.2[[1]]  #color = p-value(MI) // p(cor) is meaningless (b/c cor and p(cor) are strictly monotone)
dev.off()

# For log scale, see the supplementary code

# Geneset 1: 26th & 50th pair: Cor=0.8, MI=~0. 18th obs has something.
data[iset.Mm.c2[[1]][c(26,50)],data[iset.Mm.c2[[1]][26],]!=0 | data[iset.Mm.c2[[1]][50],]!=0]

head(Mm.c2.micor.tmp[[1]])

# 2.2.3. p-value plots
png("plot2-2-3.png", width = 480, height = 240, units = "px")
plot2.2.3 <- plot(Mm.c2.micor.tmp[[1]]$cor.p.val,Mm.c2.micor.tmp[[1]]$MI.p.val,xlab="p(cor)", ylab="p(MI)")
dev.off()

# 2.2.3.1 explanation for the upper-right cluster
head(Mm.c2.micor.tmp[[1]][Mm.c2.micor.tmp[[1]]$cor.p.val>0.8 & Mm.c2.micor.tmp[[1]]$MI.p.val>=0.5,])
Mm.c2.micor.tmp[[1]][1:30,]
table(as.numeric(data[iset.Mm.c2[[1]][1],-1]),as.numeric(data[iset.Mm.c2[[1]][7],-1]))

# 2.2.3.2 bottom-right area (Insignificant cor but sig MI)
head(Mm.c2.micor.tmp[[1]][Mm.c2.micor.tmp[[1]]$cor.p.val>0.1 & Mm.c2.micor.tmp[[1]]$MI.p.val<0.05,])
# 1st & 27th pair: cor = -.004 (p=0.9), MI = 0.014 (p=0.04)
a <- table(as.numeric(data[iset.Mm.c2[[1]][1],-1]),as.numeric(data[iset.Mm.c2[[1]][27],-1])) %>% print
# what if we omit the (3,8) cell?
a[8,4] <- 0
mi.empirical(a)
#contingency table to vectors: as.data.frame(to get counts) & cout2cases (to get repetitions)
a <- countsToCases(as.data.frame(a))
set.seed(1)
p.permutation(a$Var1, a$Var2, measure=mi.empirical.2) # MI=0.004, (p=0.33)


# 1st & 42th pair: cor = -.02 (p=0.6), MI = 0.03 (p=0.004)
a <- table(as.numeric(data[iset.Mm.c2[[1]][1],-1]),as.numeric(data[iset.Mm.c2[[1]][42],-1])) %>% print
# what if we omit the (10,16) cell?  there are 3 such cells
a[12,9] <- 0
a <- countsToCases(as.data.frame(a))
set.seed(1)
p.permutation(a$Var1, a$Var2, measure=mi.empirical.2) # MI=0.020, (p=0.054)

# p-value is heavily influenced by a single observation

plot2.2.3.1 <- ggplot(mapping=aes(as.numeric(data[iset.Mm.c2[[1]][1],-1]),as.numeric(data[iset.Mm.c2[[1]][7],-1]))) +
  geom_jitter(height=0.05, width=0.1) + xlab("Gene1") + ylab("Gene7")
plot2.2.3.2 <- ggplot(mapping=aes(as.numeric(data[iset.Mm.c2[[1]][1],-1]),as.numeric(data[iset.Mm.c2[[1]][27],-1]))) +
  geom_jitter(height=0.05, width=0.1) + xlab("Gene1") + ylab("Gene27")
png("plot2.2.3.1.png", height=300, width=480)
grid.arrange(plot2.2.3.1, plot2.2.3.2, ncol=2)
dev.off()


###################################################################################################
### 2.3 MI & cor - gene-specific plots (all pairs with the genes)
###################################################################################################

# restoring the redundant pairs 1-2 => 1-2 & 2-1
tmp <- redundant(Mm.c2.micor.tmp[[1]])
unique(tmp[,c("1","non0.1")])  #proportion of nonzeros

## Sample genes within intervals of [0%,5%], [5%,10%], [10%,100%] nonzero counts.
plot2.3.1 <- ggplot(tmp[tmp$"1" %in% c(3,13,15),], aes(MI, cor, label=pair)) + 
  facet_grid(. ~ non0.1, labeller=labeller(.rows = label_both, .cols = label_value)) +
  geom_point(aes(color=non0.2)) + ggtitle(paste0("MI x cor by zero counts - Geneset1 - Gene 3, 13, 15 with non0 prob of second vector"))
plot2.3.1
png2(plot2.3.1)
  
# Dispersions of cor: same, Dispersion of MI varies (many nonzero implies higher MI)
# Given a gene, the proportion of nonzero of the other gene determines the distribution.

## Sample genes within intervals of [0%,1%],[1%,5%], [5%,10%], [10%,15%], [15%,20%], [25%,100%]  
plot2.3.2 <- ggplot(tmp[tmp$"1" %in% c(3,13,15,44,17,52),], aes(MI, cor, label=pair)) + 
  facet_grid(. ~ non0.1, labeller=labeller(.rows = label_both, .cols = label_value)) +
  geom_point(aes(color=non0.2)) + ggtitle(paste0("MI x cor by zero counts - Geneset1 - Gene 3, 13, 15, 44, .. with non0 prob of second vector"))
plot2.3.2


## same plot with p-values shown by color
# given same sample size, p-value-of-cor ==(equiv) cor
# enough to show p-value-of-MI
plot2.3.3 <- ggplot(tmp[tmp$"1" %in% c(3,13,15),], aes(MI, cor, label=pair)) + 
  facet_grid(. ~ non0.1, labeller=labeller(.rows = label_both, .cols = label_value)) +
  geom_point(aes(color=MI.p.val)) + ggtitle(paste0("MI x cor by zero counts - Geneset1 - Gene 3, 13, 15, 44, ... with p(MI)"))
plot2.3.3
plot2.3.4 <- ggplot(tmp[tmp$"1" %in% c(3,13,15),], aes(MI, cor, label=pair)) + 
  facet_grid(I(MI.p.val<0.05) ~ non0.1, labeller=labeller(.rows = label_both, .cols = label_value)) +
  geom_point(aes(color=MI.p.val)) + ggtitle(paste0("MI x cor by zero counts - Geneset1 - Gene 3, 13, 15, 44, ... with p(MI)"))
plot2.3.4  #Those with large p-values have mostly zero correlations.
png2(plot2.3.4)


###################################################################################################
### 2.3B MI & cor - gene-specific plots (given a random genesets)
###################################################################################################

# generating 10 random genesets (with 50 genes for each geneset).
set.seed(100)
rand.set <- lapply(1:10, function(x) sample(1:dim(data)[1], 50))
rand.set.micor.1 <- lapply(1:10, function(i) pairwise.dep(data[rand.set[[i]],], list(cor, mi.empirical.2),name = c("cor","MI"), p.value=FALSE, n.sample=1000))

# First random geneset. All pair of cor & MI
plot2.3B.1 <- lapply(1:10, function(x) ggplot(transform(rand.set.micor.1[[x]], nonzero = cut(non0.min, non0.cut)), aes(MI, cor, label=pair)) +
                      facet_grid(. ~ nonzero, labeller=labeller(.rows = label_both, .cols = label_value)) +
                      geom_point(aes(color=1-non0.min)) + ggtitle(paste0("MI by zero counts - random geneset ",x, " (facets = min zero proportion)")))
png2(plot2.3B.1[[1]])

non0.cut = c(0,0.001,0.01,0.05,0.1,1) # redefined to separate strictly 0.
# First random geneset. All pair of cor & MI
# tmp <- redundant(rand.set.micor.1[[1]])
# table(unique(tmp[,c("1","non0.1")])$non0.1)  #proportion of nonzeros
# out of 50, 13 has 100% zero, 15 has > 99% zero, 14 has >95% zero, 7 have >90% zero, 1 has <90% zero.
lev <- levels(cut(tmp$non0.1, non0.cut))
set.seed(1)  # sampling a specific gene from each group
spec.gene <- lapply(1:10, function(a) {
                tmp <- redundant(rand.set.micor.1[[a]])
                spec.gene <- sapply(2:5, function(x) sample(na.omit(tmp$"1"[cut(tmp$non0.1, non0.cut)==lev[x]]),1))
                return(spec.gene)
              })
spec.gene

## Sample genes within intervals of [0%,5%], [5%,10%], [10%,100%] nonzero counts.
## within-set comparison
plot2.3B.2 <- lapply(1:10, function(a) {
                tmp <- redundant(rand.set.micor.1[[a]])
                ggplot(tmp[tmp$"1" %in% spec.gene[[a]],], aes(MI, cor, label=pair)) + 
                facet_grid(. ~ non0.1, labeller=labeller(.rows = label_both, .cols = label_value)) +
                geom_point(aes(color=non0.2)) + ggtitle(paste0("MI x cor by zero counts - random geneset",a, "- Gene", spec.gene[[a]], "with non0 prob of second vector"))
})

png("plot2-3B-2.png", height=600, width=480)
grid.arrange(plot2.3B.2[[1]],plot2.3B.2[[2]],plot2.3B.2[[3]], nrow=3)
dev.off()

## between-set comparison
rand.set.micor.2 <- pairwise.dep(data[c(rand.set[[2]], rand.set[[1]][c(18,22,43,40)]),], list(cor, mi.empirical.2),name = c("cor","MI"), p.value=FALSE, n.sample=1000)
tmp <- redundant(rand.set.micor.2)
plot2.3B.3 <- ggplot(tmp[tmp$"1" %in% c(51,52,53,54) & tmp$"2"<=50,], aes(MI, cor, label=pair)) +
  facet_grid(. ~ non0.1, labeller=labeller(.rows = label_both, .cols = label_value)) + 
  geom_point(aes(color=non0.2)) + ggtitle(paste0("MI x cor - random set1 Genes (18,22,43,40) with random set2"))

rand.set.micor.3 <- pairwise.dep(data[c(rand.set[[3]], rand.set[[1]][c(18,22,43,40)]),], list(cor, mi.empirical.2),name = c("cor","MI"), p.value=FALSE, n.sample=1000)
tmp <- redundant(rand.set.micor.3)
plot2.3B.4 <- ggplot(tmp[tmp$"1" %in% c(51,52,53,54) & tmp$"2"<=50,], aes(MI, cor, label=pair)) +
  facet_grid(. ~ non0.1, labeller=labeller(.rows = label_both, .cols = label_value)) + 
  geom_point(aes(color=non0.2)) + ggtitle(paste0("MI x cor - random set1 Genes (18,22,43,40) with random set3"))

png("plot2-3B-3.png", height=400, width=480)
grid.arrange(plot2.3B.3,plot2.3B.4, nrow=2)
dev.off()





###################################################################################################
### 2.4 MI & cor - Aggregating counts (MI)
###################################################################################################

## aggregation (Sturge's rule : mi.empirical.3)
Mm.c2.micor.bin <- pairwise.dep(data[iset.Mm.c2[[1]],], list(cor, mi.empirical.3), name = c("cor","MI"), p.value=FALSE,n.sample=200)
head(Mm.c2.micor.bin)
# restoring the redundant pairs 1-2 => 1-2 & 2-1
tmp2 <- redundant(Mm.c2.micor.bin)
tmp2 <- rbind(tmp, tmp2)
tmp2$aggr <- rep(c("nominal","aggregated"), each=dim(tmp)[1])

## Sample genes within intervals of [0%,5%], [5%,10%], [10%,100%] nonzero counts.
## Original(nominal) + Aggregated
plot2.4.1 <- ggplot(tmp2[tmp2$"1" %in% c(3,13,15),], aes(MI, cor, label=pair)) + 
  facet_grid(aggr ~ non0.1, labeller=labeller(.rows = label_both, .cols = label_value)) +
  geom_point(aes(color=non0.2)) + ggtitle(paste0("MI x cor by zero counts - Geneset1 - Gene 3, 13, 15 with non0 prob of second vector"))
plot2.4.1
png2(plot2.4.1)

## Aggregation using log-transformation
## Equal #Bins strategy (Sturge's rule + log: mi.empirical.log)
Mm.c2.micor.log <- pairwise.dep(data[iset.Mm.c2[[1]],], list(cor, mi.empirical.log), name = c("cor","MI"), p.value=FALSE)  # no p-value b/c of long time
head(Mm.c2.micor.log)
## Equal Bin-size strategy (2^k - 2^{k-1} size)
Mm.c2.micor.log2 <- pairwise.dep(data[iset.Mm.c2[[1]],], list(cor, mi.empirical.log2), name = c("cor","MI"), p.value=FALSE)  # no p-value b/c of long time
head(Mm.c2.micor.log2)
# restoring the redundant pairs 1-2 => 1-2 & 2-1
tmp2 <- redundant(Mm.c2.micor.bin)
tmp3 <- redundant(Mm.c2.micor.log); tmp3[,c("cor.p.val","MI.p.val")] <- NA
tmp4 <- redundant(Mm.c2.micor.log2); tmp4[,c("cor.p.val","MI.p.val")] <- NA
tmp3 <- rbind(tmp, tmp2, tmp3, tmp4)
tmp3$aggr <- rep(c("nominal","aggr(original bin)","aggr(log, equal #Bin)","aggr(log, equal size)"), each=dim(tmp)[1])

## Original(nominal) + Aggregated(equal #Bin) + aggregated(equal Bin size)
plot2.4.2 <- ggplot(tmp3[tmp3$"1" %in% c(3,13,15,44,17,52),], aes(MI, cor, label=pair)) + 
  facet_grid(aggr ~ non0.1, labeller=labeller(.rows = label_both, .cols = label_value)) +
  geom_point(aes(color=non0.2)) + ggtitle(paste0("MI x cor by zero counts - Geneset1 - Gene 3, 13, 15 with non0 prob of second vector"))
plot2.4.2
png2(plot2.4.2, width = 480, height = 270)


## aggregation (Scott's rule : mi.empirical.3 with type=Scott)
Mm.c2.micor.bin2 <- pairwise.dep(data[iset.Mm.c2[[1]],], list(cor, mi.empirical.3.SC), name = c("cor","MI"), p.value=FALSE,n.sample=200, type="Scott")
head(Mm.c2.micor.bin2)
# restoring the redundant pairs 1-2 => 1-2 & 2-1
tmp2 <- redundant(Mm.c2.micor.bin2)
# tmp2 <- rbind(tmp, tmp2)
tmp2$aggr <- rep("aggregated", each=dim(tmp2)[1])

## Sample genes within intervals of [0%,5%], [5%,10%], [10%,100%] nonzero counts.
## Aggregated(Scott)
plot2.4.3 <- ggplot(tmp2[tmp2$"1" %in% c(3,13,15),], aes(MI, cor, label=pair)) + 
  facet_grid(. ~ non0.1, labeller=labeller(.rows = label_both, .cols = label_value)) +
  geom_point(aes(color=non0.2)) + ggtitle(paste0("MI x cor by zero counts - Geneset1 - Gene 3, 13, 15 with non0 prob of second vector")) +
  xlim(c(0,0.45))
plot2.4.3
png2(plot2.4.3, width = 480, height = 120)


###################################################################################################
### 2.5 NMI & cor - Normalized MI
###################################################################################################
# converting MI to NMI1
Mm.c2.micor.NMI <- lapply(1:10,function(x) {Mm.c2.micor.tmp[[x]]$MI <- NMI(Mm.c2.micor.tmp[[x]]$MI); return(Mm.c2.micor.tmp[[x]])})
head(Mm.c2.micor.tmp[[1]])
head(Mm.c2.micor.NMI[[1]])

# check if NMI > cor
(abs(Mm.c2.micor.NMI[[1]]$cor) <= Mm.c2.micor.NMI[[1]]$MI)[1:10]
Mm.c2.micor.NMI[[1]][9,]

plot2.5.1 <- ggplot(transform(Mm.c2.micor.NMI[[1]], nonzero = cut(non0.min, non0.cut)), aes(MI, abs(cor), label=pair)) +
                      facet_grid(. ~ nonzero, labeller=labeller(.rows = label_both, .cols = label_value)) +
                      geom_point(aes(color=MI.p.val)) + ggtitle(paste0("NMI x |cor| by zero counts - Geneset 1")) +
  geom_abline(intercept=0, slope=1, linetype, color="red", size=0.1)
png2(plot2.5.1, width = 480, height = 150)


# NMI2 (=MI/(HX+HY)), NMI3, ..., NMI5
Mm.c2.micor.NMI25 <- pairwise.dep(data[iset.Mm.c2[[1]],], list(cor, NMI2, NMI3, NMI4, NMI5),c("cor",paste0("NMI",2:5)),p.value=FALSE)
  Mm.c2.micor.NMI25[!is.finite(Mm.c2.micor.NMI25$NMI2),]  # check if there's some error
  Mm.c2.micor.NMI25[1:30,]                                # viewing the head

plot2.5.2 <- ggplot(transform(Mm.c2.micor.NMI25, nonzero = cut(non0.min, non0.cut)), aes(NMI2, abs(cor), label=pair)) +
  facet_grid(. ~ nonzero, labeller=labeller(.rows = label_both, .cols = label_value)) +
  geom_point() + #  ggtitle(paste0("NMI2 x |cor| by zero counts - Geneset 1")) +
  geom_abline(intercept=0, slope=1, linetype, color="red", size=0.1) + xlim(c(0,1.05))
plot2.5.2B <- ggplot(transform(Mm.c2.micor.NMI25, nonzero = cut(non0.min, non0.cut)), aes(NMI3, abs(cor), label=pair)) +
  facet_grid(. ~ nonzero, labeller=labeller(.rows = label_both, .cols = label_value)) +
  geom_point() + #  ggtitle(paste0("NMI3 x |cor| by zero counts - Geneset 1")) +
  geom_abline(intercept=0, slope=1, linetype, color="red", size=0.1) + xlim(c(0,1.05))
plot2.5.2C <- ggplot(transform(Mm.c2.micor.NMI25, nonzero = cut(non0.min, non0.cut)), aes(NMI4, abs(cor), label=pair)) +
  facet_grid(. ~ nonzero, labeller=labeller(.rows = label_both, .cols = label_value)) +
  geom_point() + #  ggtitle(paste0("NMI4 x |cor| by zero counts - Geneset 1")) +
  geom_abline(intercept=0, slope=1, linetype, color="red", size=0.1) + xlim(c(0,1.05))
plot2.5.2D <- ggplot(transform(Mm.c2.micor.NMI25, nonzero = cut(non0.min, non0.cut)), aes(NMI5, abs(cor), label=pair)) +
  facet_grid(. ~ nonzero, labeller=labeller(.rows = label_both, .cols = label_value)) +
  geom_point() + #  ggtitle(paste0("NMI5 x |cor| by zero counts - Geneset 1")) +
  geom_abline(intercept=0, slope=1, linetype, color="red", size=0.1) + xlim(c(0,1.05))
png("plot2-5-2.png", height=800, width=480)
grid.arrange(plot2.5.2B, plot2.5.2, plot2.5.2D, plot2.5.2C, nrow=4)
dev.off()

  # Those with MNI4 almost equal to 1
  Mm.c2.micor.NMI25[Mm.c2.micor.NMI25$NMI4 > 0.9,]  #7-19 (cor, NMI2, 3, 4, 5) = (0.28, 0.012, 0.012, 1, 0.11)
    entropy(table(extractor(7))); entropy(table(extractor(19))); entropy(table(extractor(7),extractor(19))); mi.empirical.2(extractor(7),extractor(19))
    # HX = 0.0096, HY = 0.781, HXY = 0.781, MI = 0.0096
    table(extractor(7),extractor(19))
    plot2.5.2.1 <- ggplot(mapping=aes(extractor(7),extractor(19))) +
    geom_jitter(height=0.05, width=0.1) + xlab("Gene7") + ylab("Gene19")
    png2(plot2.5.2.1)
    table(extractor(9),extractor(57))
    table(extractor(41),extractor(17))

  # Those with MNI4 above 0.6 while min.non0 in (0.01,0.05)
  Mm.c2.micor.NMI25[Mm.c2.micor.NMI25$NMI4 > 0.6 & Mm.c2.micor.NMI25$non0.min >0.01,]  #15-38 (cor, NMI2, 3, 4, 5) = (0.12, 0.042, 0.044, 0.67, 0.17)
  entropy(table(extractor(15))); entropy(table(extractor(38))); entropy(table(extractor(15),extractor(38))); mi.empirical.2(extractor(15),extractor(38))
  # HX = 0.114, HY = 1.71, HXY = 1.75, MI = 0.077
  table(extractor(15),extractor(38))
  plot2.5.2.1 <- ggplot(mapping=aes(extractor(7),extractor(19))) +
    geom_jitter(height=0.05, width=0.1) + xlab("Gene7") + ylab("Gene19")
  table(extractor(9),extractor(57))
  table(extractor(41),extractor(17))



###################################################################################################
### 2.6 Quantifying the joint distn of MI & cor
###################################################################################################
head(Mm.c2.micor.tmp[[1]])

sapply(Mm.c2.micor.tmp,polar2)   #polar2: getting weighted mean of theta after converting Cartesian to polar

# adding theta and norm (polar3)
Mm.c2.micor.tmp <- lapply(Mm.c2.micor.tmp, polar3)

plot2.6.1 <- lapply(1:10, function(x) ggplot(na.omit(transform(Mm.c2.micor.tmp[[x]], nonzero = cut(non0.min, non0.cut))), aes(theta, norm, label=pair)) +
  facet_grid(nonzero ~ ., labeller=labeller(.rows = label_both, .cols = label_value)) +
  geom_point(aes(color=MI.p.val)) + ggtitle(paste0("Geneset ", x)))
plot2.6.1[[2]]
png("plot2-6-1.png", height=720, width=720)
grid.arrange(plot2.6.1[[1]], plot2.6.1[[2]], plot2.6.1[[3]], ncol=3)
dev.off()

tmp <- transform(Mm.c2.micor.tmp[[1]], nonzero = cut(non0.min, non0.cut))
sapply(lev, function(x) polar2(tmp[tmp$nonzero == x,]))

tmp <- transform(Mm.c2.micor.tmp[[2]], nonzero = cut(non0.min, non0.cut))
sapply(lev, function(x) polar2(tmp[tmp$nonzero == x,]))

tmp <- transform(Mm.c2.micor.tmp[[3]], nonzero = cut(non0.min, non0.cut))
sapply(lev, function(x) polar2(tmp[tmp$nonzero == x,]))

# NMI2 + polar3
Mm.c2.micor.NMI25 <- Mm.c2.micor.NMI25[,-(9:13)]
Mm.c2.micor.NMI25 <- polar3(Mm.c2.micor.NMI25, a = "cor", b = "NMI2")
plot2.6.2 <- ggplot(na.omit(transform(Mm.c2.micor.NMI25, nonzero = cut(non0.min, non0.cut))), aes(theta, norm, label=pair)) +
                      facet_grid(nonzero ~ ., labeller=labeller(.rows = label_both, .cols = label_value)) +
                      geom_point() + ggtitle("Geneset 1 - NMI2")
Mm.c2.micor.NMI25 <- polar3(Mm.c2.micor.NMI25, a = "cor", b = "NMI3")
plot2.6.3 <- ggplot(na.omit(transform(Mm.c2.micor.NMI25, nonzero = cut(non0.min, non0.cut))), aes(theta, norm, label=pair)) +
  facet_grid(nonzero ~ ., labeller=labeller(.rows = label_both, .cols = label_value)) +
  geom_point() + ggtitle("Geneset 1 - NMI3")
Mm.c2.micor.NMI25 <- polar3(Mm.c2.micor.NMI25, a = "cor", b = "NMI4")
plot2.6.4 <- ggplot(na.omit(transform(Mm.c2.micor.NMI25, nonzero = cut(non0.min, non0.cut))), aes(theta, norm, label=pair)) +
  facet_grid(nonzero ~ ., labeller=labeller(.rows = label_both, .cols = label_value)) +
  geom_point() + ggtitle("Geneset 1 - NMI4")
Mm.c2.micor.NMI25 <- polar3(Mm.c2.micor.NMI25, a = "cor", b = "NMI5")
plot2.6.5 <- ggplot(na.omit(transform(Mm.c2.micor.NMI25, nonzero = cut(non0.min, non0.cut))), aes(theta, norm, label=pair)) +
  facet_grid(nonzero ~ ., labeller=labeller(.rows = label_both, .cols = label_value)) +
  geom_point() + ggtitle("Geneset 1 - NMI5")
png("plot2-6-2.png", height=720, width=720)
grid.arrange( plot2.6.2, plot2.6.4, plot2.6.5, ncol=3)
dev.off()

sapply(list("NMI2","NMI4","NMI5"), function (b) polar2(Mm.c2.micor.NMI25, b=b))
tmp <- transform(Mm.c2.micor.NMI25, nonzero = cut(non0.min, non0.cut))
sapply(lev, function(x) polar2(tmp[tmp$nonzero == x,], b="NMI2"))
sapply(lev, function(x) polar2(tmp[tmp$nonzero == x,], b="NMI4"))
sapply(lev, function(x) polar2(tmp[tmp$nonzero == x,], b="NMI5"))


###################################################################################################
### 2.7 Parametric approach of measuring MI
###################################################################################################
### See code "functions-BivariatePoisson.R"

##############################################################################################
### 2.7.1 Parametric approach of measuring MI - BvPoi
##############################################################################################
MLE.Geneset1 <- list()
a <- Sys.time()
MLE.Geneset1$BP <- pairwise.MLE (data[iset.Mm.c2[[1]],], ML.fun = ML.BP)  ## 34 min
MLE.Geneset1$BZIP <- pairwise.MLE (data[iset.Mm.c2[[1]],], ML.fun = ML.BZIP) ## ?? min
Sys.time() - a

Mm.c2.micor.BP <- pairwise.dep(data[iset.Mm.c2[[1]],], list(cor.ML.BP, MI.ML.BP), name = c("cor.BP","MI.BP"), p.value=FALSE)
Sys.time() %>% print - a # 2 hours
Mm.c2.micor.BP$cor <- Mm.c2.micor.tmp[[1]]$cor   # adding Pearson correlation (PC)

head(Mm.c2.micor.BP)
plot2.7.1 <-  ggplot(transform(Mm.c2.micor.BP, nonzero = cut(non0.min, non0.cut)), aes(MI.BP, cor.BP, label=pair)) +
                      facet_grid(. ~ nonzero, labeller=labeller(.rows = label_both, .cols = label_value)) +
                      geom_point(aes(color=1-non0.min)) + ggtitle(paste0("parametric MI & parametric cor by zero counts - Geneset 1"))
png2(plot2.7.1, width = 720, height = 360)
plot2.7.1.2 <-  ggplot(transform(Mm.c2.micor.BP, nonzero = cut(non0.min, non0.cut)), aes(MI.BP, cor, label=pair)) +
  facet_grid(. ~ nonzero, labeller=labeller(.rows = label_both, .cols = label_value)) +
  geom_point(aes(color=1-non0.min)) + ggtitle(paste0("parametric MI & PC by zero counts - Geneset 1"))
png2(plot2.7.1.2, width = 720, height = 360)

png("plot2-7-2.png", width = 720, height = 720, units = "px")
par (mfrow=c(2,1))
plot( Mm.c2.micor.BP$MI.BP, Mm.c2.micor.tmp[[1]]$MI, xlab = "MI(BP)", ylab="Empirical MI")
abline(coef=c(0,1))
plot( Mm.c2.micor.BP$cor.BP, Mm.c2.micor.tmp[[1]]$cor, xlab = "cor(BP)", ylab="Pearson correlation")
abline(coef=c(0,1))
dev.off()

##############################################################################################
### 2.7.2 Parametric approach of measuring MI - BvZIP
##############################################################################################

