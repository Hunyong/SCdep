### Supplementary codes
### not main codes but storing these because it is wasteful to throw the relevant codes which might be needed later.


### 1.01 pairwise association macro - visualization 1 corrplot
### 1.02 pairwise association macro - visualization 2 MI by Corr plot
### 1.03 pairwise association macro - visualization 3 tSNE
###################################################################################################




###################################################################################################
### 1.01 pairwise association macro - visualization 1 correlation heatmap
###################################################################################################
library(corrplot)
# function tri2mat
tmp.mi.1 <- tri2mat(Mm.c2.mi.tmp[[1]][,2])
tmp.mi.2 <- tri2mat(Mm.c2.mi.tmp[[2]][,2])
tmp.cor.1 <- tri2mat(Mm.c2.cor.tmp[[1]][,2])

par(mfrow=c(1,2))
corrplot(tmp.mi.1, method="color")
corrplot(tmp.mi.2, method="color")
corrplot(tmp.cor.1, method="color")



###################################################################################################
### 1.02 pairwise association macro - visualization 2 MI by Corr plot
###################################################################################################

# example
plot(x = Mm.c2.mi.tmp[[1]][,2], y = Mm.c2.cor.tmp[[1]][,2])
ggplot(Mm.c2.micor.tmp[[1]], aes(MI,cor,label=pair)) +
  scale_color_continuous(name="proportion of zero") +
  geom_text(aes(color=1-non0.min)) #+ geom_jitter()

ggplot(Mm.c2.micor.tmp[[1]], aes(1,cor,label=pair)) +
  facet_grid(. ~ I(non0.min>.1) + I(non0.min>.05), labeller=labeller(.rows = label_both, .cols = label_value)) +
  geom_jitter(aes(color=1-non0.min)) + ggtitle("Correlation by zero counts")

ggplot(Mm.c2.micor.tmp[[1]], aes(1,MI,label=pair)) +
  facet_grid(. ~ I(non0.min>.1) + I(non0.min>.05), labeller=labeller(.rows = label_both, .cols = label_value)) +
  geom_jitter(aes(color=1-non0.min)) + ggtitle("MI by zero counts")

# original scale
ggplot(Mm.c2.micor.tmp[[1]], aes(MI, cor, label=pair)) +
  facet_grid(. ~ I(non0.min>.1) + I(non0.min>.05), labeller=labeller(.rows = label_both, .cols = label_value)) +
  geom_point(aes(color=1-non0.min)) + ggtitle("MI x Correlation by zero counts")

# log scale
ggplot(Mm.c2.micor.tmp.log[[1]], aes(MI,cor,label=pair)) +
  facet_grid(. ~ I(non0.min>.1) + I(non0.min>.05), labeller=labeller(.rows = label_both, .cols = label_value)) +
  geom_point(aes(color=1-non0.min)) + ggtitle("MI x Correlation by zero counts (log scale)")
# geom_text(aes(color=1-non0.min))

# MI by cor by density(zero count)
library("scatterplot3d")
scatterplot3d(Mm.c2.micor.tmp[[1]][,c(6,2,3)])



###################################################################################################
### 1.03 pairwise association macro - visualization 3 tSNE
###################################################################################################
library(tsne)   # https://cran.r-project.org/web/packages/tsne/tsne.pdf
library(ggplot2)

# example (NA filled with 0)
tmp.mi.1[is.na(tmp.mi.1)] <- 0
tmp.cor.1[is.na(tmp.cor.1)] <- 0
tmp.mi.1.tsne <-  tsne(tmp.mi.1); tmp.mi.1.tsne <- as.data.frame(tmp.mi.1.tsne)
tmp.cor.1.tsne <-  tsne(tmp.cor.1); tmp.cor.1.tsne <- as.data.frame(tmp.cor.1.tsne)
plot(tmp.mi.1.tsne); text(tmp.mi.1.tsne)
plot(tmp.cor.1.tsne); text(tmp.cor.1.tsne)
ggplot(as.data.frame(tmp.mi.1.tsne), aes(V1,V2,label=rownames(tmp.mi.1.tsne))) +
  geom_text() #+ geom_jitter()
ggplot(as.data.frame(tmp.cor.1.tsne), aes(V1,V2,label=rownames(tmp.cor.1.tsne))) +
  geom_text() #+ geom_jitter()