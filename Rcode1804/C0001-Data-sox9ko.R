###################################################################################################
### Section 1. Data step   by Dr. Di Wu
### Data: data <- sox9ko.pan.800 <- "sox9ko_paneth_800_cell_counts.txt"
### Gene sets: Mm.c2 (unofficial names), Mm.c2.symb (official names), 
###            iset.Mm.c2 (locations corresponding to data)
###################################################################################################

source("https://bioconductor.org/biocLite.R")
biocLite("biomaRt")
library(biomaRt)
  
# 1. data
  data <-read.delim("sox9ko_paneth_800_cell_counts.txt" , sep="\t")
  dim(data)

# 2. gene.emsl.mmu
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

# 3. Mm.c2, Mm.c2.symb: Genesets and symbols
  load("mouse_c2_v4.rdata")
  #Mm.c2

  Mm.c2.symb<-vector("list", length(Mm.c2))
  for (iset in 1:length(Mm.c2)){
    Mm.c2.symb[[iset]]<-unique(gene.emsl.mmu[gene.emsl.mmu[,2]%in%Mm.c2[[iset]], 1] )}
  names(Mm.c2.symb)<-names(Mm.c2)

  Mm.c2.symb[[1]]
  
# 4. iset.Mm.c2  
  iset.Mm.c2 = list()
  for (k in 1:length(Mm.c2))
  {iset.Mm.c2[[k]]<-which(  (data$Gene)  %in%   as.character(Mm.c2.symb[[k]][Mm.c2.symb[[k]]!=""])    )}
  names(iset.Mm.c2) <- names(Mm.c2)
  
  length(Mm.c2) # 4722 Gene sets

# IF fails, try saveRDS and readRDS.
  saveRDS(iset.Mm.c2, "iset.Mm.c2.rds")
