#This script contains codes for bioinformatic analysis  for the pre-print 
#Quantitative landscapes reveal trajectories of cell-state transitions associated with drug resistance in melanoma
#Pillai, Chen et al., 2022


#################################################################

#Figure 2
#Pseudo-time analysis

#################################################################
library(monocle)
library(AUCell)
library(ggpubr)

##Read in required data
#Read-in a file containing genesets for Hyperpigmented, Invasive, NCSC and Transitory phenotypes (as defined in Rambow et al., Cell 2018)
#Each row has one geneset
geneset <- read.delim("clusters.txt")
geneSets <- list(genesets[1,],genesets[2,],genesets[3,],genesets[4,])
#Read in the expression data, rows contain genes and columns contain samples 
df <- read.csv("GSE115978.csv",row.names = 1)
#Read in the annotation data for GSE115978
annotations <- read.csv("GSE115978.csv",row.names = 1)
#Generating feature data
fd <- as.data.frame(rownames(df))
rownames(fd) <- fd$`rownames(df)`
names(fd) <- c("gene_short_name")
fd <- new("AnnotatedDataFrame", data = fd)
#Generating Phenotype data - kmeans clustering, we do not use this information for downstream analysis
clusters <- kmeans(df,5)
rownames(clusters) <- names(df)
names(clusters) <- c("Phenotype")
pd <- new("AnnotatedDataFrame", data = names(df))
##Pseudo-time analysis using monocle
cds <- newCellDataSet(as.matrix(df), phenoData = pd,featureData = fd, expressionFamily=negbinomial.size()) #Creating and normalizing data
cds <-  estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
#Selecting network genes for ordering
cds_ordering_genes <-  c("NFIC", "AHR", "JUN", "SMAD3", "KLF4", "TBX3", "NR3C1","MITF","FOS", "SMAD4","TFAP2A", "NR2F1", "MAFB", "ETV5", "STAT5A",
                         "FOXF1","TFE3") 
#Ordering of cells in dataset based on genes
cds <-setOrderingFilter(cds, ordering_genes = cds_ordering_genes)
#Dimensional reduction using DDRTree method
cds <- reduceDimension(cds, method = 'DDRTree')
cds <- orderCells(cds, reverse=TRUE) 
#Plotting of Fig. 2Di
plot_cell_trajectory(cds, color_by = "Pseudotime")+
  guides(fill=guide_legend(title="Pseudotime"))+
  theme(text=element_text(size=16))+
  scale_color_gradient2(low = "darkviolet", high="red", mid = "yellow", midpoint=(min(cds$Pseudotime) + (max(cds$Pseudotime)-min(cds$Pseudotime))/2))
ggsave("Pseudotime.png")
#Plotting of Fig. 2Di
plot_cell_trajectory(cds, color_by = annotations[names(df),1])+
  guides(fill=guide_legend(title="Pseudotime"))+
  theme(text=element_text(size=16))+rremove("legend")
ggsave("ITH.png")

#AUCell scoring of samples
cells_rankings <- AUCell_buildRankings(as.matrix(df))
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
AUC <- cells_AUC@assays@data$AUC
AUC <- round(AUC, digits = 2)
#Plotting AUCell score plots for each geneset
for (j in 1:nrow(AUC)) {
  plot_cell_trajectory(cds, color_by = AUC[j,])+
    guides(fill=guide_legend(title= "A"))+
    theme(text=element_text(size=16),legend.text=element_text(size=12))+
    scale_color_gradient2(low = "darkviolet", high="red", mid = "yellow", midpoint=(min(AUC[j,]) + (max(AUC[j,])-min(AUC[j,]))/2))
  ggsave(paste0(ds[i],"_",rownames(AUC)[j],"_pt.png"))
}


#################################################################

#Figure 3
#Pie charts and SDI calculation

#################################################################
library(vegan)
library(ggplot2)

#Plot each arm of dimensional reduction plot
plot_cell_trajectory(cds, color_by = "State") %>%
  scale_fill_discrete(labels = c("Invasive", "Proliferative", "Hyperdifferentiated"))
#Identify cells lying in each arm - i.e. each resistant path or proliferative arm
assignment <- as.data.frame(as.numeric(cds$State))
rownames(assignment) <-  rownames(pd@data)
#Combine meta-data on tumor of origin with phenotype of sample
data <- cbind(annotations[rownames(assignment),], assignment)
names(data) <- c("Tumor","Trajectory")
data <-table(data)
#Calculate Shannon's Diversity Index
div_index <- diversity(data)
#Calculate fraction of samples in each phenotype in each tumor
data <- as.data.frame(t(apply(data,1, function(x) x/sum(x))))
rownames(data) <- apply(data.frame(rownames(data), div_index), 1, paste, collapse=", SDI = ")
data <- data[complete.cases(data),]
data$Tumor <- rownames(data)
#Plot pie-charts
data <- reshape2::melt(data, id.vars='Tumor')
ggplot(data, aes(x="", y=value, fill=variable))+
  theme(axis.line=element_blank(), axis.text.x=element_blank(),axis.text.y=element_blank(),
        axis.ticks=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank())+
  geom_bar(stat = "identity")+
  coord_polar("y", start = 0)+
  facet_wrap(~Tumor, ncol = 7, labeller = label_both)+
  scale_fill_discrete(labels = c("Invasive", "Proliferative", "Hyperdifferentiated"))
ggsave("Piecharts.png")

