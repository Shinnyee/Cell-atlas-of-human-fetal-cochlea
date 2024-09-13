rm(list=ls())
library(Seurat)
library(rtracklayer)
library(tibble)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(Hmisc)
library(cowplot)
library(zoo)
library(ggplot2)
library(loomR)
library(scater)
library(stringr)
library(corrplot)
library(matrixStats)
library(eulerr)
library(viridis)
library(pvclust)
library(parallel)
library(dendextend)
library(future)
plan("multicore", workers = 4)
options(future.globals.maxSize= 62914560000) # 60GB (60000*1024^2)
####################################################################################################
###################################################################################################
#####################################################################################

library(SingleCellExperiment)
library(Matrix)
HUMAN_cochlea_combined_filtered <- readRDS("C:/Users/Dell/Desktop/WORKPLACE/R/HUMAN_cochlea_combined_filtered.rds")
#for mouse, we used datasets of both P25 OUR RNA-SEQ AND P28 JZ PUBLIC RNA-SEQ
seurat_mouse <- readRDS("C:/Users/Dell/Desktop/WORKPLACE/R/seurat_mouse_coch.rds")
table(seurat_mouse$subclass2)
#rename subclasses
Idents(seurat_mouse) <- "subclass2"
table(Idents(seurat_mouse))
new.cluster.ids <- c("TBC","BC","DC_OPC","MC","PC","Fibro3","Fibro1/2",
                     "HC","Fibro4","IPh_IBC","SGN","IC","OB/M",
                     "He","Fib","GC","SMC")
names(new.cluster.ids) <- levels(seurat_mouse)
seurat_mouse<- RenameIdents(seurat_mouse, new.cluster.ids)
DimPlot(seurat_mouse, reduction = "umap",label = TRUE,repel = TRUE)
seurat_mouse$subclasses_label <- Idents(seurat_mouse)
p1=DimPlot(seurat_mouse, reduction = "umap",label = TRUE,repel = TRUE,group.by = "subclasses_label")
p2=DimPlot(seurat_mouse, reduction = "umap",label = TRUE,repel = TRUE,group.by = "subclass2")
p1+p2
#rename subtype
DefaultAssay(seurat_mouse) <- "integrated"
seurat_mouse <- RunUMAP(seurat_mouse, dims = 1:30)
table(Idents(seurat_mouse))
seurat_mouse <- FindNeighbors(seurat_mouse,dims = 1:50)
seurat_mouse <- FindClusters(seurat_mouse,resolution = 1)
DimPlot(seurat_mouse, label = TRUE)
p1=DimPlot(seurat_mouse, reduction = "umap",label = TRUE,repel = TRUE,group.by = "subclasses_label")
p2=DimPlot(seurat_mouse, reduction = "umap",label = TRUE,repel = TRUE,group.by = "seurat_clusters")
p1+p2
table(seurat_mouse$seurat_clusters,seurat_mouse$subclasses_label)

a <- table(seurat_mouse$subclasses_label,seurat_mouse$seurat_clusters)
a <- as.data.frame.array(a)
write.csv(a,file = "cell_proportion_of_mouse_cluster_label.csv")

Idents(seurat_mouse) <- "seurat_clusters"
table(Idents(seurat_mouse))
new.cluster.ids <- c("BC_1","TBC_1","DC_OPC_1","PC_1","TBC_2","Fibro3","MC_1",
                     "Fibro1/2","OHC","PC_2","Fibro4","IPh_IBC_1","IC",
                     "SGN_1","OB/M","BC_2","SGN_2","DC_OPC_2","MC_2","DC_OPC_3",
                     "HeC_1","DC_OPC_4","IPh_IBC_2","IHC","PC_3","Fib_1","Fib_2",
                     "GC","HeC_2","IPh_IBC_3","SGN_3","SMC","SGN_4")
names(new.cluster.ids) <- levels(seurat_mouse)
seurat_mouse<- RenameIdents(seurat_mouse, new.cluster.ids)
DimPlot(seurat_mouse, reduction = "umap",label = TRUE,repel = TRUE)
seurat_mouse$cluster_label <- Idents(seurat_mouse)
p1=DimPlot(seurat_mouse, reduction = "umap",label = TRUE,repel = TRUE,group.by = "cluster_label")
p2=DimPlot(seurat_mouse, reduction = "umap",label = TRUE,repel = TRUE,group.by = "subclasses_label")
p1+p2

DefaultAssay(seurat_mouse) <- "SCT"
Idents(seurat_mouse) <- factor(Idents(seurat_mouse), 
                               levels = c("OHC","IHC","DC_OPC_1","DC_OPC_2","DC_OPC_3","DC_OPC_4",
                                          "IPh_IBC_1","IPh_IBC_2", "IPh_IBC_3",
                                          "TBC_1","TBC_2",
                                          "PC_1","PC_2","PC_3",
                                          "HeC_1","HeC_2",
                                          "SGN_1","SGN_2","SGN_3","SGN_4",
                                          "GC",
                                          "IC",
                                          "MC_1","MC_2",
                                          "BC_1","BC_2",
                                          "Fib_1","Fib_2",
                                          "Fibro1/2",
                                          "Fibro3", 
                                          "Fibro4",
                                          "SMC",
                                          "OB/M"
                                          
                                          
                                          
                                          
                               ))
markers.to.plot<-c("TMC1","SLC26A5","OTOF","MYO7A","SLC17A8","POU4F3",
                   "TSC22D1","EPYC","FBXO2","TMOD1","FGFR3","GATA3","GJB2",'DPYSL2',
                   "DGKB","CEACAM16","SLC1A3","S100A6","OTOG","OTOGL",
                   "EMILIN2",	"AXIN2",	"TJP1",	"ATP1A2","OTOR",
                   "ENAH","CLU","NUDT4","PTGDS",
                   "ANXA1","KCNQ5",
                   "NEFH","NEFL","SNAP25","PVALB","TUBB3","CALB2","NTNG1",
                   "MBP",	"PMP22",	"PLP1",
                   "MPZ",
                   "DCT","TYR","MET",	"KCNJ13",	"KCNJ10",
                   "DCLK1","ESRRB","KCNQ1",	"KCNE1","DCN",	"ITIH5",
                   "EMCN","ARHGAP6","CLDN11",
                   
                   "CREB5","ITGA8",	"IFIT3",
                   "VEPH1","CHST9",	"POSTN",	"LUM","PDGFRB","MAP1B","NEBL",
                   "COCH","SLC4A10","SLC4A11",	"COL9A1",	"COL9A2",
                   "SLC7A11","SLC13A3","KCNK2",
                   "ATP2B2","SLC26A7","MEIS1","DKK2",
                   "CEMIP","THSD4",
                   "RUNX2","COL1A2","RGS5","COL5A2",
                   "CD163","MEF2C","CX3CR1",	"C1QA"
                   
)
DotPlot(seurat_mouse, features = markers.to.plot, dot.scale = 8,
        cols  =c("white", "#ad9300")) + RotatedAxis()

markers.to.plot<-c("TMC1","MYO7A","GJB2","GATA3","TSC22D1","EPYC","OTOG","OTOGL","EMILIN2",
                   "RARRES1","NEFL","SNAP25","MBP","PLP1",
                   "MPZ","PMP22","DCT","TYR","DCLK1","ESRRB","TJP1",
                   "CLDN11", "ITGA8",	"CREB5",   "VWF","MECOM",
                   "VEPH1","CHST9","COCH","SLC4A10","SLC7A11","SLC13A3",
                   "ATP2B2","SLC26A7",
                   "CEMIP","THSD4","RUNX2","COL1A2","CD163","MEF2C",
                   "IGFBP7","COL4A1",
                   "HBG2","HBM","TOP2A","HMGB2")
DotPlot(seurat_mouse, features = markers.to.plot, dot.scale = 8,
        cols  =c("white", "#ad9300")) + RotatedAxis()
dev.off()
saveRDS(seurat_mouse,file = "seurat_mouse_new.rds")
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
seurat_human=HUMAN_cochlea_combined_filtered
rm(HUMAN_cochlea_combined_filtered)
table(Idents(seurat_human))
#rename subclasses
seurat_human$subclasses_label <- seurat_human$subclasses
p1=DimPlot(seurat_human,group.by = "subclasses_label",label = TRUE)
p2=DimPlot(seurat_human,group.by = "subclasses",label = TRUE)
p3=DimPlot(seurat_human,group.by = "seurat_clusters",label = TRUE)
p1+p2
p1+p3
table(seurat_human$subclasses_label, seurat_human$seurat_clusters)


dev.off()
# re-clustering
DefaultAssay(seurat_human) <- "integrated"
seurat_human <- RunUMAP(seurat_human, dims = 1:100)
table(Idents(seurat_human))
seurat_human <- FindNeighbors(seurat_human,dims = 1:100)
seurat_human <- FindClusters(seurat_human,resolution = 2.6)
#rename subclasses
p1=DimPlot(seurat_human,group.by = "subclasses_label",label = TRUE)
p2=DimPlot(seurat_human,group.by = "subclasses",label = TRUE)
p3=DimPlot(seurat_human,group.by = "seurat_clusters",label = TRUE)
p1+p2
p1+p3

dev.off()
DefaultAssay(seurat_human) <- "SCT"
markers.to.plot<-c("TMC1","PCP4",
                   "GJB2","GATA3","TSC22D1","EPYC","OTOG","OTOGL",
                   "EMILIN2","RARRES1",
                   "NEFL","SNAP25",
                   "MBP","PLP1",
                   "MPZ","PMP22",
                   "DCT","TYR",
                   "DCLK1","ESRRB",
                   "TJP1","CLDN11",
                   "ITGA8",	"CREB5", 
                   "VWF","MECOM",
                   "VEPH1","CHST9",
                   "COCH","SLC4A10",
                   "SLC7A11","SLC13A3"
                   ,"ATP2B2","SLC26A7",
                   "CEMIP","THSD4",
                   "RUNX2","COL1A2",
                   "CD163","MEF2C",
                   "IGFBP7","COL4A1",
                   "HBG2","HBM",
                   "TOP2A","HMGB2")
DotPlot(seurat_human, features = markers.to.plot, dot.scale = 8,
        cols  =c("white", "#ad9300")) + RotatedAxis()
a <- table(seurat_human$subclasses_label,seurat_human$seurat_clusters)
a <- as.data.frame.array(a)
write.csv(a,file = "cell_proportion_of_human_cluster_label.csv")
#rename subclasses
Idents(seurat_human) <- "seurat_clusters"
table(Idents(seurat_human))
new.cluster.ids <- c("SMC","TBC","BC","GC","CEC","BC","GC",
                     "CEC","GC","IC","Fibro1","Fibro1","Fibro1",
                     "SMC","GC","IC","CEC","SMC","CEC","GC","Fibro2",
                     "TBC","SMC","CEC","Fibro1","CEC","Fibro1","BC",
                     "Fibro2","Fibro1","ECs","PVM/M","DC_PC_HeC","BC",
                     "Immun","IPh_IBC","Fibro2","Fibro4","SGN","EBs","Fibro1",
                     "NPCs","Fibro3","Fibro1","CEC","SMC","IC","MC","Fib",
                     "DC_PC_HeC","DC_PC_HeC","BC","SGN","SGN","HC"
                     
                     
                     
                     )
names(new.cluster.ids) <- levels(seurat_human)
seurat_human<- RenameIdents(seurat_human, new.cluster.ids)
DimPlot(seurat_human, reduction = "umap",label = TRUE,repel = TRUE)
seurat_human$subclasses_label <- Idents(seurat_human)
p1=DimPlot(seurat_human, reduction = "umap",label = TRUE,repel = TRUE,group.by = "subclasses_label")
p2=DimPlot(seurat_human, reduction = "umap",label = TRUE,repel = TRUE,group.by = "subclasses")
p1+p2
#rename subtype
Idents(seurat_human) <- "seurat_clusters"
table(Idents(seurat_human))
new.cluster.ids <- c("SMC_1","TBC_1","BC_1","GC_1","CEC_1","BC_2","GC_2",
                     "CEC_2","GC_3","IC_1","Fibro1_1","Fibro1_2","Fibro1_3",
                     "SMC_2","GC_4","IC_2","CEC_3","SMC_3","CEC_4","GC_5","Fibro2_1",
                     "TBC_2","SMC_4","CEC_5","Fibro1_4","CEC_6","Fibro1_5","BC_3",
                     "Fibro2_2","Fibro1_6","ECs","PVM/M","DC_PC_HeC_1","BC_4",
                     "Immun","IPh_IBC","Fibro2_3","Fibro4","SGN_1","EBs","Fibro1_7",
                     "NPCs","Fibro3","Fibro1_8","CEC_7","SMC_5","IC_3","MC","Fib",
                     "DC_PC_HeC_2","DC_PC_HeC_3","BC_5","SGN_2","SGN_3","HC")
names(new.cluster.ids) <- levels(seurat_human)
seurat_human<- RenameIdents(seurat_human, new.cluster.ids)
DimPlot(seurat_human, reduction = "umap",label = TRUE,repel = TRUE)
seurat_human$cluster_label <- Idents(seurat_human)
p1=DimPlot(seurat_human, reduction = "umap",label = TRUE,repel = TRUE,group.by = "cluster_label")
p2=DimPlot(seurat_human, reduction = "umap",label = TRUE,repel = TRUE,group.by = "subclasses_label")
p1+p2
saveRDS(seurat_human,file = "seurat_human_new.rds")
#####################################################################################################

DefaultAssay(seurat_human) <- "SCT"
p5 <- DimPlot(seurat_human, reduction = "umap",label = TRUE,group.by = "subclasses_label")
p6 <- DimPlot(seurat_human, reduction = "umap", label = TRUE,group.by = "seurat_clusters")
p7 <- DimPlot(seurat_human, reduction = "umap", label = TRUE,group.by = "cluster_label")
p5|p6|p7


Idents(seurat_human) <- factor(Idents(seurat_human), 
                                     levels = c( "HC","DC_PC_HeC_1","DC_PC_HeC_2","DC_PC_HeC_3",
                                                 "IPh_IBC",
                                                 "TBC_1","TBC_2",
                                                 "SGN_1","SGN_2","SGN_3",
                                                 "GC_1","GC_2","GC_3","GC_4","GC_5",
                                                 "IC_1","IC_2","IC_3",
                                                 "MC","BC_1","BC_2","BC_3","BC_4","BC_5",
                                                 "PVM/M",
                                                 "CEC_1","CEC_2","CEC_3","CEC_4","CEC_5","CEC_6","CEC_7",
                                                 "Fib","Fibro1_1","Fibro1_2","Fibro1_3","Fibro1_4","Fibro1_5","Fibro1_6","Fibro1_7","Fibro1_8",
                                                 "Fibro2_1","Fibro2_2","Fibro2_3",
                                                 "Fibro3","Fibro4",
                                                 "SMC_1","SMC_2","SMC_3", "SMC_4","SMC_5",
                                                 "Immun","ECs", "EBs","NPCs"
                                                 
                                     ))
markers.to.plot<-c("TMC1","PCP4",
                   "GJB2","GATA3","TSC22D1","EPYC","OTOG","OTOGL",
                   "EMILIN2","RARRES1",
                   "NEFL","SNAP25",
                   "MBP","PLP1",
                   "MPZ","PMP22",
                   "DCT","TYR",
                   "DCLK1","ESRRB",
                   "TJP1","CLDN11",
                   "ITGA8",	"CREB5", 
                   "VWF","MECOM",
                   "VEPH1","CHST9",
                   "COCH","SLC4A10",
                   "SLC7A11","SLC13A3"
                   ,"ATP2B2","SLC26A7",
                   "CEMIP","THSD4",
                   "RUNX2","COL1A2",
                   "CD163","MEF2C",
                   "IGFBP7","COL4A1",
                   "HBG2","HBM",
                   "TOP2A","HMGB2"
                   
)
DotPlot(seurat_human, features = markers.to.plot, dot.scale = 8,
        cols  =c("white", "#ad9300")) + RotatedAxis()
dev.off()
#####################################################################################
#####################################################################################
################################human cochlear clusters were used for hiearchical clustering
### cluster_label were used
#####################################################################################
#####################################################################################
seurat_human_new <- readRDS("C:/Users/Dell/Desktop/WORKPLACE/R/seurat_human_new.rds")
table(seurat_human_new$cluster_label,seurat_human_new$orig.ident)
library(pvclust)
library(parallel)
library(dendextend)
library(future)
plan("multisession", workers = 12)
options(future.globals.maxSize= 62914560000) # 60GB (60000*1024^2)

library(Seurat)
library(dplyr)
library(cowplot)
library(zoo)
library(ggplot2)
library(scater)
library(stringr)
library(gplots)
library(matrixStats)
DefaultAssay(seurat_human) <- "SCT"
neurons_avg <- AverageExpression(seurat_human)
neurons_avg_data <- neurons_avg[["SCT"]]
#gene filtering at least 20% expressed in any one cluster.
NoisyGenes <- function(object, min.pct=0.2, clusters.use) {
  clusters <- clusters.use
  genes.use = rownames(object)
  object.raw.data <- as.matrix(GetAssayData(object, slot="counts"))
  pct.matrix = matrix(data=NA, nrow=length(genes.use), ncol=length(clusters))
  rownames(pct.matrix) <- genes.use
  colnames(pct.matrix) <- clusters
  thresh.min=0
  for (i in clusters){
    cells.cluster <- WhichCells(object=object, idents=i)
    data.cluster <- object.raw.data[,colnames(object.raw.data) %in% cells.cluster]
    pct.cluster <- round(apply(object.raw.data[genes.use, cells.cluster, drop = F],1,function(x)return(length(x[x>thresh.min])/length(x))),3)
    pct.matrix[,i] <- pct.cluster
  }
  pct.max <- rowMaxs(pct.matrix)
  names(pct.max) <- genes.use
  noisy.genes <- names(pct.max[pct.max < min.pct])
  return(noisy.genes)
}
# identify noisy genes. the third argument of the function are the clusters. 
#remember to check that the Idents in Seurat are the clusters
noisy.liz <- NoisyGenes(seurat_human, 0.5, levels(Idents(seurat_human)))#0.4;0.5

neurons_avg_data_filtered <- neurons_avg_data[!rownames(neurons_avg_data) %in% noisy.liz,]

dist_func <- function(x){x <- as.dist(1-cor(x, method="s"))}

pvclust.avg_dist <- pvclust(neurons_avg_data_filtered, 
                            method.dist=dist_func, method.hclust="ward.D2", nboot=10000, parallel=T)

plot(pvclust.avg_dist)
pvrect(pvclust.avg_dist, alpha=0.90)
dev.off()

#plotting dendrogram colore by significance
dend <- as.dendrogram(pvclust.avg_dist)
pdf("../plotting/plots/pvclust.TFs.dist_dendextend_10000_black_au_new.pdf", width=60, heigh=45)
dend %>% pvclust_show_signif_gradient(pvclust.avg_dist, signif_col_fun = colorRampPalette(c("lightgrey","grey","black")), signif_type="au") %>%
  plot(main = "neuronal clusters branch confidence")
dev.off()

dendrogram_names <- labels(dend)
dend <- as.dendrogram(pvclust.avg_dist)
dend %>%
  pvclust_show_signif(pvclust.avg_dist, signif_value = c("black", "grey"), show_type = "col") %>%
  plot(main = "Cluster dendrogram with AU/BP values (%)")
# pvrect(result, alpha=0.95)
pvrect2(pvclust.avg_dist, alpha=0.90)
dev.off()
plot(rotate(dend),
     main =
       "Rotated"
)
dev.off()
###make barplots according to cluster ordering, divided by orig.ident and region
library(dittoSeq)
dittoBarPlot(object = seurat_human_new,var = "orig.ident",group.by = "cluster_label")
Idents(seurat_human_new) <- factor(seurat_human_new$cluster_label,levels =
                                          c("HC","Fibro3","MC","Fib","DC_PC_HeC",
                                            "IPh_IBC","CEC_4","CEC_6",
                                            "CEC_1","CEC_2","IC_1","IC_2","IC_3",
                                            "Immun","NPCs","GC_2","GC_1",
                                            "GC_4","GC_3","CEC_5","CEC_3","EBs","SGN",
                                            "SMC_4","ECs","SMC_3","SMC_1","SMC_2",
                                            "Fibro4","Fibro2_1","Fibro2_2","BC_1",
                                            "PVM/M","TBC","Fibro1_3","Fibro1_4",
                                            "BC_2","Fibro1_5","Fibro1_1","Fibro1_2"))
table(Idents(seurat_human_new))
seurat_human_new$cluster_label2 = Idents(seurat_human_new)

dittoBarPlot(object = seurat_human_new,var = "orig.ident",group.by = "cluster_label2")

Idents(seurat_human_new) <- "subclasses_label"
table(Idents(seurat_human_new))
#Stria Vascularis; Modiolus;Spiral Ligament; Cochlear Epithelium;
new.cluster.ids <- c("Stria Vascularis","Spiral Ligament","Modiolus","Cochlear Epithelium",
                     "Spiral Ligament","Spiral Ligament","Stria Vascularis","Stria Vascularis",
                     "Cochlear Epithelium","Modiolus","Others","Stria Vascularis","Others",
                     "Cochlear Epithelium","Others","Spiral Ligament","Spiral Ligament",
                     "Others","Stria Vascularis","Spiral Ligament","Cochlear Epithelium"
                     )
names(new.cluster.ids) <- levels(seurat_human_new)
seurat_human_new<- RenameIdents(seurat_human_new, new.cluster.ids)
DimPlot(seurat_human_new, reduction = "umap",label = TRUE,repel = TRUE)
seurat_human_new$classes <- Idents(seurat_human_new)
dittoBarPlot(object = seurat_human_new,var = "classes",group.by = "cluster_label2")
dev.off()


########################################################################
#####################################################################################
###########cumulative plot and dendrogram based on TFs##############################
#####################################################################################
#####################################################################################
seurat_human_new <- readRDS("C:/Users/Dell/Desktop/WORKPLACE/R/seurat_human_new.rds")
seurat_mouse_new <- readRDS("C:/Users/Dell/Desktop/WORKPLACE/R/seurat_mouse_new.rds")

##for mouse species

library(Seurat)
library(matrixStats)
library(ggplot2)

Idents(seurat_mouse_new) <- "cluster_label"
DimPlot(seurat_mouse_new, reduction = "umap",label = TRUE)
DefaultAssay(seurat_mouse_new) <- "SCT"
object <- seurat_mouse_new
table(Idents(seurat_mouse_new))
clusters <-levels(Idents(seurat_mouse_new))
genes.use = rownames(object)
object.raw.data <- as.matrix(GetAssayData(object, slot="data",assay = "SCT"))
pct.matrix = matrix(data=NA, nrow=length(genes.use), ncol=length(clusters))
rownames(pct.matrix) <- genes.use
colnames(pct.matrix) <- clusters
thresh.min=0
for (i in clusters){
  cells.cluster <- WhichCells(object=object, idents=i)
  data.cluster <- object.raw.data[,colnames(object.raw.data) %in% cells.cluster]
  pct.cluster <- round(apply(object.raw.data[genes.use, cells.cluster, drop = F],1,function(x)return(length(x[x>thresh.min])/length(x))),3)
  pct.matrix[,i] <- pct.cluster
}
pct.df <- as.data.frame(pct.matrix)
Human_TFs <- read.csv("Human_TFs.csv", header=T, row.names = 1)
Human_TFs$gene_family <- "transcription factors"
Human_TFs <- Human_TFs[,c(1,3,2)]
names(Human_TFs) <- c("ncbi_gene_symbol","gene_family","TF_family")

Human_GPCRs <- read.delim("GPCRs_gene_symbols.txt", header=T, sep = "\t")
head(Human_GPCRs)
Human_GPCRs <- Human_GPCRs[,c(1:2)]
Human_GPCRs$gene_family <- "GPCR"
Human_GPCRs$TF_family <- NA

Human_ICs <- read.delim("ion-channels_gene_symbols.txt", header=T, sep = "\t")
head(Human_ICs)
Human_ICs <- Human_ICs[,c(1:2)]
Human_ICs$gene_family <- "ion channels"
Human_ICs$TF_family <- NA

Human_CAMs <- read.delim("CAMs_gene_symbols.txt", header=T, sep = "\t")
head(Human_CAMs)
Human_CAMs <- Human_CAMs[,c(1:2)]
Human_CAMs$gene_family <- "cell adhesion molecules"
Human_CAMs$TF_family <- NA

Human_PGs <- read.delim("proteoglycans_gene_symbols.txt", header=T, sep = "\t")
head(Human_PGs)
Human_PGs <- Human_PGs[,c(1:2)]
Human_PGs$gene_family <- "proteoglycans"
Human_PGs$TF_family <- NA

Human_ribo <- read.delim("ribosome_gene_symbols.txt", header=T, sep = "\t")
head(Human_ribo)
Human_ribo <- Human_ribo[,c(1:2)]
Human_ribo$gene_family <- "ribosome"
Human_ribo$TF_family <- NA


#transcription factors
TF_intersection <- intersect(rownames(pct.df),Human_TFs$ncbi_gene_symbol)
pct.df_TFs <- pct.df[rownames(pct.df) %in% TF_intersection,]
Human_TFs_cropped <- Human_TFs[Human_TFs$ncbi_gene_symbol %in% TF_intersection,]

Human_TFs_cluster_counts <- matrix(data=NA,nrow=length(rownames(pct.df_TFs)),ncol=2)
Human_TFs_cluster_counts <- as.data.frame(Human_TFs_cluster_counts)
names(Human_TFs_cluster_counts) <- c("gene","cluster_count")
for (i in 1:length(rownames(pct.df_TFs))){
  gene_row <- pct.df_TFs[i,]
  cluster_count <- length(gene_row[gene_row > 0.2])
  Human_TFs_cluster_counts$gene[i] <- rownames(gene_row)
  Human_TFs_cluster_counts$cluster_count[i] <- cluster_count
}

Human_TFs_cluster_counts_matched <- Human_TFs_cropped[order(match(Human_TFs_cropped$ncbi_gene_symbol, 
                                                                  Human_TFs_cluster_counts$gene)),]
Human_TFs_cluster_counts$gene_family <- Human_TFs_cluster_counts_matched$gene_family
Human_TFs_cluster_counts$TF_family <- Human_TFs_cluster_counts_matched$TF_family
TF_counts <- Human_TFs_cluster_counts


#GPCRs
GPCR_intersection <- intersect(rownames(pct.df),Human_GPCRs$ncbi_gene_symbol)
pct.df_GPCRs <- pct.df[rownames(pct.df) %in% GPCR_intersection,]
Human_GPCRs_cropped <- Human_GPCRs[Human_GPCRs$ncbi_gene_symbol %in% GPCR_intersection,]

Human_GPCRs_cluster_counts <- matrix(data=NA,nrow=length(rownames(pct.df_GPCRs)),ncol=2)
Human_GPCRs_cluster_counts <- as.data.frame(Human_GPCRs_cluster_counts)
names(Human_GPCRs_cluster_counts) <- c("gene","cluster_count")

for (i in 1:length(rownames(pct.df_GPCRs))){
  gene_row <- pct.df_GPCRs[i,]
  cluster_count <- length(gene_row[gene_row > 0.05])
  Human_GPCRs_cluster_counts$gene[i] <- rownames(gene_row)
  Human_GPCRs_cluster_counts$cluster_count[i] <- cluster_count
}
Human_GPCRs_cluster_counts_matched <- Human_GPCRs_cropped[order(match(Human_GPCRs_cropped$ncbi_gene_symbol, 
                                                                      Human_GPCRs_cluster_counts$gene)),]
Human_GPCRs_cluster_counts$gene_family <- Human_GPCRs_cluster_counts_matched$gene_family
GPCR_counts <- Human_GPCRs_cluster_counts


#ion channels
IC_intersection <- intersect(rownames(pct.df),Human_ICs$ncbi_gene_symbol)
pct.df_ICs <- pct.df[rownames(pct.df) %in% IC_intersection,]
Human_ICs_cropped <- Human_ICs[Human_ICs$ncbi_gene_symbol %in% IC_intersection,]

Human_ICs_cluster_counts <- matrix(data=NA,nrow=length(rownames(pct.df_ICs)),ncol=2)
Human_ICs_cluster_counts <- as.data.frame(Human_ICs_cluster_counts)
names(Human_ICs_cluster_counts) <- c("gene","cluster_count")

for (i in 1:length(rownames(pct.df_ICs))){
  gene_row <- pct.df_ICs[i,]
  cluster_count <- length(gene_row[gene_row > 0.1])
  Human_ICs_cluster_counts$gene[i] <- rownames(gene_row)
  Human_ICs_cluster_counts$cluster_count[i] <- cluster_count
}


Human_ICs_cluster_counts_matched <- Human_ICs_cropped[order(match(Human_ICs_cropped$ncbi_gene_symbol,
                                                                  Human_ICs_cluster_counts$gene)),]
Human_ICs_cluster_counts$gene_family <- Human_ICs_cluster_counts_matched$gene_family
IC_counts <- Human_ICs_cluster_counts

#cell adhesion molecules
CAM_intersection <- intersect(rownames(pct.df),Human_CAMs$ncbi_gene_symbol)
pct.df_CAMs <- pct.df[rownames(pct.df) %in% CAM_intersection,]
Human_CAMs_cropped <- Human_CAMs[Human_CAMs$ncbi_gene_symbol %in% CAM_intersection,]

Human_CAMs_cluster_counts <- matrix(data=NA,nrow=length(rownames(pct.df_CAMs)),ncol=2)
Human_CAMs_cluster_counts <- as.data.frame(Human_CAMs_cluster_counts)
names(Human_CAMs_cluster_counts) <- c("gene","cluster_count")

for (i in 1:length(rownames(pct.df_CAMs))){
  gene_row <- pct.df_CAMs[i,]
  cluster_count <- length(gene_row[gene_row > 0.05])
  Human_CAMs_cluster_counts$gene[i] <- rownames(gene_row)
  Human_CAMs_cluster_counts$cluster_count[i] <- cluster_count
}

Human_CAMs_cluster_counts_matched <- Human_CAMs_cropped[order(match(Human_CAMs_cropped$ncbi_gene_symbol, 
                                                                    Human_CAMs_cluster_counts$gene)),]
Human_CAMs_cluster_counts$gene_family <- Human_CAMs_cluster_counts_matched$gene_family
CAM_counts <- Human_CAMs_cluster_counts
#proteoglycans
PG_intersection <- intersect(rownames(pct.df),Human_PGs$ncbi_gene_symbol)
pct.df_PGs <- pct.df[rownames(pct.df) %in% PG_intersection,]
Human_PGs_cropped <- Human_PGs[Human_PGs$ncbi_gene_symbol %in% PG_intersection,]

Human_PGs_cluster_counts <- matrix(data=NA,nrow=length(rownames(pct.df_PGs)),ncol=2)
Human_PGs_cluster_counts <- as.data.frame(Human_PGs_cluster_counts)
names(Human_PGs_cluster_counts) <- c("gene","cluster_count")
for (i in 1:length(rownames(pct.df_PGs))){
  gene_row <- pct.df_PGs[i,]
  cluster_count <- length(gene_row[gene_row > 0.05])
  Human_PGs_cluster_counts$gene[i] <- rownames(gene_row)
  Human_PGs_cluster_counts$cluster_count[i] <- cluster_count
}

Human_PGs_cluster_counts_matched <- Human_PGs_cropped[order(match(Human_PGs_cropped$ncbi_gene_symbol, 
                                                                  Human_PGs_cluster_counts$gene)),]
Human_PGs_cluster_counts$gene_family <- Human_PGs_cluster_counts_matched$gene_family
PG_counts <- Human_PGs_cluster_counts


#ribosomes
ribo_intersection <- intersect(rownames(pct.df),Human_ribo$ncbi_gene_symbol)
pct.df_ribos <- pct.df[rownames(pct.df) %in% ribo_intersection,]
Human_ribos_cropped <- Human_ribo[Human_ribo$ncbi_gene_symbol %in% ribo_intersection,]
Human_ribos_cropped <- Human_ribos_cropped[unique(Human_ribos_cropped$ncbi_gene_symbol),]

Human_ribos_cluster_counts <- matrix(data=NA,nrow=length(rownames(pct.df_ribos)),ncol=2)
Human_ribos_cluster_counts <- as.data.frame(Human_ribos_cluster_counts)
names(Human_ribos_cluster_counts) <- c("gene","cluster_count")

for (i in 1:length(rownames(pct.df_ribos))){
  gene_row <- pct.df_ribos[i,]
  cluster_count <- length(gene_row[gene_row > 0.025])
  Human_ribos_cluster_counts$gene[i] <- rownames(gene_row)
  Human_ribos_cluster_counts$cluster_count[i] <- cluster_count
}


Human_ribos_cluster_counts_matched <- Human_ribos_cropped[order(match(Human_ribos_cropped$ncbi_gene_symbol, 
                                                                      Human_ribos_cluster_counts$gene)),]
Human_ribos_cluster_counts$gene_family <- Human_ribos_cluster_counts_matched$gene_family
Human_ribos_cluster_counts$ribo_family <- Human_ribos_cluster_counts_matched$ribo_family
ribo_counts <- Human_ribos_cluster_counts


#Homeodomain 
countsHD <- TF_counts[TF_counts$TF_family=="Homeodomain",]
countsHD <- countsHD[order(countsHD$cluster_count, decreasing = F),]
countsHD <- countsHD[countsHD$cluster_count>0,]
tHD <- table(countsHD$cluster_count)
xHD <- as.numeric(names(tHD))
yHD <- as.numeric(tHD)
plot(xHD,(cumsum(yHD)/max(cumsum(yHD))), cex=0.5, type="o", pch=19)
#bHLH 
countsbHLH <- TF_counts[TF_counts$TF_family=="bHLH",]
countsbHLH <- countsbHLH[order(countsbHLH$cluster_count, decreasing = F),]
countsbHLH <- countsbHLH[countsbHLH$cluster_count>0,]
tbHLH <- table(countsbHLH$cluster_count)
xbHLH <- as.numeric(names(tbHLH))
ybHLH <- as.numeric(tbHLH)
plot(xbHLH,(cumsum(ybHLH)/max(cumsum(ybHLH))), cex=0.5, type="o", pch=19)
#bZIP 
countsZIP <- TF_counts[TF_counts$TF_family=="bZIP",]
countsZIP <- countsZIP[order(countsZIP$cluster_count, decreasing = F),]
countsZIP <- countsZIP[countsZIP$cluster_count>0,]
tZIP <- table(countsZIP$cluster_count)
xZIP <- as.numeric(names(tZIP))
yZIP <- as.numeric(tZIP)
plot(xZIP,(cumsum(yZIP)/max(cumsum(yZIP))), cex=0.5, type="o", pch=19)

#Forkhead 
countsFH <- TF_counts[TF_counts$TF_family=="Forkhead",]
countsFH <- countsFH[order(countsFH$cluster_count, decreasing = F),]
countsFH <- countsFH[countsFH$cluster_count>0,]
tFH <- table(countsFH$cluster_count)
xFH <- as.numeric(names(tFH))
yFH <- as.numeric(tFH)
plot(xFH,(cumsum(yFH)/max(cumsum(yFH))), cex=0.5, type="o", pch=19)

#all TFs 
counts <- TF_counts[order(TF_counts$cluster_count, decreasing = F),]
countsTFs <- counts[counts$cluster_count>0,]
tTFs <- table(countsTFs$cluster_count)
xTFs <- as.numeric(names(tTFs))
yTFs <- as.numeric(tTFs)
plot(xTFs,(cumsum(yTFs)/max(cumsum(yTFs))), cex=0.5, type="o", pch=19)

#Nuclear receptor 
countsNR <- TF_counts[TF_counts$TF_family=="Nuclear receptor",]
countsNR <- countsNR[order(countsNR$cluster_count, decreasing = F),]
countsNR <- countsNR[countsNR$cluster_count>0,]
tNR <- table(countsNR$cluster_count)
xNR <- as.numeric(names(tNR))
yNR <- as.numeric(tNR)
plot(xNR,(cumsum(yNR)/max(cumsum(yNR))), cex=0.5, type="o", pch=19)

#C2H2 ZF 
countsC2F <- TF_counts[TF_counts$TF_family=="C2H2 ZF",]
countsC2F <- countsC2F[order(countsC2F$cluster_count, decreasing = F),]
countsC2F <- countsC2F[countsC2F$cluster_count>0,]
tC2F <- table(countsC2F$cluster_count)
xC2F <- as.numeric(names(tC2F))
yC2F <- as.numeric(tC2F)
plot(xC2F,(cumsum(yC2F)/max(cumsum(yC2F))), cex=0.5, type="o", pch=19)

plot(xTFs,(cumsum(yTFs)/max(cumsum(yTFs))), cex=0.5, type="o", pch=19, xlab = "Number of cochlear epithelium cell types expressing", ylab = "Cumulative fraction")
lines(xNR,(cumsum(yNR)/max(cumsum(yNR))), cex=0.5, type="o", pch=19, col="magenta")
lines(xC2F,(cumsum(yC2F)/max(cumsum(yC2F))), cex=0.5, type="o", pch=19, col="yellow")
lines(xHD,(cumsum(yHD)/max(cumsum(yHD))), cex=0.5, type="o", pch=19, col="darkgreen")
lines(xbHLH,(cumsum(ybHLH)/max(cumsum(ybHLH))), cex=0.5, type="o", pch=19, col="blue")
lines(xFH,(cumsum(yFH)/max(cumsum(yFH))), cex=0.5, type="o", pch=19, col="red")
lines(xZIP,(cumsum(yZIP)/max(cumsum(yZIP))), cex=0.5, type="o", pch=19, col="orange")
#all ion channels 
counts_ICs <- IC_counts[order(IC_counts$cluster_count, decreasing = F),]
countsICs <- counts_ICs[counts_ICs$cluster_count>0,]
tICs <- table(countsICs$cluster_count)
xICs <- as.numeric(names(tICs))
yICs <- as.numeric(tICs)
plot(xICs,(cumsum(yICs)/max(cumsum(yICs))), cex=0.5, type="o", pch=19)

#all GPCRs 
counts_GPCRs <- GPCR_counts[order(GPCR_counts$cluster_count, decreasing = F),]
countsGPCRs <- counts_GPCRs[counts_GPCRs$cluster_count>0,]
tGPCRs <- table(countsGPCRs$cluster_count)
xGPCRs <- as.numeric(names(tGPCRs))
yGPCRs <- as.numeric(tGPCRs)
plot(xGPCRs,(cumsum(yGPCRs)/max(cumsum(yGPCRs))), cex=0.5, type="o", pch=19)

#all CAMs 
counts_CAMs <- CAM_counts[order(CAM_counts$cluster_count, decreasing = F),]
countsCAMs <- counts_CAMs[counts_CAMs$cluster_count>0,]
tCAMs <- table(countsCAMs$cluster_count)
xCAMs <- as.numeric(names(tCAMs))
yCAMs <- as.numeric(tCAMs)
plot(xCAMs,(cumsum(yCAMs)/max(cumsum(yCAMs))), cex=0.5, type="o", pch=19)

#all proteoglycans
counts_PGs <- PG_counts[order(PG_counts$cluster_count, decreasing = F),]
countsPGs <- counts_PGs[counts_PGs$cluster_count>0,]
tPGs <- table(countsPGs$cluster_count)
xPGs <- as.numeric(names(tPGs))
yPGs <- as.numeric(tPGs)
plot(xPGs,(cumsum(yPGs)/max(cumsum(yPGs))), cex=0.5, type="o", pch=19)
#all ribosomal genes 
counts_ribos <- ribo_counts[order(ribo_counts$cluster_count, decreasing = F),]
countsribos <- counts_ribos[counts_ribos$cluster_count>0,]
tRibos <- table(countsribos$cluster_count)
xRibos <- as.numeric(names(tRibos))
yRibos <- as.numeric(tRibos)
plot(xRibos,(cumsum(yRibos)/max(cumsum(yRibos))), cex=0.8, type="o", pch=19)
plot(xTFs,(cumsum(yTFs)/max(cumsum(yTFs))), cex=0.8, type="o", pch=19, 
     xlab = "Number of neuron types expressing",
     ylab = "Cumulative fraction",xlim=c(0,40),ylim=c(0,1), col="black")
lines(xC2F,(cumsum(yC2F)/max(cumsum(yC2F))), cex=0.8, type="o", pch=19, col="magenta")
lines(xHD,(cumsum(yHD)/max(cumsum(yHD))), cex=0.8, type="o", pch=19, col="darkgreen")
lines(xGPCRs,(cumsum(yGPCRs)/max(cumsum(yGPCRs))), cex=0.8, type="o", pch=19, col="blue")
lines(xICs,(cumsum(yICs)/max(cumsum(yICs))), cex=0.8, type="o", pch=19, col="red")
lines(xCAMs,(cumsum(yCAMs)/max(cumsum(yCAMs))), cex=0.8, type="o", pch=19, col="orange")
lines(xPGs,(cumsum(yPGs)/max(cumsum(yPGs))), cex=0.8, type="o", pch=19, col="darkblue")
lines(xRibos,(cumsum(yRibos)/max(cumsum(yRibos))), cex=0.8, type="o", pch=19, col="darkred")

legend("bottomright", legend=c("all TFs ", "zinc finger-C2H2",
                               "homeodomain ","GPCRs ",
                               "ion channels ","cell adhesion molecules",
                               "proteoglycans ","ribosomal genes "),
       col=c("black", "magenta","darkgreen","blue","red","orange","darkblue","darkred"), lty=2, cex=0.5)




dev.off()
############################################################################################

##for human species

library(Seurat)
library(matrixStats)
library(ggplot2)

Idents(seurat_human) <- "cluster_label"
DimPlot(seurat_human, reduction = "umap",label = TRUE)
DefaultAssay(seurat_human) <- "SCT"
object <- seurat_human
table(Idents(seurat_human))
clusters <-levels(Idents(seurat_human))
genes.use = rownames(object)
object.raw.data <- as.matrix(GetAssayData(object, slot="data",assay = "SCT"))
pct.matrix = matrix(data=NA, nrow=length(genes.use), ncol=length(clusters))
rownames(pct.matrix) <- genes.use
colnames(pct.matrix) <- clusters
thresh.min=0
for (i in clusters){
  cells.cluster <- WhichCells(object=object, idents=i)
  data.cluster <- object.raw.data[,colnames(object.raw.data) %in% cells.cluster]
  pct.cluster <- round(apply(object.raw.data[genes.use, cells.cluster, drop = F],1,function(x)return(length(x[x>thresh.min])/length(x))),3)
  pct.matrix[,i] <- pct.cluster
}
pct.df <- as.data.frame(pct.matrix)
Human_TFs <- read.csv("Human_TFs.csv", header=T, row.names = 1)
Human_TFs$gene_family <- "transcription factors"
Human_TFs <- Human_TFs[,c(1,3,2)]
names(Human_TFs) <- c("ncbi_gene_symbol","gene_family","TF_family")

Human_GPCRs <- read.delim("GPCRs_gene_symbols.txt", header=T, sep = "\t")
head(Human_GPCRs)
Human_GPCRs <- Human_GPCRs[,c(1:2)]
Human_GPCRs$gene_family <- "GPCR"
Human_GPCRs$TF_family <- NA

Human_ICs <- read.delim("ion-channels_gene_symbols.txt", header=T, sep = "\t")
head(Human_ICs)
Human_ICs <- Human_ICs[,c(1:2)]
Human_ICs$gene_family <- "ion channels"
Human_ICs$TF_family <- NA

Human_CAMs <- read.delim("CAMs_gene_symbols.txt", header=T, sep = "\t")
head(Human_CAMs)
Human_CAMs <- Human_CAMs[,c(1:2)]
Human_CAMs$gene_family <- "cell adhesion molecules"
Human_CAMs$TF_family <- NA

Human_PGs <- read.delim("proteoglycans_gene_symbols.txt", header=T, sep = "\t")
head(Human_PGs)
Human_PGs <- Human_PGs[,c(1:2)]
Human_PGs$gene_family <- "proteoglycans"
Human_PGs$TF_family <- NA

Human_ribo <- read.delim("ribosome_gene_symbols.txt", header=T, sep = "\t")
head(Human_ribo)
Human_ribo <- Human_ribo[,c(1:2)]
Human_ribo$gene_family <- "ribosome"
Human_ribo$TF_family <- NA


#transcription factors
TF_intersection <- intersect(rownames(pct.df),Human_TFs$ncbi_gene_symbol)
pct.df_TFs <- pct.df[rownames(pct.df) %in% TF_intersection,]
Human_TFs_cropped <- Human_TFs[Human_TFs$ncbi_gene_symbol %in% TF_intersection,]

Human_TFs_cluster_counts <- matrix(data=NA,nrow=length(rownames(pct.df_TFs)),ncol=2)
Human_TFs_cluster_counts <- as.data.frame(Human_TFs_cluster_counts)
names(Human_TFs_cluster_counts) <- c("gene","cluster_count")
for (i in 1:length(rownames(pct.df_TFs))){
  gene_row <- pct.df_TFs[i,]
  cluster_count <- length(gene_row[gene_row > 0.2])
  Human_TFs_cluster_counts$gene[i] <- rownames(gene_row)
  Human_TFs_cluster_counts$cluster_count[i] <- cluster_count
}

Human_TFs_cluster_counts_matched <- Human_TFs_cropped[order(match(Human_TFs_cropped$ncbi_gene_symbol, 
                                                                  Human_TFs_cluster_counts$gene)),]
Human_TFs_cluster_counts$gene_family <- Human_TFs_cluster_counts_matched$gene_family
Human_TFs_cluster_counts$TF_family <- Human_TFs_cluster_counts_matched$TF_family
TF_counts <- Human_TFs_cluster_counts


#GPCRs
GPCR_intersection <- intersect(rownames(pct.df),Human_GPCRs$ncbi_gene_symbol)
pct.df_GPCRs <- pct.df[rownames(pct.df) %in% GPCR_intersection,]
Human_GPCRs_cropped <- Human_GPCRs[Human_GPCRs$ncbi_gene_symbol %in% GPCR_intersection,]

Human_GPCRs_cluster_counts <- matrix(data=NA,nrow=length(rownames(pct.df_GPCRs)),ncol=2)
Human_GPCRs_cluster_counts <- as.data.frame(Human_GPCRs_cluster_counts)
names(Human_GPCRs_cluster_counts) <- c("gene","cluster_count")

for (i in 1:length(rownames(pct.df_GPCRs))){
  gene_row <- pct.df_GPCRs[i,]
  cluster_count <- length(gene_row[gene_row > 0.05])
  Human_GPCRs_cluster_counts$gene[i] <- rownames(gene_row)
  Human_GPCRs_cluster_counts$cluster_count[i] <- cluster_count
}
Human_GPCRs_cluster_counts_matched <- Human_GPCRs_cropped[order(match(Human_GPCRs_cropped$ncbi_gene_symbol, 
                                                                      Human_GPCRs_cluster_counts$gene)),]
Human_GPCRs_cluster_counts$gene_family <- Human_GPCRs_cluster_counts_matched$gene_family
GPCR_counts <- Human_GPCRs_cluster_counts


#ion channels
IC_intersection <- intersect(rownames(pct.df),Human_ICs$ncbi_gene_symbol)
pct.df_ICs <- pct.df[rownames(pct.df) %in% IC_intersection,]
Human_ICs_cropped <- Human_ICs[Human_ICs$ncbi_gene_symbol %in% IC_intersection,]

Human_ICs_cluster_counts <- matrix(data=NA,nrow=length(rownames(pct.df_ICs)),ncol=2)
Human_ICs_cluster_counts <- as.data.frame(Human_ICs_cluster_counts)
names(Human_ICs_cluster_counts) <- c("gene","cluster_count")

for (i in 1:length(rownames(pct.df_ICs))){
  gene_row <- pct.df_ICs[i,]
  cluster_count <- length(gene_row[gene_row > 0.2])
  Human_ICs_cluster_counts$gene[i] <- rownames(gene_row)
  Human_ICs_cluster_counts$cluster_count[i] <- cluster_count
}


Human_ICs_cluster_counts_matched <- Human_ICs_cropped[order(match(Human_ICs_cropped$ncbi_gene_symbol,
                                                                  Human_ICs_cluster_counts$gene)),]
Human_ICs_cluster_counts$gene_family <- Human_ICs_cluster_counts_matched$gene_family
IC_counts <- Human_ICs_cluster_counts

#cell adhesion molecules
CAM_intersection <- intersect(rownames(pct.df),Human_CAMs$ncbi_gene_symbol)
pct.df_CAMs <- pct.df[rownames(pct.df) %in% CAM_intersection,]
Human_CAMs_cropped <- Human_CAMs[Human_CAMs$ncbi_gene_symbol %in% CAM_intersection,]

Human_CAMs_cluster_counts <- matrix(data=NA,nrow=length(rownames(pct.df_CAMs)),ncol=2)
Human_CAMs_cluster_counts <- as.data.frame(Human_CAMs_cluster_counts)
names(Human_CAMs_cluster_counts) <- c("gene","cluster_count")

for (i in 1:length(rownames(pct.df_CAMs))){
  gene_row <- pct.df_CAMs[i,]
  cluster_count <- length(gene_row[gene_row > 0.05])
  Human_CAMs_cluster_counts$gene[i] <- rownames(gene_row)
  Human_CAMs_cluster_counts$cluster_count[i] <- cluster_count
}

Human_CAMs_cluster_counts_matched <- Human_CAMs_cropped[order(match(Human_CAMs_cropped$ncbi_gene_symbol, 
                                                                    Human_CAMs_cluster_counts$gene)),]
Human_CAMs_cluster_counts$gene_family <- Human_CAMs_cluster_counts_matched$gene_family
CAM_counts <- Human_CAMs_cluster_counts
#proteoglycans
PG_intersection <- intersect(rownames(pct.df),Human_PGs$ncbi_gene_symbol)
pct.df_PGs <- pct.df[rownames(pct.df) %in% PG_intersection,]
Human_PGs_cropped <- Human_PGs[Human_PGs$ncbi_gene_symbol %in% PG_intersection,]

Human_PGs_cluster_counts <- matrix(data=NA,nrow=length(rownames(pct.df_PGs)),ncol=2)
Human_PGs_cluster_counts <- as.data.frame(Human_PGs_cluster_counts)
names(Human_PGs_cluster_counts) <- c("gene","cluster_count")
for (i in 1:length(rownames(pct.df_PGs))){
  gene_row <- pct.df_PGs[i,]
  cluster_count <- length(gene_row[gene_row > 0.05])
  Human_PGs_cluster_counts$gene[i] <- rownames(gene_row)
  Human_PGs_cluster_counts$cluster_count[i] <- cluster_count
}

Human_PGs_cluster_counts_matched <- Human_PGs_cropped[order(match(Human_PGs_cropped$ncbi_gene_symbol, 
                                                                  Human_PGs_cluster_counts$gene)),]
Human_PGs_cluster_counts$gene_family <- Human_PGs_cluster_counts_matched$gene_family
PG_counts <- Human_PGs_cluster_counts


#ribosomes
ribo_intersection <- intersect(rownames(pct.df),Human_ribo$ncbi_gene_symbol)
pct.df_ribos <- pct.df[rownames(pct.df) %in% ribo_intersection,]
Human_ribos_cropped <- Human_ribo[Human_ribo$ncbi_gene_symbol %in% ribo_intersection,]
Human_ribos_cropped <- Human_ribos_cropped[unique(Human_ribos_cropped$ncbi_gene_symbol),]

Human_ribos_cluster_counts <- matrix(data=NA,nrow=length(rownames(pct.df_ribos)),ncol=2)
Human_ribos_cluster_counts <- as.data.frame(Human_ribos_cluster_counts)
names(Human_ribos_cluster_counts) <- c("gene","cluster_count")

for (i in 1:length(rownames(pct.df_ribos))){
  gene_row <- pct.df_ribos[i,]
  cluster_count <- length(gene_row[gene_row > 0.1])
  Human_ribos_cluster_counts$gene[i] <- rownames(gene_row)
  Human_ribos_cluster_counts$cluster_count[i] <- cluster_count
}


Human_ribos_cluster_counts_matched <- Human_ribos_cropped[order(match(Human_ribos_cropped$ncbi_gene_symbol, 
                                                                      Human_ribos_cluster_counts$gene)),]
Human_ribos_cluster_counts$gene_family <- Human_ribos_cluster_counts_matched$gene_family
Human_ribos_cluster_counts$ribo_family <- Human_ribos_cluster_counts_matched$ribo_family
ribo_counts <- Human_ribos_cluster_counts


#Homeodomain 
countsHD <- TF_counts[TF_counts$TF_family=="Homeodomain",]
countsHD <- countsHD[order(countsHD$cluster_count, decreasing = F),]
countsHD <- countsHD[countsHD$cluster_count>0,]
tHD <- table(countsHD$cluster_count)
xHD <- as.numeric(names(tHD))
yHD <- as.numeric(tHD)
plot(xHD,(cumsum(yHD)/max(cumsum(yHD))), cex=0.5, type="o", pch=19)
#bHLH 
countsbHLH <- TF_counts[TF_counts$TF_family=="bHLH",]
countsbHLH <- countsbHLH[order(countsbHLH$cluster_count, decreasing = F),]
countsbHLH <- countsbHLH[countsbHLH$cluster_count>0,]
tbHLH <- table(countsbHLH$cluster_count)
xbHLH <- as.numeric(names(tbHLH))
ybHLH <- as.numeric(tbHLH)
plot(xbHLH,(cumsum(ybHLH)/max(cumsum(ybHLH))), cex=0.5, type="o", pch=19)
#bZIP 
countsZIP <- TF_counts[TF_counts$TF_family=="bZIP",]
countsZIP <- countsZIP[order(countsZIP$cluster_count, decreasing = F),]
countsZIP <- countsZIP[countsZIP$cluster_count>0,]
tZIP <- table(countsZIP$cluster_count)
xZIP <- as.numeric(names(tZIP))
yZIP <- as.numeric(tZIP)
plot(xZIP,(cumsum(yZIP)/max(cumsum(yZIP))), cex=0.5, type="o", pch=19)

#Forkhead 
countsFH <- TF_counts[TF_counts$TF_family=="Forkhead",]
countsFH <- countsFH[order(countsFH$cluster_count, decreasing = F),]
countsFH <- countsFH[countsFH$cluster_count>0,]
tFH <- table(countsFH$cluster_count)
xFH <- as.numeric(names(tFH))
yFH <- as.numeric(tFH)
plot(xFH,(cumsum(yFH)/max(cumsum(yFH))), cex=0.5, type="o", pch=19)

#all TFs 
counts <- TF_counts[order(TF_counts$cluster_count, decreasing = F),]
countsTFs <- counts[counts$cluster_count>0,]
tTFs <- table(countsTFs$cluster_count)
xTFs <- as.numeric(names(tTFs))
yTFs <- as.numeric(tTFs)
plot(xTFs,(cumsum(yTFs)/max(cumsum(yTFs))), cex=0.5, type="o", pch=19)

#Nuclear receptor 
countsNR <- TF_counts[TF_counts$TF_family=="Nuclear receptor",]
countsNR <- countsNR[order(countsNR$cluster_count, decreasing = F),]
countsNR <- countsNR[countsNR$cluster_count>0,]
tNR <- table(countsNR$cluster_count)
xNR <- as.numeric(names(tNR))
yNR <- as.numeric(tNR)
plot(xNR,(cumsum(yNR)/max(cumsum(yNR))), cex=0.5, type="o", pch=19)

#C2H2 ZF 
countsC2F <- TF_counts[TF_counts$TF_family=="C2H2 ZF",]
countsC2F <- countsC2F[order(countsC2F$cluster_count, decreasing = F),]
countsC2F <- countsC2F[countsC2F$cluster_count>0,]
tC2F <- table(countsC2F$cluster_count)
xC2F <- as.numeric(names(tC2F))
yC2F <- as.numeric(tC2F)
plot(xC2F,(cumsum(yC2F)/max(cumsum(yC2F))), cex=0.5, type="o", pch=19)

plot(xTFs,(cumsum(yTFs)/max(cumsum(yTFs))), cex=0.5, type="o", pch=19, xlab = "Number of cochlear epithelium cell types expressing", ylab = "Cumulative fraction")
lines(xNR,(cumsum(yNR)/max(cumsum(yNR))), cex=0.5, type="o", pch=19, col="magenta")
lines(xC2F,(cumsum(yC2F)/max(cumsum(yC2F))), cex=0.5, type="o", pch=19, col="yellow")
lines(xHD,(cumsum(yHD)/max(cumsum(yHD))), cex=0.5, type="o", pch=19, col="darkgreen")
lines(xbHLH,(cumsum(ybHLH)/max(cumsum(ybHLH))), cex=0.5, type="o", pch=19, col="blue")
lines(xFH,(cumsum(yFH)/max(cumsum(yFH))), cex=0.5, type="o", pch=19, col="red")
lines(xZIP,(cumsum(yZIP)/max(cumsum(yZIP))), cex=0.5, type="o", pch=19, col="orange")
#all ion channels 
counts_ICs <- IC_counts[order(IC_counts$cluster_count, decreasing = F),]
countsICs <- counts_ICs[counts_ICs$cluster_count>0,]
tICs <- table(countsICs$cluster_count)
xICs <- as.numeric(names(tICs))
yICs <- as.numeric(tICs)
plot(xICs,(cumsum(yICs)/max(cumsum(yICs))), cex=0.5, type="o", pch=19)

#all GPCRs 
counts_GPCRs <- GPCR_counts[order(GPCR_counts$cluster_count, decreasing = F),]
countsGPCRs <- counts_GPCRs[counts_GPCRs$cluster_count>0,]
tGPCRs <- table(countsGPCRs$cluster_count)
xGPCRs <- as.numeric(names(tGPCRs))
yGPCRs <- as.numeric(tGPCRs)
plot(xGPCRs,(cumsum(yGPCRs)/max(cumsum(yGPCRs))), cex=0.5, type="o", pch=19)

#all CAMs 
counts_CAMs <- CAM_counts[order(CAM_counts$cluster_count, decreasing = F),]
countsCAMs <- counts_CAMs[counts_CAMs$cluster_count>0,]
tCAMs <- table(countsCAMs$cluster_count)
xCAMs <- as.numeric(names(tCAMs))
yCAMs <- as.numeric(tCAMs)
plot(xCAMs,(cumsum(yCAMs)/max(cumsum(yCAMs))), cex=0.5, type="o", pch=19)

#all proteoglycans
counts_PGs <- PG_counts[order(PG_counts$cluster_count, decreasing = F),]
countsPGs <- counts_PGs[counts_PGs$cluster_count>0,]
tPGs <- table(countsPGs$cluster_count)
xPGs <- as.numeric(names(tPGs))
yPGs <- as.numeric(tPGs)
plot(xPGs,(cumsum(yPGs)/max(cumsum(yPGs))), cex=0.5, type="o", pch=19)
#all ribosomal genes 
counts_ribos <- ribo_counts[order(ribo_counts$cluster_count, decreasing = F),]
countsribos <- counts_ribos[counts_ribos$cluster_count>0,]
tRibos <- table(countsribos$cluster_count)
xRibos <- as.numeric(names(tRibos))
yRibos <- as.numeric(tRibos)
plot(xRibos,(cumsum(yRibos)/max(cumsum(yRibos))), cex=0.8, type="o", pch=19)
plot(xTFs,(cumsum(yTFs)/max(cumsum(yTFs))), cex=0.8, type="o", pch=19, 
     xlab = "Number of neuron types expressing",
     ylab = "Cumulative fraction",xlim=c(0,60),ylim=c(0,1), col="black")
lines(xC2F,(cumsum(yC2F)/max(cumsum(yC2F))), cex=0.8, type="o", pch=19, col="magenta")
lines(xHD,(cumsum(yHD)/max(cumsum(yHD))), cex=0.8, type="o", pch=19, col="darkgreen")
lines(xGPCRs,(cumsum(yGPCRs)/max(cumsum(yGPCRs))), cex=0.8, type="o", pch=19, col="blue")
lines(xICs,(cumsum(yICs)/max(cumsum(yICs))), cex=0.8, type="o", pch=19, col="red")
lines(xCAMs,(cumsum(yCAMs)/max(cumsum(yCAMs))), cex=0.8, type="o", pch=19, col="orange")
lines(xPGs,(cumsum(yPGs)/max(cumsum(yPGs))), cex=0.8, type="o", pch=19, col="darkblue")
lines(xRibos,(cumsum(yRibos)/max(cumsum(yRibos))), cex=0.8, type="o", pch=19, col="darkred")

legend("bottomright", legend=c("all TFs", "zinc finger-C2H2 ","homeodomain","GPCRs ","ion channels ","cell adhesion molecules ","proteoglycans ","ribosomal genes "),
       col=c("black", "magenta","darkgreen","blue","red","orange","darkblue","darkred"), lty=1, cex=0.35)


dev.off()
#####################################################################################
#####################################################################################
########################################################################################
#######++++++++++++++++++HIEARCHICAL CLUSTERING BY TFS+++++++++++++++++++++++++++++++++
#####################################################################################
#####################################################################################
#FOR MOUSE

library(pvclust)
library(parallel)
library(dendextend)
library(future)
plan("multicore", workers = 4)
options(future.globals.maxSize= 62914560000) # 60GB (60000*1024^2)
library(Seurat)
library(dplyr)
library(cowplot)
library(zoo)
library(ggplot2)
library(scater)
library(stringr)
library(gplots)
library(matrixStats)
sce=seurat_mouse
sce <- SCTransform(sce, method = "glmGamPoi")
Idents(sce) <- "cluster_label"

DefaultAssay(sce)<-"SCT"
neurons_avg <- AverageExpression(sce)
neurons_avg_data <- neurons_avg[["SCT"]]
NoisyGenes <- function(object, min.pct=0.2, clusters.use) {
  clusters <- clusters.use
  genes.use = rownames(object)
  object.raw.data <- as.matrix(GetAssayData(object, slot="counts"))
  pct.matrix = matrix(data=NA, nrow=length(genes.use), ncol=length(clusters))
  rownames(pct.matrix) <- genes.use
  colnames(pct.matrix) <- clusters
  thresh.min=0
  for (i in clusters){
    cells.cluster <- WhichCells(object=object, idents=i)
    data.cluster <- object.raw.data[,colnames(object.raw.data) %in% cells.cluster]
    pct.cluster <- round(apply(object.raw.data[genes.use, cells.cluster, drop = F],1,function(x)return(length(x[x>thresh.min])/length(x))),3)
    pct.matrix[,i] <- pct.cluster
  }
  pct.max <- rowMaxs(pct.matrix)
  names(pct.max) <- genes.use
  noisy.genes <- names(pct.max[pct.max < min.pct])
  return(noisy.genes)
}
Human_TFs <- read.csv("human_tf_gene_list.csv",row.names = 1)
#Human_TFs <- read.csv("human_tf_gene_list2.csv")

noisy.liz <- NoisyGenes(sce, 0.2, levels(Idents(sce)))#0.65
neurons_avg_data_filtered <- neurons_avg_data[!rownames(neurons_avg_data) %in% noisy.liz,]
neurons_avg_TFs <- neurons_avg_data_filtered[rownames(neurons_avg_data_filtered) %in% Human_TFs$HGNC.symbol,]
tfs <- as.data.frame(neurons_avg_TFs)
write.csv(tfs,file="mouse_tfs_for_hierarchical_clustering.csv")

############################################ 244 TFs remain##################################################
pheatmap::pheatmap(neurons_avg_TFs,scale = "row")
dist_func <- function(x){x <- as.dist(1-cor(x, method="s"))}
pvclust.TFs_dist <- pvclust(neurons_avg_TFs, method.dist=dist_func, 
                            method.hclust="ward.D2", nboot=10000, parallel=T)
plot(pvclust.TFs_dist)
dev.off()
pvrect(pvclust.TFs_dist, alpha=0.9)
dend <- as.dendrogram(pvclust.TFs_dist)
dend %>% pvclust_show_signif_gradient(pvclust.TFs_dist, signif_col_fun = colorRampPalette(c("lightgrey","grey","black")), signif_type="au") %>%
  plot(main = "neuronal clusters branch confidence")
dend <- rotate(dend)
dendrogram_names <- labels(dend)
dend <- as.dendrogram(pvclust.TFs_dist)
dend %>%
  pvclust_show_signif(pvclust.TFs_dist, signif_value = c("black", "black"), show_type = "col") %>%
  plot(main = "Cluster dendrogram with AU/BP values (%)")
pvrect2(pvclust.TFs_dist, alpha=0.90)
dendrogram_names <- labels(dend)
levels(sce) <- dendrogram_names
genes <- c("IKZF2","ZNF516","ISL1","ZMAT4","SIX2",#OHC
           "HDAC9","ZNF410","BPNT1","TBX2","BTG2",#IHC
           "NME2","EIF3K","ID3","SOX4","RPL7",#DCOPC_2
           "DEPTOR","ATOH8","IRX6","BCL6","ZNF37A",#DCOPC_3
           "MXD1","PER1","HIF1A","JUN","HES1",#IPHIBC_1
           "CTNNBIP1","ID1","EGR1","HEY2","HEY1",#IPHIBC_2
           "SMAD9","POLK","BACH1","TSC22D1","SMAD1",#DCOPC_1
           "TSC22D2","SOX6","MEIS2","ZNF609","DACH1",#HEC
           "PLAGL1","POT1","EDA","TNFSF4","CUX1",#IPHIBC_5
           "ARHGAP35","PRDM11","RPS6KA5","PHTF1","CUL5",#IPHIBC_3
           "TRPS1","ETV5","PRKD2","RORB","ZNF185",#IPHIBC_4
           "ARNT2","ZNF326","RAPGEF4","IKBKB","LGR4",#PC_1
           "EPAS1","ZBTB16","SMO","TAF4B","KLF13",#PC_3
           "SOX9","AEBP1","BHLHE41","POU3F3","POU3F4",#PC_2
           "RPL7A","RPS27A","PRKCH","JAK2","MEF2A",#TBC_1
           "IRX5","IRX3","CREB5","EPHA5","MET"#TBC_2
)
tfs$gene <- rownames(tfs)
DefaultAssay(sce) <- "SCT"
DotPlot(sce, features = genes,cols = c("white","#284aaa"))+theme(axis.text.x = element_text(angle = 90, hjust = 1))
DotPlot(sce, features = tfs$gene,cols = c("white","#284aaa")
        )+theme(axis.text.x = element_text(angle = 90, hjust = 1))
pdf(file="1067_TFS_mouse.pdf", width=250, height=15)
DotPlot(sce, features = tfs$gene,
        dot.scale = 8,cols  =c("white","#284aaa")
) + RotatedAxis()
dev.off(
  
)
genes <- c("NRG1","ISL1","RUNX1","TLE4","POU3F3","DDR2","PLPP3","CREB5","NR2F2",
           "IRX5","TBX18","POU3F4","PLEKHA4","GULP1","AFF3","CAMK1D","ZFHX4","EBF1",
           "WNT5A","ZEB2","SPI1","MBNL1","POU2F2","SATB2","PTGIS","PAX3","CERS6","TRAK2",
           "LMX1A","SMYD3","MEIS2","TSHZ2","PAX2","MAML3","GATA3","FUS","NME2",
           "NR2F1","NCOA2","PBX1","ATF6","PBX3","TOX","CHD7","MLLT3","ARHGAP35",
           "IKZF2","SALL1","ZBTB7A","TBX2","JUP","ZNF385A","POU4F3"
)
pdf(file="SHARED_TFS_COE.pdf", width=24, height=6)
DotPlot(sce, features = genes,
        dot.scale = 8,cols  =c("white","#284aaa")
)+ RotatedAxis()
dev.off(
  
)
##############################################################################
##############################################################################
##############################################################################
#FOR human
library(pvclust)
library(parallel)
library(dendextend)
library(future)
plan("multicore", workers = 4)
options(future.globals.maxSize= 62914560000) # 60GB (60000*1024^2)
library(Seurat)
library(dplyr)
library(cowplot)
library(zoo)
library(ggplot2)
library(scater)
library(stringr)
library(gplots)
library(matrixStats)
sce=seurat_human
sce <- SCTransform(sce, method = "glmGamPoi")
Idents(sce) <- "cluster_label"

DefaultAssay(sce)<-"SCT"
neurons_avg <- AverageExpression(sce)
neurons_avg_data <- neurons_avg[["SCT"]]
NoisyGenes <- function(object, min.pct=0.2, clusters.use) {
  clusters <- clusters.use
  genes.use = rownames(object)
  object.raw.data <- as.matrix(GetAssayData(object, slot="counts"))
  pct.matrix = matrix(data=NA, nrow=length(genes.use), ncol=length(clusters))
  rownames(pct.matrix) <- genes.use
  colnames(pct.matrix) <- clusters
  thresh.min=0
  for (i in clusters){
    cells.cluster <- WhichCells(object=object, idents=i)
    data.cluster <- object.raw.data[,colnames(object.raw.data) %in% cells.cluster]
    pct.cluster <- round(apply(object.raw.data[genes.use, cells.cluster, drop = F],1,function(x)return(length(x[x>thresh.min])/length(x))),3)
    pct.matrix[,i] <- pct.cluster
  }
  pct.max <- rowMaxs(pct.matrix)
  names(pct.max) <- genes.use
  noisy.genes <- names(pct.max[pct.max < min.pct])
  return(noisy.genes)
}
Human_TFs <- read.csv("human_tf_gene_list.csv",row.names = 1)
#Human_TFs <- read.csv("human_tf_gene_list2.csv")

noisy.liz <- NoisyGenes(sce, 0.2, levels(Idents(sce)))#0.75
neurons_avg_data_filtered <- neurons_avg_data[!rownames(neurons_avg_data) %in% noisy.liz,]
neurons_avg_TFs <- neurons_avg_data_filtered[rownames(neurons_avg_data_filtered) %in% Human_TFs$HGNC.symbol,]
tfs <- as.data.frame(neurons_avg_TFs)
write.csv(tfs,file="human_tfs_for_hierarchical_clustering.csv")
###################################################249 TFs remain###############################################
pheatmap::pheatmap(neurons_avg_TFs,scale = "row")
dist_func <- function(x){x <- as.dist(1-cor(x, method="s"))}
pvclust.TFs_dist <- pvclust(neurons_avg_TFs, method.dist=dist_func, 
                            method.hclust="ward.D2", nboot=10000, parallel=T)
plot(pvclust.TFs_dist)
dev.off()
pvrect(pvclust.TFs_dist, alpha=0.9)
dend <- as.dendrogram(pvclust.TFs_dist)
dend %>% pvclust_show_signif_gradient(pvclust.TFs_dist, signif_col_fun = colorRampPalette(c("lightgrey","grey","black")), signif_type="au") %>%
  plot(main = "neuronal clusters branch confidence")
dend <- rotate(dend)
dendrogram_names <- labels(dend)
dend <- as.dendrogram(pvclust.TFs_dist)
dend %>%
  pvclust_show_signif(pvclust.TFs_dist, signif_value = c("black", "black"), show_type = "col") %>%
  plot(main = "Cluster dendrogram with AU/BP values (%)")
pvrect2(pvclust.TFs_dist, alpha=0.90)
dendrogram_names <- labels(dend)
levels(sce) <- dendrogram_names
genes <- c("IKZF2","ZNF516","ISL1","ZMAT4","SIX2",#OHC
           "HDAC9","ZNF410","BPNT1","TBX2","BTG2",#IHC
           "NME2","EIF3K","ID3","SOX4","RPL7",#DCOPC_2
           "DEPTOR","ATOH8","IRX6","BCL6","ZNF37A",#DCOPC_3
           "MXD1","PER1","HIF1A","JUN","HES1",#IPHIBC_1
           "CTNNBIP1","ID1","EGR1","HEY2","HEY1",#IPHIBC_2
           "SMAD9","POLK","BACH1","TSC22D1","SMAD1",#DCOPC_1
           "TSC22D2","SOX6","MEIS2","ZNF609","DACH1",#HEC
           "PLAGL1","POT1","EDA","TNFSF4","CUX1",#IPHIBC_5
           "ARHGAP35","PRDM11","RPS6KA5","PHTF1","CUL5",#IPHIBC_3
           "TRPS1","ETV5","PRKD2","RORB","ZNF185",#IPHIBC_4
           "ARNT2","ZNF326","RAPGEF4","IKBKB","LGR4",#PC_1
           "EPAS1","ZBTB16","SMO","TAF4B","KLF13",#PC_3
           "SOX9","AEBP1","BHLHE41","POU3F3","POU3F4",#PC_2
           "RPL7A","RPS27A","PRKCH","JAK2","MEF2A",#TBC_1
           "IRX5","IRX3","CREB5","EPHA5","MET"#TBC_2
)
tfs$gene <- rownames(tfs)
DefaultAssay(sce) <- "SCT"
DotPlot(sce, features = genes,cols = c("white","#284aaa"))+theme(axis.text.x = element_text(angle = 90, hjust = 1))
DotPlot(sce, features = tfs$gene,cols = c("white","#284aaa")
)+theme(axis.text.x = element_text(angle = 90, hjust = 1))
pdf(file="988_TFS_human.pdf", width=250, height=15)
DotPlot(sce, features = tfs$gene,
        dot.scale = 8,cols  =c("white","#284aaa")
) + RotatedAxis()
dev.off(
  
)
genes <- c("EGLN1","EPAS1","RAPGEF4","CFLAR","MECOM","HIVEP3","PAX3","MITF",
           "MET","ZEB2","KAT2B","RELN","GAS7","ZNF536","PLEKHA4","REL","PLSCR1",
           "TLR2","HMGB2","ZNF451","EZH2","TOP2A","PAWR","ZFHX3","NOTCH3","POU3F1",
           "MAF","MKRN1","NCOA4","SATB2","MEF2C","RUNX2","ZNF326","WNT5A","LGR4",
           "PLPP3","DDR2","TBX18","ZNF385D","CREB5","BNC2","ZIC1","CAMTA1","FOS",
           "MYT1L","ZNF804A","EPHA5","TLE4","PBX3","TSHZ3","PBX1","ZHX2","GATA3",
           "SOX6","ADAMTS17","ZHX3","BMP7","HDAC9","NRG1","ZMAT4","RORB","ZFP64",
           "IKZF2","LMX1A","ALK","MEIS1","SGK3","CUX2","DACH1","ESRRG"

           
)

DotPlot(sce, features = genes,
        dot.scale = 8,cols  =c("white","#284aaa")
)+ RotatedAxis()
dev.off(
  
)


##########################################################################################
###########################################################################################
########===========================SCTNORMALIZATION & INTEGRATION===========================
###############################################################################################
################################################################################################
seurat_human <- readRDS("C:/Users/Dell/Desktop/WORKPLACE/R/seurat_human_new.rds")
seurat_mouse <- readRDS("C:/Users/Dell/Desktop/WORKPLACE/R/seurat_mouse_new.rds")
######EXCLUDE SUBCLASSES THAT DUE TO TECHNIQUE REASONS WHICH COULD NOT BE CAPTURED.
###CEC
DefaultAssay(seurat_human) <- "RNA"
DefaultAssay(seurat_mouse) <- "RNA"
seurat_human$species <- "human"
seurat_mouse$species <- "mouse"
table(seurat_human$orig.ident,seurat_human$species)
table(seurat_mouse$orig.ident,seurat_mouse$species)
table(seurat_human$subclasses_label)
table(seurat_mouse$subclasses_label)

#for human species, we exclude cec,ec,pvm/m,eb,cc
Idents(seurat_human) <- "subclasses_label"
seurat_human <- subset(seurat_human,idents=c("CEC","ECs","PVM/M","EBs","CCs"), invert=TRUE)
Idents(seurat_human) <- "cluster_label"
table(seurat_human$cluster_label)
seurat_human <- subset(seurat_human,idents=c("CEC_1","CEC_2","CEC_4","CEC_5",
                                             "CEC_6","CEC_7","ECs","EBs","CCs"), 
                       invert=TRUE)
seurat_human <- subset(seurat_human,idents=c("PVM/M"), 
                       invert=TRUE)
table(seurat_human$cluster_label)
p1=DimPlot(seurat_human,label = TRUE)
p1

#for mouse species, we exclude PC
Idents(seurat_mouse) <- "subclasses_label"
seurat_mouse <- subset(seurat_mouse,idents=c("PC"), invert=TRUE)
Idents(seurat_mouse) <- "cluster_label"
table(seurat_mouse$cluster_label)
seurat_mouse <- subset(seurat_mouse,idents=c("PC_1"), invert=TRUE)
table(seurat_mouse$cluster_label)
p2=DimPlot(seurat_mouse,label = TRUE)
p1+p2
DefaultAssay(seurat_mouse) <- "SCT"
FeaturePlot(seurat_mouse,features = c("COCH","SLC4A10"),label = TRUE)

#RENAME SUBCLASSES_LABEL AND CLUSTER_LABEL
Idents(seurat_mouse) <- "subclasses_label"
table(Idents(seurat_mouse))
new.cluster.ids <- c("TBC","BC","DC_PC_HeC","MC","Fibro3","Fibro1","HC",
                     "Fibro2","IPh_IBC","SGN","IC","Immun","Fibro4",
                     "Fib","GC","SMC")
names(new.cluster.ids) <- levels(seurat_mouse)
seurat_mouse<- RenameIdents(seurat_mouse, new.cluster.ids)
DimPlot(seurat_mouse, reduction = "umap",label = TRUE,repel = TRUE)
seurat_mouse$subclasses_label.new <- Idents(seurat_mouse)

Idents(seurat_mouse) <- "cluster_label"
table(Idents(seurat_mouse))
new.cluster.ids <- c("BC_1","TBC_1","DC_PC_HeC_1","TBC_2","Fibro3","MC_1","Fibro1",
                     "OHC","Fibro2","IPh_IBC_1","IC","SGN_1","Immun",
                     "BC_2","SGN_2","DC_PC_HeC_2","MC_2","DC_PC_HeC_3",
                     "Fibro4_1","DC_PC_HeC_4","IPh_IBC_2","IHC","Fib_1","Fib_2",
                     "GC","Fibro4_2","IPh_IBC_3","SGN_3","SMC","SGN_4")
names(new.cluster.ids) <- levels(seurat_mouse)
seurat_mouse<- RenameIdents(seurat_mouse, new.cluster.ids)
DimPlot(seurat_mouse, reduction = "umap",label = TRUE,repel = TRUE)
seurat_mouse$cluster_label.new <- Idents(seurat_mouse)

seurat_human$subclasses_label.new <- seurat_human$subclasses_label
table(seurat_human$cluster_label,seurat_human$seurat_clusters)

a <- table(seurat_human$cluster_label,seurat_human$seurat_clusters)
a <- as.data.frame.array(a)
write.csv(a,file = "cell_proportion_of_clusetrs_human.csv")
p1=DimPlot(seurat_human,group.by = "subclasses_label",label = TRUE)
p2=DimPlot(seurat_human,group.by = "seurat_clusters",label = TRUE)
p1+p2
Idents(seurat_human) <- "seurat_clusters"
table(Idents(seurat_human))
new.cluster.ids <- c("SMC_1","TBC_1","BC_1","GC_1","BC_2","GC_2","GC_3",
                     "IC_1","Fibro1_1","Fibro1_2","Fibro1_3","SMC_2","GC_4",
                     "IC_2","SMC_3","GC_5","Fibro2_1","TBC_2","SMC_4","Fibro1_4",
                     "Fibro1_5","BC_3","Fibro2_2","Fibro1_6","DC_PC_HeC_1","BC_4",
                     "Immun","IPh_IBC","Fibro2_3","Fibro4","SGN_1","Fibro1_7",
                     "Fibro3","Fibro1_8","SMC_5","IC_3","MC","Fib","DC_PC_HeC_2",
                     "DC_PC_HeC_3","BC_5","SGN_2","SGN_3","HC")
names(new.cluster.ids) <- levels(seurat_human)
seurat_human<- RenameIdents(seurat_human, new.cluster.ids)
DimPlot(seurat_human, reduction = "umap",label = TRUE,repel = TRUE)
seurat_human$cluster_label.new <- Idents(seurat_human)

Idents(seurat_human) <- "seurat_clusters"
table(Idents(seurat_human))
new.cluster.ids <- c("SMC","TBC","BC","GC","BC","GC","GC",
                     "IC","Fibro1","Fibro1","Fibro1","SMC","GC",
                     "IC","SMC","GC","Fibro2","TBC","SMC","Fibro1",
                     "Fibro1","BC","Fibro2","Fibro1","DC_PC_HeC","BC",
                     "Immun","IPh_IBC","Fibro2","Fibro4","SGN","Fibro1",
                     "Fibro3","Fibro1","SMC","IC","MC","Fib","DC_PC_HeC",
                     "DC_PC_HeC","BC","SGN","SGN","HC")
names(new.cluster.ids) <- levels(seurat_human)
seurat_human<- RenameIdents(seurat_human, new.cluster.ids)
DimPlot(seurat_human, reduction = "umap",label = TRUE,repel = TRUE)
seurat_human$subclasses_label.new <- Idents(seurat_human)

#######################################################
#weighted average downsample... 200 cells max per pre-integrated cluster for species with most subclasses
human_subclass <- table(seurat_human$subclasses_label.new,
                        seurat_human$cluster_label.new)
mouse_subclass <- table(seurat_mouse$subclasses_label.new, 
                        seurat_mouse$cluster_label.new)

all_subclasses <- data.frame(subclass = unique(c(rownames(human_subclass), 
                                                 rownames(mouse_subclass))))
all_subclasses$human <- 0
all_subclasses$mouse <- 0

for(i in 1:nrow(all_subclasses)){ #calculate how many clusters per subclass
  all_subclasses$human[i] <- length(which(human_subclass[which(rownames(human_subclass) == all_subclasses$subclass[i]), ] > 0))
  
  all_subclasses$mouse[i] <- length(which(mouse_subclass[which(rownames(mouse_subclass) == all_subclasses$subclass[i]), ] > 0))
}
all_subclasses

all_subclasses$max_cells <- 0
for(i in 1:nrow(all_subclasses)){ #calculate max cells per subclass
  all_subclasses$max_cells[i] <- all_subclasses[i, 1 + which.max(all_subclasses[i, 2:3])] * 200
}
all_subclasses
all_subclasses$human_cells_per_cl <- 0
all_subclasses$mouse_cells_per_cl <- 0

for(i in 1:nrow(all_subclasses)){
  all_subclasses$human_cells_per_cl[i] <- round(all_subclasses$max_cells[i] / all_subclasses$human[i])
  
  all_subclasses$mouse_cells_per_cl[i] <- round(all_subclasses$max_cells[i] / all_subclasses$mouse[i])
}
all_subclasses

#subset human data
cells_to_keep <- NA
for(i in 1:nrow(all_subclasses)){
  clusters_to_subset <- unique(seurat_human$cluster_label.new[which(seurat_human$subclasses_label.new %in% all_subclasses$subclass[i])])
  for(p in 1:length(clusters_to_subset)){
    if(length(which(seurat_human$cluster_label.new == clusters_to_subset[p])) > all_subclasses$human_cells_per_cl[i]){
      cells_to_keep <- c(cells_to_keep, sample(which(seurat_human$cluster_label.new == clusters_to_subset[p]), all_subclasses$human_cells_per_cl[i]))
    }
    if(length(which(seurat_human$cluster_label.new == clusters_to_subset[p])) <= all_subclasses$human_cells_per_cl[i]){
      cells_to_keep <- c(cells_to_keep, sample(which(seurat_human$cluster_label.new == clusters_to_subset[p])))
    }
  }
}
table(seurat_human$cluster_label.new)
cells_to_keep <- cells_to_keep[-1]
seurat_human
seurat_human$sample_id <- seq(1,30249)
Idents(seurat_human) <- seurat_human$sample_id
seurat_human <- subset(seurat_human, idents = cells_to_keep)
table(seurat_human$cluster_label.new)
table(seurat_human$seurat_clusters)

#subset mouse data
cells_to_keep <- NA
for(i in 1:nrow(all_subclasses)){
  clusters_to_subset <- unique(seurat_mouse$cluster_label.new[which(seurat_mouse$subclasses_label.new %in% all_subclasses$subclass[i])])
  for(p in 1:length(clusters_to_subset)){
    if(length(which(seurat_mouse$cluster_label.new == clusters_to_subset[p])) > all_subclasses$mouse_cells_per_cl[i]){
      cells_to_keep <- c(cells_to_keep, sample(which(seurat_mouse$cluster_label.new == clusters_to_subset[p]), all_subclasses$mouse_cells_per_cl[i]))
    }
    if(length(which(seurat_mouse$cluster_label.new == clusters_to_subset[p])) <= all_subclasses$mouse_cells_per_cl[i]){
      cells_to_keep <- c(cells_to_keep, sample(which(seurat_mouse$cluster_label.new == clusters_to_subset[p])))
    }
  }
}
table(seurat_mouse$cluster_label.new)
cells_to_keep <- cells_to_keep[-1]
seurat_mouse
seurat_mouse$sample_id <- seq(1,14855)
Idents(seurat_mouse) <- seurat_mouse$sample_id
seurat_mouse <- subset(seurat_mouse, idents = cells_to_keep)
table(seurat_mouse$cluster_label.new)
table(seurat_mouse$seurat_clusters)

###############################SCT norm and integration  ###############################
all.data <- merge(x = seurat_human, y = seurat_mouse)

combined.list <- SplitObject(all.data, split.by = "species")
rm(all.data)
for (i in 1:length(combined.list)) {
  combined.list[[i]] <- SCTransform(combined.list[[i]], vars.to.regress = "nCount_RNA",
                                    verbose = TRUE, method="glmGamPoi")
}
saveRDS(combined.list,file = "combined_list_with_weight_average_downsample.rds")

################################################################################
################################################################################
########===============WITHOUT WEIGHT AVERAGE DOWNSAMPLE=======================
################################################################################
###############################################################################
all.data <- merge(x = seurat_human, y = seurat_mouse)

combined.list <- SplitObject(all.data, split.by = "species")
rm(all.data)
for (i in 1:length(combined.list)) {
  combined.list[[i]] <- SCTransform(combined.list[[i]], vars.to.regress = "nCount_RNA",
                                    verbose = TRUE, method="glmGamPoi")
}
saveRDS(combined.list,file = "combined_list.rds")
############################################################################################
library(scrattch.hicat)
table(seurat_human$cluster_label.new)
table(seurat_human$seurat_clusters)
Var.genes.human         <- select_markers(combined.list$human@assays$SCT@counts,
                                          combined.list$human$seurat_clusters, n.markers = 100)
Var.genes.human.markers <- Var.genes.human$markers
table(seurat_mouse$cluster_label.new)
table(seurat_mouse$seurat_clusters)
Var.genes.mouse         <- select_markers(combined.list$mouse@assays$SCT@counts, 
                                          combined.list$mouse$seurat_clusters, n.markers = 100)
Var.genes.mouse.markers <- Var.genes.mouse$markers

total.Var.genes <- unique(c(Var.genes.human$markers,Var.genes.mouse$markers))
total.Var.genes <- total.Var.genes[which(total.Var.genes %in% rownames(combined.list$human@assays$SCT@counts))]
total.Var.genes <- total.Var.genes[which(total.Var.genes %in% rownames(combined.list$mouse@assays$SCT@counts))]
write.csv(total.Var.genes,file = "features used for cca anchors with weight average downsample.csv")
features=total.Var.genes
##########
#features <- SelectIntegrationFeatures(object.list =combined.list, nfeatures = 3000)
#library(future)
#options(future.globals.maxSize = 10000 * 1024^2)
combined.list <- PrepSCTIntegration(object.list =combined.list,
                                    anchor.features = features)
cochlea.anchors <- FindIntegrationAnchors(object.list = combined.list,
                                          normalization.method = "SCT", 
                                          anchor.features = features,dims = 1:40,
                                          reduction="cca")
#rm(combined.list)
cochlea.combined.sct <- IntegrateData(anchorset = cochlea.anchors, 
                                      normalization.method = "SCT",dims = 1:40)
DefaultAssay(cochlea.combined.sct) <- "integrated"
cochlea.combined.sct <- RunPCA(cochlea.combined.sct, features = features, npcs=200)
ElbowPlot(cochlea.combined.sct , ndims = 200)
cochlea.combined.sct <- FindNeighbors(cochlea.combined.sct, 
                                      reduction = "pca", dims = 1:50,nn.eps = 0)
cochlea.combined.sct <- FindClusters(cochlea.combined.sct)
saveRDS(cochlea.anchors,file="cochlea_anchors_integrated_with_weight_average_downsample.rds")

cochlea.combined.sct <- RunUMAP(cochlea.combined.sct, reduction = "pca", dims = 1:50)
#cochlea.combined.sct <- RunTSNE(cochlea.combined.sct, reduction = "pca", dims = 1:30)
cochlea.combined.sct$seurat_clusters.new <- as.integer(cochlea.combined.sct$seurat_clusters)
# Visualization
library(patchwork)
p1 <- DimPlot(cochlea.combined.sct, reduction = "umap", group.by = "species")
p2 <- DimPlot(cochlea.combined.sct, reduction = "umap", label = TRUE)
p3 <- DimPlot(cochlea.combined.sct, reduction = "umap", label = TRUE,group.by = "subclasses_label.new")
p1+p2+p3
p4 <- DimPlot(cochlea.combined.sct, reduction = "umap", label = TRUE,
              group.by = "subclasses_label",split.by = "species")
p4
DefaultAssay(cochlea.combined.sct) <- "SCT"
VlnPlot(cochlea.combined.sct,features =c("OTOF","SLC26A5","TMC1"),
        split.by = "species",stack = TRUE,flip = TRUE )
FeaturePlot(cochlea.combined.sct,features =c("OTOF","SLC26A5","TMC1"),
        label = TRUE,split.by = "species")
FeaturePlot(cochlea.combined.sct,features =c("SLC26A5"),
            label = TRUE,split.by = "species")
dev.off()
Idents(cochlea.combined.sct) <- "seurat_clusters"
cochlea.combined.sct2 <- subset(cochlea.combined.sct,idents="28",invert=TRUE)


DefaultAssay(cochlea.combined.sct2) <- "integrated"
cochlea.combined.sct2 <- RunPCA(cochlea.combined.sct2, features = features, npcs=100)
ElbowPlot(cochlea.combined.sct2 , ndims = 100)
cochlea.combined.sct2 <- FindNeighbors(cochlea.combined.sct2, 
                                      reduction = "pca", dims = 1:50,nn.eps = 0)
cochlea.combined.sct2 <- FindClusters(cochlea.combined.sct2)


cochlea.combined.sct2 <- RunUMAP(cochlea.combined.sct2, reduction = "pca", dims = 1:50)
#cochlea.combined.sct <- RunTSNE(cochlea.combined.sct, reduction = "pca", dims = 1:30)
cochlea.combined.sct2$seurat_clusters.new <- as.integer(cochlea.combined.sct2$seurat_clusters)
# Visualization
library(patchwork)
p1 <- DimPlot(cochlea.combined.sct2, reduction = "umap", group.by = "species")
p2 <- DimPlot(cochlea.combined.sct2, reduction = "umap", label = TRUE)
p3 <- DimPlot(cochlea.combined.sct2, reduction = "umap", label = TRUE,group.by = "subclasses_label.new")
p1+p2+p3
DefaultAssay(cochlea.combined.sct2) <- "SCT"
VlnPlot(cochlea.combined.sct2,features =c("OTOF","SLC26A5","TMC1"),
        split.by = "species",stack = TRUE,flip = TRUE )
FeaturePlot(cochlea.combined.sct2,features =c("OTOF","SLC26A5","TMC1"),
            split.by = "species",pt.size = 1)
FeaturePlot(cochlea.combined.sct2,features =c("SLC26A5"),
            label = TRUE,split.by = "species")
dev.off()

table(cochlea.combined.sct2$subclasses_label.new,cochlea.combined.sct2$seurat_clusters.new)
#####============================name integrated_subclasses or _clusters
DimPlot(cochlea.combined.sct2, reduction = "umap",label = TRUE)
DimPlot(cochlea.combined.sct2, reduction = "umap",label = TRUE,group.by = "subclasses_label.new")

p1 <- DimPlot(cochlea.combined.sct2, reduction = "umap", split.by = "species",label = TRUE,group.by = "subclasses_label.new")
p2 <- DimPlot(cochlea.combined.sct2, reduction = "umap", label = TRUE,group.by = "seurat_clusters")
p1+p2
p3 <- DimPlot(cochlea.combined.sct2, reduction = "umap",label = TRUE,group.by = "subclasses_label.new")
p4 <- DimPlot(cochlea.combined.sct2, reduction = "umap", label = TRUE,group.by = "seurat_clusters")
p3+p4
p1+p2+p3
p4
p5 <- DimPlot(cochlea.combined.sct2, reduction = "umap", split.by = "species",label = TRUE,group.by = "cluster_label.new")
p6 <- DimPlot(cochlea.combined.sct2, reduction = "umap", label = TRUE,group.by = "seurat_clusters",split.by = "species")
p5+p6
dev.off()
cochlea.combined.sct2$subclasses.species <- paste(cochlea.combined.sct2$subclasses_label.new,cochlea.combined.sct2$species,
                                                 sep = "_")
a <- table(cochlea.combined.sct2$subclasses.species,cochlea.combined.sct2$seurat_clusters.new)
a <- as.data.frame.array(a)
write.csv(a,file = "cell_proportion_of_integrated_cluster_new.csv")
table(seurat_human$cluster_label.new)
table(seurat_mouse$cluster_label.new)
saveRDS(cochlea.combined.sct2,file="Cochlea_integrated.rds")

###rename integrated_subclasses
Idents(cochlea.combined.sct2) <- "seurat_clusters"
table(Idents(cochlea.combined.sct2))
new.cluster.ids <- c("BC","Fibro1","TBC","MC","BC","GC","TBC","SMC",
                     "DC_PC_HeC","IC","Fibro4","DC_PC_HeC","BC","Immun",
                     "DC_PC_HeC","SMC","SGN","SGN","Fibro3","Fibro2","Fibro1","IPh_IBC",
                     "Fib","IPh_IBC","GC","IPh_IBC","Fibro3","mixture","Fibro1","HC",
                     "MC","IC","BC","Fibro3","IPh_IBC","SMC","DC_PC_HeC","SGN","HC",
                     "SGN","SGN","IPh_IBC")
names(new.cluster.ids) <- levels(cochlea.combined.sct2)
cochlea.combined.sct2<- RenameIdents(cochlea.combined.sct2, new.cluster.ids)
DimPlot(cochlea.combined.sct2, reduction = "umap",label = TRUE,repel = TRUE)
cochlea.combined.sct2$integrated_subclasses <- Idents(cochlea.combined.sct2)
###rename integrated_clusters
Idents(cochlea.combined.sct2) <- "seurat_clusters"
table(Idents(cochlea.combined.sct2))
new.cluster.ids <- c("BC_1","Fibro1_1","TBC_1","MC_1","BC_2","GC_1","TBC_2","SMC_1",
                     "DC_PC_HeC_1","IC_1","Fibro4","DC_PC_HeC_2","BC_3","Immun",
                     "DC_PC_HeC_3","SMC_2","SGN_1","SGN_2","Fibro3_1","Fibro2",
                     "Fibro1_2","IPh_IBC_1",
                     "Fib","IPh_IBC_2","GC_2","IPh_IBC_3","Fibro3_2","mixture","Fibro1_3","HC_1",
                     "MC_2","IC_2","BC_4","Fibro3_3","IPh_IBC_4","SMC_3","DC_PC_HeC_4",
                     "SGN_3","HC_2",
                     "SGN_4","SGN_5","IPh_IBC_5")
names(new.cluster.ids) <- levels(cochlea.combined.sct2)
cochlea.combined.sct2<- RenameIdents(cochlea.combined.sct2, new.cluster.ids)
DimPlot(cochlea.combined.sct2, reduction = "umap",label = TRUE,repel = TRUE)
DimPlot(cochlea.combined.sct2, reduction = "umap",label = TRUE,repel = TRUE,
        group.by = "integrated_subclasses")
cochlea.combined.sct2$integrated_clusters <- Idents(cochlea.combined.sct2)
## REMOVE MIXTURE
Idents(cochlea.combined.sct2) <- "integrated_subclasses"
cochlea.combined.sct2 <- subset(cochlea.combined.sct2,idents="mixture",invert=TRUE)
###rename classes
Idents(cochlea.combined.sct2) <- "integrated_subclasses"
table(Idents(cochlea.combined.sct2))
#Stria Vascularis; Modiolus;Spiral Ligament; Cochlear Epithelium;
new.cluster.ids <- c("Stria Vascularis","Spiral Ligament","Cochlear Epithelium",
                     "Stria Vascularis","Modiolus","Spiral Ligament","Cochlear Epithelium",
                     "Stria Vascularis","Spiral Ligament","Immun","Modiolus","Spiral Ligament",
                     "Spiral Ligament","Cochlear Epithelium","Spiral Ligament","Cochlear Epithelium"
  )
names(new.cluster.ids) <- levels(cochlea.combined.sct2)
cochlea.combined.sct2<- RenameIdents(cochlea.combined.sct2, new.cluster.ids)
DimPlot(cochlea.combined.sct2, reduction = "umap",label = TRUE,repel = TRUE)
cochlea.combined.sct2$integrated_classes <- Idents(cochlea.combined.sct2)

dev.off()
saveRDS(cochlea.combined.sct2,file="Cochlea_integrated.rds")
Cochlea_integrated <- readRDS("C:/Users/Dell/Desktop/WORKPLACE/R/Cochlea_integrated.rds")
cochlea.combined.sct=Cochlea_integrated
a <- table(cochlea.combined.sct2$subclasses.species,cochlea.combined.sct2$integrated_clusters)
a <- as.data.frame.array(a)
write.csv(a,file = "cell_proportion_of_integrated_cluster_rename.csv")
###change colors
# for species human: #2EC794   mouse: #612EC7
data_to_plot <- data.frame(cochlea.combined.sct2@reductions$umap@cell.embeddings)
cochlea.combined.sct2$new_id <- colnames(cochlea.combined.sct2)
data_to_plot$species <- cochlea.combined.sct2$species[match(rownames(data_to_plot), cochlea.combined.sct2$new_id)]
data_to_plot$plot_order <- 0
for(i in 1:nrow(data_to_plot)){
  data_to_plot$plot_order[i] <- sample(1:100000, 1)
  print(i)
}
data_to_plot <- arrange(data_to_plot, plot_order)

ggplot(data_to_plot,
       aes(x = UMAP_1, y = UMAP_2, color = species)) +
  geom_point(size = 0.5)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_color_manual(values = c("#308ce1", "#f28eac"))
dev.off()

# for COCHLEA DIVISION Stria Vascularis #7B68EE; Modiolus #3CB371;Spiral Ligament #FFDEAD; Cochlear Epithelium #F08080


##C0C0C0
data_to_plot <- data.frame(cochlea.combined.sct2@reductions$umap@cell.embeddings)
cochlea.combined.sct2$new_id <- colnames(cochlea.combined.sct2)
data_to_plot$class <- cochlea.combined.sct2$integrated_classes[match(rownames(data_to_plot), cochlea.combined.sct2$new_id)]
table(data_to_plot$class)
data_to_plot$plot_order <- 0
for(i in 1:nrow(data_to_plot)){
  data_to_plot$plot_order[i] <- sample(1:100000, 1)
  print(i)
}
data_to_plot <- arrange(data_to_plot, plot_order)

ggplot(data_to_plot,
       aes(x = UMAP_1, y = UMAP_2, color = class)) +
  geom_point(size = 0.5)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_color_manual(values = c("#536ac7", "#90e5cb", "#e5e08f", "#f586b7","gray"))


#FOR INTEGRATED_CLUSTERS
# umap cluster colors
table(cochlea.combined.sct2$integrated_clusters)
levels(Idents(cochlea.combined.sct2))
Idents(cochlea.combined.sct2) <- "integrated_clusters"
my_cols <- c("BC_1"="#FFA07A","Fibro1_1"="#FFD700","TBC_1"="#FA8072","MC_1"="#FFFACD",
             "BC_2"="#FFA500","GC_1"="#FF8C00",
             "TBC_2"="#FFEFD5","SMC_1"="#32CD32","DC_PC_HeC_1"="#a30ad9",
             "IC_1"="#EE82EE","Fibro4"="#CD5C5C","DC_PC_HeC_2"="#FF69B4",
             "BC_3"="#AFEEEE","Immun"="#43b3ae","DC_PC_HeC_3"="#FFC0CB",
             "SMC_2"="#F4A460",
             "SGN_1"="#90EE90","SGN_2"="#E0FFFF","Fibro3_1"="#EEE8AA",
             "Fibro2"="#DAA520","Fibro1_2"="#F0E68C",
             "IPh_IBC_1"="#98FB98",
             "Fib"="#2E8B57","IPh_IBC_2"="#DC143C","GC_2"="#ADFF2F","IPh_IBC_3"="#FFF8DC",
             "Fibro3_2"="#7FFFD4","Fibro1_3"="#3CB371",
             "HC_1"="#8B0000","MC_2"="#8FBC8F","IC_2"="#DEB887",
             "BC_4"="#FFD700",
             "Fibro3_3"="#00FF00",
             "SMC_3"="#800080","DC_PC_HeC_4"="#A52A2A","SGN_3"="#1fe9f2",
             "HC_2"="#800000","SGN_4"="#abf2fa","SGN_5"="#36913a","IPh_IBC_5"="#b445a7","IPh_IBC_4"="#230aba"
             
)
new.cluster.ids <- c("TBC_1","BC_1","BC_2","SMC_1","Fibro1/2_1","GC_1","TBC_2","IC_1",
                     "GC_2","SMC_2","MC_1","DC_PC_HeC_1","GC_3","DC_PC_HeC_2","Fibro1/2_2",
                     "SMC_3","Fibro1/2_3","BC_3","Fibro1/2_4","IC_2","Fibro1/2_5",
                     "Fibro1/2_6","Fibro1/2_7","Immun","IPh_IBC_1","Fibro1/2_8","GC_4",
                     "Fibro4","Fibro3_1","SGN_1","SGN_2","GC_5","Fibro3_2","Mixture","HC_1",
                     "Fibro3_3","BC_4","SMC_4","IC_3","Fib","IPh_IBC_2","DC_PC_HeC_3","HC_2",
                     "DC_PC_HeC_4","IPh_IBC_3","SGN_3","SGN_4","HC_3","DC_PC_HeC_5","IPh_IBC_4",
                     "MC_2","SGN_5","SGN_6")


my_cols2 <- my_cols[order(as.integer(names(my_cols)))]
scales::show_col(my_cols2)
DimPlot(cochlea.combined.sct2,
        cols = my_cols2, label=TRUE , repel=TRUE,reduction = "umap",pt.size = 0.5
        ,split.by = "species",group.by = "integrated_clusters")

DimPlot(cochlea.combined.sct2,
        cols = my_cols2, label=TRUE , repel=TRUE,reduction = "umap",pt.size = 0.5
        ,group.by = "integrated_clusters")

p2=DimPlot(cochlea.combined.sct2,
         label=TRUE , repel=TRUE,reduction = "umap",pt.size = 0.5
        ,group.by = "integrated_subclasses")
p2
data_to_plot$species <- cochlea.combined.sct2$species[match(rownames(data_to_plot), cochlea.combined.sct2$new_id)]
data_to_plot$plot_order <- 0
for(i in 1:nrow(data_to_plot)){
  data_to_plot$plot_order[i] <- sample(1:100000, 1)
  print(i)
}
data_to_plot <- arrange(data_to_plot, plot_order)
p3=ggplot(data_to_plot,
       aes(x = UMAP_1, y = UMAP_2, color = species)) +
  geom_point(size = 0.5)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_color_manual(values = c("#308ce1", "#f28eac"))
p3
library(patchwork)
p3+p2
dev.off()
#####################################################################################
#####################################################################################
#####################################################################################
################################################################################################
library(conos)
library(pagoda2)
devtools::install_github("huqiwen0313/speciesTree")
library(speciesTree)
# building the Dendrogram with leaves colores by species' mixing
# get the data from integrated assay
Idents(cochlea.combined.sct2) <- "integrated_clusters"
sce=cochlea.combined.sct2
sce=Cochlea_integrated
Idents(sce)<- "integrated_clusters"
table(sce$integrated_clusters)
# renaming clusters starting from 1 instead of 0
DefaultAssay(sce) <- "integrated"
new.ids <- seq(1, length(levels(Idents(sce))), by=1)
names(new.ids) <- levels(Idents(sce))
sce <- RenameIdents(sce, new.ids)
species_annot <- sce@meta.data$species
meta_df <- sce@meta.data
names(species_annot) <- rownames(sce@meta.data)
meta_df$speciesTree <- species_annot

sce <- AddMetaData(sce, meta_df)
comb.integrated=sce
rm(sce)
expression.matrix <- GetAssayData(comb.integrated, slot="data", assay="integrated")
#Idents(comb.integrated) <- "integrated_clusters"
meta_clusters <- as.integer(Idents(comb.integrated))# used for color-base tree
meta_clusters <- Idents(comb.integrated)# optional
#names(meta_clusters) <- names(Idents(comb.integrated))
upperlevelinfo = NULL
species_annot <- comb.integrated@meta.data$species
names(species_annot) <- rownames(comb.integrated@meta.data)

# obtaining distance matrix and performing hierarchical clustering
d <- cluster.matrix.expression.distances(t(expression.matrix), groups=meta_clusters, dist="cor", 
                                         useVariablegenes=FALSE,  use.scaled.data=TRUE)

dendr <- hclust(as.dist(d), method='ward.D2')
dend <- as.dendrogram(dendr)
# coloring dendrogram leaves according to the proportion of cells from each species in each cluster
dendr <- TransferDend(dend, renameCluster = FALSE, cls.groups = meta_clusters)
cls.groups <- dendr$new.groups
dend <- dendr$dendrogram
leafcontent <- dendr$leafcontent
stability.measurements = NULL
dend <- AddTreeAttribute(dend, species_annot, leafcontent)
dend <- dendSetWidthBysize(dend, scale = 8)
colorpallete <- colorRampPalette(c("blue", "grey", "grey",  "grey", "red"))(101)

scales::show_col(colorpallete)
upperLevelnodes = NULL
fac <- as.factor(species_annot)
totalCells <- table(fac)
cc2col <- function(cc, rate=15, base=0.001){
  cc <- round((cc[2:3]/totalCells)/sum(cc[2:3]/totalCells), 2)
  cv <- cc
  cv <- dexp(cv, rate)
  cv <- cv/rate * (1-base)
  col <- adjustcolor(rgb(cv[1],cv[2], 1), offset = c(0.1, 0.1, 0.1,0.1))
  return(col)
}

cbm <- function(d,fac) {
  if(is.leaf(d)) {
    #lc <- fac[leafContent[[as.numeric(attr(d,'label'))]]]
    cc <- attr(d, "cc")
    col <- cc2col(cc)
    attr(d,"edgePar") <- c(attr(d,"edgePar"),list(col=col))
    return(d);
  } else {
    oa <- attributes(d);
    d <- lapply(d,cbm,fac=fac);
    attributes(d) <- oa;
    cc <- attr(d, "cc")
    col <- cc2col(cc)
    attr(d,"edgePar") <- c(attr(d,"edgePar"),list(col=col))
    return(d);
  }
  
  
}
dend <- cbm(dend,species_annot)
tree <- dend


plot(tree)
dev.off()
table(comb.integrated$integrated_clusters,comb.integrated$seurat_clusters.new)
levels(comb.integrated$integrated_clusters)

############################################# Confusion heatmap  #############################################################################
##############################################################################################################################################
library(Seurat)
library(dplyr)
library(Matrix)
library(matrixStats)
library(feather)
sample.combined = comb.integrated





sample.combined$label_for_heatmap <- paste(sample.combined$species, 
                                           sample.combined$cluster_label.new, sep = "_")

reorder_matrix <- function(matrix1, by.rows = TRUE) {
  if (by.rows == TRUE) {
    conf.order <- order(apply(matrix1, 1, which.max))
    matrix1.reordered <- matrix1[conf.order, ]
  } else {
    conf.order <- order(apply(matrix1, 2, which.max))
    matrix1.reordered <- matrix1[, conf.order]
  }
}

compare_cl <- function(cl, ref.cl,
                       plot.title = NA, plot.silent = TRUE,
                       heat.colors = colorRampPalette(c("white", "grey70", "black"))(100),
                       row.cl.num = min(length(unique(cl)),
                                        length(unique(ref.cl)))) {
  library(grid)
  library(pheatmap)
  
  conf1 <- table(cl, ref.cl)
  conf1 <- sweep(conf1, 1, rowSums(conf1), "/")
  conf2 <- reorder_matrix(conf1)
  
  # Cluster co-occurence
  cl.prop.cocl <- apply(conf1, 2, function(x) {
    grid1 <- expand.grid(x, x)
    min.prop <- apply(grid1, 1, min)
  })
  
  cl.prop.cocl.total <- apply(cl.prop.cocl, 1, sum)
  cl.prop.cocl.m <- matrix(cl.prop.cocl.total, nrow(conf1), nrow(conf1),
                           dimnames = list(rownames(conf1), rownames(conf1)))
  
  ph1 <- pheatmap(conf2, cutree_rows = row.cl.num, clustering_method = "ward.D2",
                  # annotation_row = ref.cl.anno[, -grep("cluster_label", colnames(ref.cl.anno))],
                  color = heat.colors, fontsize = 6,
                  main = plot.title, silent = plot.silent)
  return(list(conf = conf2, cocl = cl.prop.cocl.m, ph = ph1))
}

# Heatmap palette
#heat.colors <- colorRampPalette(c("grey99", "orange", "red"))(100)
#heat.colors <- colorRampPalette(c("white", "#95f0f0", "#00fcfc"))(100)
heat.colors <- colorRampPalette(c("white", "grey70", "black"))(100)
#"white", "grey70", "black"
# Compare clustering
ref.cl <- sample.combined$label_for_heatmap
cca.cl <- sample.combined$integrated_clusters
compare.species <- unique(sample.combined$species)
dend.order <- labels(dend)
write.csv(dend.order,file = "dendrogram_order_cochlea.csv")
cl.conf <- compare_cl(ref.cl, cca.cl)
cocl <- cl.conf$conf
a <- as.data.frame.array(cl.conf[["conf"]])
write.csv(a, file = "cell_proportion_integrated_clusters.csv")
a<- read.csv("cell_proportion_integrated_clusters.csv",row.names = 1)# after re-arrangement

library(fpc)
clus.method <- "single"

clus.num <- pamk(cocl, 1:(min(nrow(cocl), ncol(cocl)) - 1))$nc

ph1 <- pheatmap(a, clustering_method = clus.method,
                cutree_cols = clus.num, cutree_rows = clus.num,
                color = heat.colors,
                fontsize = 6)

pheatmap(a, cluster_rows = FALSE, cluster_cols = FALSE,
         color = heat.colors,
         fontsize = 15) 
dev.off()
################################################################################################
################################################################################################
################################################################################################
#METANEIGHBOR ANALYSIS##########################################################################
################################################################################################
library(Seurat)
library(gridExtra)
library(ggplot2)
library(sctransform)
library(tidyverse)
library(dplyr)
library(Matrix)
library(matrixStats)
library(gplots)
library(ggplot2)
library(feather)
library(SingleCellExperiment)
options(stringsAsFactors = FALSE)
sce <- readRDS("C:/Users/Dell/Desktop/WORKPLACE/R/Cochlea_integrated.rds")
sce=cochlea.combined.sct2
sce=Cochlea_integrated
DefaultAssay(sce) <- "SCT"
sce$class_label = "Coch"
table(sce$class_label)
sce$sample_id <- colnames(sce)
sce$study_id <- sce$species# for each replicate
#sce$study_id <- "human"# as uniform
sce$joint_cluster_label <- sce$integrated_clusters
library(scater)
sce <- as.SingleCellExperiment(sce)
sce
library(MetaNeighbor)
dim(sce)
head(colData(sce))
table(sce$integrated_clusters,sce$study_id)
global_hvgs <- variableGenes(dat=sce,exp_labels=sce$study_id)
length(global_hvgs)
global_hvgs2=global_hvgs[1:250]
aurocs <- MetaNeighborUS(var_genes = global_hvgs2,
                         
                         dat = sce,
                         
                         study_id = sce$study_id,
                         
                         cell_type = sce$integrated_subclasses,#CHANGE TO INTEGRATED_CLUSTERS
                         
                         fast_version=TRUE)
plotBPlot(head(aurocs, 50))
library(viridis)
pheatmap::pheatmap(aurocs,color = viridis(100),border_color = FALSE)
pheatmap::pheatmap(aurocs,border_color = FALSE)
plotHeatmap(aurocs,cex=1)
dev.off()
top_hits = topHits(aurocs,
                   
                   dat = sce,
                   
                   study_id = sce$study_id,
                   
                   cell_type = sce$integrated_clusters,
                   
                   threshold = 0.9)
best_hits = MetaNeighborUS(var_genes = global_hvgs,
                           dat = sce,
                           study_id = sce$study_id,
                           cell_type = sce$integrated_clusters,
                           fast_version = TRUE,
                           one_vs_best = TRUE, symmetric_output = FALSE)
plotHeatmap(best_hits, cex = 0.8)
dev.off()
#####functional characterization of integrated clusters

go_sets = readRDS("go_human.rds")##produced in figure1 draft
dim(sce)
known_genes = rownames(sce)
go_sets = lapply(go_sets, function(gene_set) {
  gene_set[gene_set %in% known_genes]
})
min_size = 10
max_size = 100
go_set_size = sapply(go_sets, length)
go_sets = go_sets[go_set_size >= min_size &
                    go_set_size <= max_size]
length(go_sets)

aurocs = MetaNeighbor(dat = sce,
                      experiment_labels = sce$study_id,
                      celltype_labels = sce$subclasses_label,
                      genesets = go_sets,
                      fast_version = TRUE, bplot = FALSE, batch_size = 50)
table(sce$seurat_clusters.new,sce$integrated_clusters)
aurocs <- as.data.frame(aurocs)
plotDotPlot(dat = sce,
            experiment_labels = sce$study_id,
            celltype_labels = sce$subclasses_label,
            gene_set = go_sets[["GO:0050982|detection of mechanical stimulus|BP"]])
write.csv(aurocs,file = "GO_INTEGRATED_DATA_SUBCLASSES_LABEL.csv")
sce$seurat_clusters.new.integrated_cluster <- paste(sce$seurat_clusters.new,sce$integrated_clusters,sep = "_")
aurocs = MetaNeighbor(dat = sce,
                      experiment_labels = sce$study_id,
                      celltype_labels = sce$seurat_clusters.new.integrated_cluster,
                      genesets = go_sets,
                      fast_version = TRUE, bplot = FALSE, batch_size = 50)
write.table(aurocs, "functional_aurocs_integrated_clusters.txt")
aurocs2 <- as.data.frame(aurocs)
write.csv(aurocs2,file = "GO_INTEGRATED_DATA_integrated_cluster.csv")
plotBPlot(head(aurocs, 100))
gs_size = sapply(go_sets, length)
aurocs_df = data.frame(go_term = rownames(aurocs), aurocs)
aurocs_df$average = rowMeans(aurocs)
aurocs_df$n_genes = gs_size[rownames(aurocs)]
write.table(aurocs_df, "functional_aurocs_integrated_clusters2.txt")
write.csv(aurocs_df, "functional_aurocs_integrated_clusters2.csv")
head(aurocs_df[order(aurocs_df$average, decreasing = TRUE),],10)
small_sets = aurocs_df[aurocs_df$n_genes < 20,]
head(small_sets[order(small_sets$average, decreasing = TRUE),],10)
plotDotPlot(dat = sce,
            experiment_labels = sce$study_id,
            celltype_labels = sce$seurat_clusters.new,
            gene_set = go_sets[["GO:0002093|auditory receptor cell morphogenesis|BP"]])

plotDotPlot(dat = sce,
            experiment_labels = sce$study_id,
            celltype_labels = sce$integrated_subclasses,
            gene_set = go_sets[["GO:0051963|regulation of synapse assembly|BP"]])

###############################################################################
Cochlea_integrated <- readRDS("C:/Users/Dell/Desktop/WORKPLACE/R/Cochlea_integrated.rds")
table(Cochlea_integrated$seurat_clusters.new,Cochlea_integrated$cluster_label)
Cochlea_integrated$cluster_label.species <- paste(Cochlea_integrated$cluster_label,
                                                  Cochlea_integrated$species,sep = "_")
table(Cochlea_integrated$cluster_label.species)
Cochlea_integrated$integrated_clusters.species <- paste(Cochlea_integrated$integrated_clusters,
                                                        Cochlea_integrated$species,sep = "_")
Idents(Cochlea_integrated) <- "species"
Cochlea_integrated_mouse <- subset(Cochlea_integrated,idents="mouse")
Cochlea_integrated_human <- subset(Cochlea_integrated,idents="human")
mouse_regions <- table(Cochlea_integrated_mouse$integrated_clusters, 
                       Cochlea_integrated_mouse$cluster_label.species)
mouse_regions <- as.data.frame(mouse_regions)
mouse_regions$species <- "mouse"
colnames(mouse_regions) <- c("integrated_clusters","region","frequency","species")
human_regions <- table(Cochlea_integrated_human$integrated_clusters, 
                       Cochlea_integrated_human$cluster_label.species)
human_regions <- as.data.frame(human_regions)
human_regions$species <- "human"
colnames(human_regions) <- c("integrated_clusters","region","frequency","species")
alluvial_table <- rbind(human_regions,mouse_regions)

mouse_regions$mouse_regions <- mouse_regions$region
mouse_regions$human_regions <- NA
human_regions$human_regions <- human_regions$region
human_regions$mouse_regions <- NA
alluvial_table_species <- rbind(human_regions,mouse_regions)
write.csv(alluvial_table_species, file = "cochlea_integrated_CCA_alluvial_table.csv")

alluvial_table_species <- as.data.frame(alluvial_table_species)
alluvial_table2 <- alluvial_table[,c(1,2,3)]
head(alluvial_table2)
library(ggsankey)
#remotes::install_github("davidsjoberg/ggsankey")
library(ggplot2)
library(cols4all)
#BiocManager::install("cols4all")

nodes <- data.frame(name = c(as.character(alluvial_table_species$region),
                             as.character(alluvial_table_species$integrated_clusters)) %>% unique()) #unique去重
head(nodes)
alluvial_table_species$IDintegrated_clusters <- match(alluvial_table_species$integrated_clusters, 
                                                      nodes$name) -1 #ID必须从0开始索引，所以统一减去1
alluvial_table_species$IDregion <- match(alluvial_table_species$region, nodes$name) -1
alluvial_table_species$IDmouse_region <- match(alluvial_table_species$mouse_regions, nodes$name) -1
alluvial_table_species$IDhuman_region <- match(alluvial_table_species$human_regions, nodes$name) -1
head(alluvial_table_species)

hum_source <- match(alluvial_table_species$IDhuman_region, alluvial_table_species$IDregion)-1
hum_source <- hum_source[!is.na(hum_source)]
mus_source <- (match(alluvial_table_species$IDmouse_region, alluvial_table_species$IDregion))-1
mus_source <- mus_source[!is.na(mus_source)]

value <- as.numeric(alluvial_table_species$frequency)
target <- as.numeric(match(alluvial_table_species$integrated_clusters,  nodes$name))-1

target_cropped <- target[1:length(hum_source)]
target_cropped2 <- target[c(length(hum_source)+1):length(target)]
source_new <-  as.numeric(c(hum_source, target_cropped2))
target_new <-  as.numeric(c(target_cropped, mus_source))
data <- data.frame(value = value,target=target_new,source=source_new)
data_filtered <- data[data$value > 0,]

clusters <- list(names = data.frame(alluvial_table_species$IDregion),data=data_filtered)

saveRDS(clusters, "clusters_network3D_mouse_human_regions.rds")

library(networkD3)
clusters <- readRDS("clusters_network3D_mouse_human_regions.rds")
nodes <- data.frame(name = c(as.character(alluvial_table_species$IDregion),
                             as.character(alluvial_table_species$IDintegrated_clusters)) %>% unique()) #unique去重
head(nodes)
sankeyNetwork(Links = clusters$data, Nodes = nodes, Source = "source", Target = "target",
              Value = "value", NodeID = "name", units = "cells", fontSize = 8, 
              nodeWidth = 40, nodePadding = 10, height = 800, width = 800, 
              iterations = 10000, sinksRight = F)


p5 <- sankeyNetwork(Links = clusters$data,
                    Nodes = nodes,
                    Source = "source",units = "cells",
                    Target = "target",
                    Value = "value",
                    NodeID = "name",
                    sinksRight = F,
                    nodeWidth = 40, #节点宽度
                    fontSize = 15,
                    nodePadding = 10) #节点间间隙
p5

dev.off()
#######################################################################################
#######################################################################################
##########Gene ontology analysis of conserved cluster markers
## use gProfiler
#For the gene ontology analysis of maker genes we first calculated conserved marker genes for 
#integrated_subclasses using Seurat’s FindConservedMarkers 
#(grouping.var = "species", assay="SCT", only.pos= TRUE). We removed genes that were not 
#detected in at least 20% of a single cluster in both species from the marker genes list. We used 
#gprofiler2 (69) to calculate enriched GO-terms (gost: user_threshold = 0.01, correction_method = 
#"g_SCS"), as background we used the intersection of expressed one-to-one orthologs between 
#mouse and human
sce <- readRDS("C:/Users/Dell/Desktop/WORKPLACE/R/Cochlea_integrated.rds")
table(sce$integrated_subclasses,sce$species)
table(sce$species,sce$subclasses_label.new)
Idents(sce) <- "integrated_subclasses"
table(sce$species,sce$integrated_subclasses)
Idents(sce) <- "integrated_clusters"
table(sce$species,sce$integrated_clusters)
brain_integrated_markers_conc <- matrix(,0,14)

cluster_list <- c("BC_1","Fibro1_1","TBC_1","MC_1","BC_2","GC_1","TBC_2","SMC_1",
                  "DC_PC_HeC_1","IC_1","Fibro4","DC_PC_HeC_2","BC_3",
                  "Immun","DC_PC_HeC_3","SMC_2","SGN_1","SGN_2","Fibro3_1",
                  "Fibro2","Fibro1_2","IPh_IBC_1","Fib","IPh_IBC_2","GC_2",
                  "Fibro3_2","Fibro1_3","HC_1","MC_2","IC_2","BC_4",
                  "Fibro3_3","IPh_IBC_4","SMC_3","DC_PC_HeC_4","SGN_3","HC_2",
                  "SGN_4","SGN_5","IPh_IBC_5")# switch to integrated subclasses
cluster_list
sce <- PrepSCTFindMarkers(sce)
library(metap)
for ( i in 1:length(cluster_list)){
  tmp <- i
  brain_integrated_markers <- FindConservedMarkers(sce, 
                                                   ident.1=cluster_list[tmp], grouping.var = "species", 
                                                   verbose = FALSE, assay="SCT", only.pos= TRUE)
  brain_integrated_markers$cluster <- ""

  
  brain_integrated_markers$cluster <- cluster_list[tmp]
  brain_integrated_markers$gene <- rownames(brain_integrated_markers)
  brain_integrated_markers_conc <- rbind(brain_integrated_markers_conc, brain_integrated_markers)
}
write.csv(brain_integrated_markers_conc, file = paste0("Data_human-mouse_conserved_markers_table_integrated_clusters.csv"))
library(Rfast)
NoisyGenes <- function(object, min.pct=0.2, clusters.use) {
  clusters <- clusters.use
  genes.use = rownames(object)
  object.raw.data <- as.matrix(GetAssayData(object, slot="counts"))
  pct.matrix = matrix(data=NA, nrow=length(genes.use), ncol=length(clusters))
  rownames(pct.matrix) <- genes.use
  colnames(pct.matrix) <- clusters
  thresh.min=0
  for (i in clusters){
    cells.cluster <- WhichCells(object=object, idents=i)
    data.cluster <- object.raw.data[,colnames(object.raw.data) %in% cells.cluster]
    pct.cluster <- round(apply(object.raw.data[genes.use, cells.cluster, drop = F],1,function(x)return(length(x[x>thresh.min])/length(x))),3)
    pct.matrix[,i] <- pct.cluster
  }
  pct.max <- rowMaxs(pct.matrix)
  names(pct.max) <- genes.use
  noisy.genes <- names(pct.max[pct.max < min.pct])
  return(noisy.genes)
}
noisy.liz <- NoisyGenes(sce, 0.2, levels(Idents(sce)))

bi_unique_conserved_markers <- unique(brain_integrated_markers_conc$gene )


diffexpr <- bi_unique_conserved_markers

diffexpr <- diffexpr[!diffexpr %in% noisy.liz]
write.csv(as.data.frame(diffexpr), file = "Data_human-mouse_conserved_markers_table_new_remove_replicates.csv")
library(gprofiler2)

gostres <- gost(query = diffexpr,
                organism = "hsapiens", ordered_query = FALSE,
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                measure_underrepresentation = FALSE, evcodes = FALSE,
                user_threshold = 0.01, correction_method = "g_SCS",
                domain_scope = "annotated", custom_bg = NULL,
                numeric_ns = "",  as_short_link = FALSE,sources = c("GO:MF", "GO:BP","GO:CC"))

gostplot(gostres, capped = FALSE, interactive = TRUE)
result_df <- gostres$result
result_df$GO_term <- paste(result_df$term_id,result_df$term_name)

result_df_top10 <- result_df %>% group_by(source) %>% top_n(10, -p_value)
result_df_top10$GO_term <- factor(result_df_top10$GO_term, 
                                  levels = result_df_top10$GO_term[order(-log10(result_df_top10$p_value))])
ggplot(result_df_top10, aes(x=GO_term, y=-log10(p_value), fill=source)) +     geom_bar(stat="identity")  + theme(axis.text.x = element_text(angle = 270,hjust=0))
dev.off()
ggplot(result_df_top10, aes(x=GO_term, y=-log10(p_value), fill=source)) +     geom_bar(stat="identity") + theme_bw() + theme(axis.text.x = element_text(angle = 270,hjust=0))

dev.off()

result_mf_top10 <- result_df_top10[result_df_top10$source=="GO:MF",]
result_cc_top10 <- result_df_top10[result_df_top10$source=="GO:CC",]
result_bp_top10 <- result_df_top10[result_df_top10$source=="GO:BP",]

ggplot(result_mf_top10, aes(x=GO_term, y=-log10(p_value))) +     geom_bar(stat="identity") + theme_bw() + theme(axis.text.x = element_text(angle = 270,hjust=0))
dev.off()
ggplot(result_cc_top10, aes(x=GO_term, y=-log10(p_value))) +     geom_bar(stat="identity") + theme_bw() + theme(axis.text.x = element_text(angle = 270,hjust=0))
dev.off()
ggplot(result_bp_top10, aes(x=GO_term, y=-log10(p_value))) +     geom_bar(stat="identity") + theme_bw() + theme(axis.text.x = element_text(angle = 270,hjust=0))
dev.off()

tf_gene_list <- read.csv("human_tf_gene_list.csv",row.names = 1)

brain_integrated_markers_conc$HGNC.symbol <- brain_integrated_markers_conc$gene
conserved_tf <- merge(x=brain_integrated_markers_conc,y=tf_gene_list, by="HGNC.symbol", all=FALSE)

write.csv(conserved_tf,file = "conserved_tf_across_human_and_mouse_new.csv")

DimPlot(sce,reduction = "umap",split.by = "species")
#################hierarchical clustering of integrated clusters based on TFs of conserved marker genes

library(pvclust)
library(parallel)
library(dendextend)
library(future)
plan("multicore", workers = 4)
options(future.globals.maxSize= 62914560000) # 60GB (60000*1024^2)
library(Seurat)
library(dplyr)
library(cowplot)
library(zoo)
library(ggplot2)
library(scater)
library(stringr)
library(gplots)
library(matrixStats)
sce <- readRDS("C:/Users/Dell/Desktop/WORKPLACE/R/Cochlea_integrated.rds")
#sce <- SCTransform(sce, method = "glmGamPoi")
Idents(sce) <- "integrated_subclasses"

DefaultAssay(sce)<-"SCT"
neurons_avg <- AverageExpression(sce)
neurons_avg_data <- neurons_avg[["SCT"]]
NoisyGenes <- function(object, min.pct=0.2, clusters.use) {
  clusters <- clusters.use
  genes.use = rownames(object)
  object.raw.data <- as.matrix(GetAssayData(object, slot="counts"))
  pct.matrix = matrix(data=NA, nrow=length(genes.use), ncol=length(clusters))
  rownames(pct.matrix) <- genes.use
  colnames(pct.matrix) <- clusters
  thresh.min=0
  for (i in clusters){
    cells.cluster <- WhichCells(object=object, idents=i)
    data.cluster <- object.raw.data[,colnames(object.raw.data) %in% cells.cluster]
    pct.cluster <- round(apply(object.raw.data[genes.use, cells.cluster, drop = F],1,function(x)return(length(x[x>thresh.min])/length(x))),3)
    pct.matrix[,i] <- pct.cluster
  }
  pct.max <- rowMaxs(pct.matrix)
  names(pct.max) <- genes.use
  noisy.genes <- names(pct.max[pct.max < min.pct])
  return(noisy.genes)
}
conserved_TFs <- read.csv("conserved_tf_across_human_and_mouse_for_dendrogram_new.csv")
#Human_TFs <- read.csv("human_tf_gene_list2.csv")

noisy.liz <- NoisyGenes(sce, 0.2, levels(Idents(sce)))#0.25
neurons_avg_data_filtered <- neurons_avg_data[!rownames(neurons_avg_data) %in% noisy.liz,]
neurons_avg_TFs <- neurons_avg_data_filtered[rownames(neurons_avg_data_filtered) %in% conserved_TFs$gene,]
tfs <- as.data.frame(neurons_avg_TFs)

# 
neurons_avg_TFs2 <- ScaleData(neurons_avg_TFs)
pheatmap::pheatmap(neurons_avg_TFs2,scale = "row")
dist_func <- function(x){x <- as.dist(1-cor(x, method="s"))}
pvclust.TFs_dist <- pvclust(neurons_avg_TFs, method.dist=dist_func, 
                            method.hclust="ward.D2", nboot=10000, parallel=T)

plot(pvclust.TFs_dist)
pvrect(pvclust.TFs_dist, alpha=0.9)
dev.off()

dend <- as.dendrogram(pvclust.TFs_dist)
dend %>% pvclust_show_signif_gradient(pvclust.TFs_dist, signif_col_fun = colorRampPalette(c("lightgrey","grey","black")), signif_type="au") %>%
  plot(main = "neuronal clusters branch confidence")
dend <- rotate(dend)
dendrogram_names <- labels(dend)
dend <- as.dendrogram(pvclust.TFs_dist)
dend %>%
  pvclust_show_signif(pvclust.TFs_dist, signif_value = c("black", "black"), show_type = "col") %>%
  plot(main = "Cluster dendrogram with AU/BP values (%)")
pvrect2(pvclust.TFs_dist, alpha=0.90)
dendrogram_names <- labels(dend)
levels(sce) <- dendrogram_names
genes <- c("IKZF2","ZNF516","ISL1","ZMAT4","SIX2",#OHC
           "HDAC9","ZNF410","BPNT1","TBX2","BTG2",#IHC
           "NME2","EIF3K","ID3","SOX4","RPL7",#DCOPC_2
           "DEPTOR","ATOH8","IRX6","BCL6","ZNF37A",#DCOPC_3
           "MXD1","PER1","HIF1A","JUN","HES1",#IPHIBC_1
           "CTNNBIP1","ID1","EGR1","HEY2","HEY1",#IPHIBC_2
           "SMAD9","POLK","BACH1","TSC22D1","SMAD1",#DCOPC_1
           "TSC22D2","SOX6","MEIS2","ZNF609","DACH1",#HEC
           "PLAGL1","POT1","EDA","TNFSF4","CUX1",#IPHIBC_5
           "ARHGAP35","PRDM11","RPS6KA5","PHTF1","CUL5",#IPHIBC_3
           "TRPS1","ETV5","PRKD2","RORB","ZNF185",#IPHIBC_4
           "ARNT2","ZNF326","RAPGEF4","IKBKB","LGR4",#PC_1
           "EPAS1","ZBTB16","SMO","TAF4B","KLF13",#PC_3
           "SOX9","AEBP1","BHLHE41","POU3F3","POU3F4",#PC_2
           "RPL7A","RPS27A","PRKCH","JAK2","MEF2A",#TBC_1
           "IRX5","IRX3","CREB5","EPHA5","MET"#TBC_2
)
tfs$gene <- rownames(tfs)
DefaultAssay(sce) <- "SCT"
DotPlot(sce, features = genes,cols = c("white","#284aaa"))+theme(axis.text.x = element_text(angle = 90, hjust = 1))
DotPlot(sce, features = tfs$gene,cols = c("white","#284aaa")
)+theme(axis.text.x = element_text(angle = 90, hjust = 1))
pdf(file="200_TFS_conserved.pdf", width=100, height=15)
DotPlot(sce, features = tfs$gene,
        dot.scale = 8,cols  =c("white","#284aaa")
) + RotatedAxis()
dev.off(
  
)
genes <- c("EGLN1","EPAS1","RAPGEF4","CFLAR","MECOM","HIVEP3","PAX3","MITF",
           "MET","ZEB2","KAT2B","RELN","GAS7","ZNF536","PLEKHA4","REL","PLSCR1",
           "TLR2","HMGB2","ZNF451","EZH2","TOP2A","PAWR","ZFHX3","NOTCH3","POU3F1",
           "MAF","MKRN1","NCOA4","SATB2","MEF2C","RUNX2","ZNF326","WNT5A","LGR4",
           "PLPP3","DDR2","TBX18","ZNF385D","CREB5","BNC2","ZIC1","CAMTA1","FOS",
           "MYT1L","ZNF804A","EPHA5","TLE4","PBX3","TSHZ3","PBX1","ZHX2","GATA3",
           "SOX6","ADAMTS17","ZHX3","BMP7","HDAC9","NRG1","ZMAT4","RORB","ZFP64",
           "IKZF2","LMX1A","ALK","MEIS1","SGK3","CUX2","DACH1","ESRRG"
           
           
)

DotPlot(sce, features = genes,
        dot.scale = 8,cols  =c("white","#284aaa")
)+ RotatedAxis()
dev.off(
  
)

DimPlot(sce,reduction = "umap",split.by = "species",label = TRUE)
DefaultAssay(sce) <- "SCT"
FeaturePlot(sce,features = "HDAC9",split.by = "species",label = TRUE)

DefaultAssay(sce) <- "SCT"
sce$class_label = "Coch"
table(sce$class_label)
sce$sample_id <- colnames(sce)
sce$study_id <- sce$species# for each replicate
#sce$study_id <- "human"# as uniform
sce$joint_cluster_label <- sce$integrated_clusters
library(scater)
sce2 <- as.SingleCellExperiment(sce)
sce2
library(MetaNeighbor)
dim(sce2)
head(colData(sce2))
table(sce2$integrated_clusters,sce2$study_id)

noisy.liz <- NoisyGenes(sce, 0.2, levels(Idents(sce)))#0.25
neurons_avg_data_filtered <- neurons_avg_data[!rownames(neurons_avg_data) %in% noisy.liz,]
neurons_avg_TFs <- neurons_avg_data_filtered[rownames(neurons_avg_data_filtered) %in% conserved_TFs$gene,]
tfs <- as.data.frame(neurons_avg_TFs)
tfs$gene <- rownames(tfs)

#we used conserved tfs as hvgs for hierarchical clustering
global_hvgs2=tfs$gene
aurocs <- MetaNeighborUS(var_genes = global_hvgs2,
                         
                         dat = sce2,
                         
                         study_id = sce$study_id,
                         
                         cell_type = sce$integrated_clusters,#CHANGE TO INTEGRATED_CLUSTERS
                         
                         fast_version=TRUE)
dist_func <- function(x){x <- as.dist(1-cor(x, method="s"))}
pvclust.TFs_dist <- pvclust(aurocs, method.dist=dist_func, 
                            method.hclust="ward.D2", nboot=10000, parallel=T)

plot(pvclust.TFs_dist)
pvrect(pvclust.TFs_dist, alpha=0.9)
dend <- as.dendrogram(pvclust.TFs_dist)
dend %>% pvclust_show_signif_gradient(pvclust.TFs_dist, signif_col_fun = colorRampPalette(c("lightgrey","grey","black")), signif_type="au") %>%
  plot(main = "neuronal clusters branch confidence")
dend <- rotate(dend)
dendrogram_names <- labels(dend)
dend <- as.dendrogram(pvclust.TFs_dist)
dend %>%
  pvclust_show_signif(pvclust.TFs_dist, signif_value = c("black", "black"), show_type = "col") %>%
  plot(main = "Cluster dendrogram with AU/BP values (%)")
dev.off()
library(viridis)
pheatmap::pheatmap(aurocs,color = viridis(100),border_color = FALSE)
pheatmap::pheatmap(aurocs,border_color = FALSE)
plotHeatmap(aurocs,cex=1)
dev.off()
################################################################################################
sce <- readRDS("C:/Users/Dell/Desktop/WORKPLACE/R/Cochlea_integrated.rds")

Idents(sce) <- "species"
DimPlot(sce,reduction = "umap",split.by = "species",label = TRUE,cols = c("#308ce1", "#f28eac"))
DefaultAssay(sce) <- "RNA"
sce <- SCTransform(sce, vars.to.regress = "nCount_RNA")
Idents(sce) <- "integrated_subclasses"
FeaturePlot(sce,features = "ZBTB20",split.by = "species",label = FALSE,pt.size = 1)
FeaturePlot(sce,features = "DACH1",split.by = "species",label = FALSE,pt.size = 1)
VlnPlot(sce, features = "DACH1",split.by = "species",assay = "SCT")
FeaturePlot(sce,features = "IKZF2",split.by = "species",label = FALSE,pt.size = 1)
FeaturePlot(sce,features = "RORB",split.by = "species",label = FALSE,pt.size = 1)
FeaturePlot(sce,features = "ESRRG",split.by = "species",label = FALSE,pt.size = 1)
FeaturePlot(sce,features = "MYT1L",split.by = "species",label = FALSE,pt.size = 1)
FeaturePlot(sce,features = "PROX1",split.by = "species",label = FALSE,pt.size = 1)
FeaturePlot(sce,features = "GATA3",split.by = "species",label = FALSE,pt.size = 1)
FeaturePlot(sce,features = "ZEB2",split.by = "species",label = FALSE,pt.size = 1)
FeaturePlot(sce,features = "PRDM16",split.by = "species",label = FALSE,pt.size = 1)
FeaturePlot(sce,features = "RFX3",split.by = "species",label = FALSE,pt.size = 1)
dev.off()







