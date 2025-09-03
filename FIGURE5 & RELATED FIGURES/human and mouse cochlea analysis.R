rm(list=ls())
#dev.off()
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
#library(loomR)
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
library(SingleCellExperiment)
library(Matrix)
#####################################################
#transfer H5AD into SEURAT OBJECT

library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(patchwork)
library(zellkonverter)
#BiocManager::install("zellkonverter")
#BiocManager::install("loomR",force = TRUE)
library(loomR)
library(Seurat)
library(SeuratDisk)
library(data.table)
adata_loom <- connect(filename = "human_cochlea_late_2nd_trimeste.loom",
                      mode = "r+",skip.validate = TRUE)
matrix=adata_loom[["matrix"]][,]
matrix=t(matrix)
dim(matrix)
gene = adata_loom$row.attrs$var_names[]
barcode = adata_loom$col.attrs$obs_names[]

meta_data = read.csv('human_cochlea_late_2nd_trimeste_obs.csv',row.names = 1) # as form as dataframe format
meta_feature = read.csv('human_cochlea_late_2nd_trimeste_var.csv',row.names = 1)

colnames(matrix)= barcode
row.names(matrix)= gene
#x_scvi = adata_loom$col.attrs$X_scVI[,]
x_umap=read.csv("human_cochlea_late_2nd_trimester_umap.csv",row.names = 1)


seurat_object_human= CreateSeuratObject(counts = matrix,meta.data = meta_data,
                                  project = 'human_cochlea_loom',
                                  min.cells = 0, 
                                  min.features = 0)
seurat_object_human@assays[["RNA"]]@meta.features <- meta_feature

rownames(x_umap) = barcode
colnames(x_umap) = c('UMAP_1','UMAP_2')


seurat_object_human$cellid <- rownames(seurat_object_human@meta.data)
head(seurat_object_human$cellid)
rownames(x_umap)=seurat_object_human$cellid


seurat_object_human <-FindVariableFeatures(seurat_object_human)
seurat_object_human <- ScaleData(seurat_object_human)
seurat_object_human <- RunPCA(seurat_object_human)
seurat_object_human<- RunUMAP(seurat_object_human,dims = 1:30)
seurat_object_human <- FindNeighbors(seurat_object_human,dims = 1:30)
seurat_object_human <- FindClusters(seurat_object_human,resolution = 0.8)
DimPlot(seurat_object_human,reduction = "umap",group.by = "cell_type")

seurat_object_human@reductions[["umap"]]@cell.embeddings=as.matrix(x_umap)
DimPlot(seurat_object_human,reduction = "umap",group.by = "cell_type")

sce_human = as.SingleCellExperiment(seurat_object_human)
colnames(x_umap) = c('umap_1','umap_2')
sce_human@int_colData@listData[["reducedDims"]]@listData[["umap"]] = x_umap
#sce1@int_colData@listData[["reducedDims"]]@listData[["scANVI"]] = x_scANVI
seurat_object_human <- as.Seurat(sce_human)
#seurat_object <- RunUMAP(seurat_object,reduction = "scANVI",dims = 1:20)
DimPlot(seurat_object_human,reduction = "umap",group.by = "cell_type")
DimPlot(seurat_object_human,reduction = "umap",group.by = "sample")
dev.off()
library(future)
plan("multicore", workers = 4)
options(future.globals.maxSize= 62914560000) # 60GB (60000*1024^2)
####################################################################################################
############for mouse species###########################################################
#####################################################################################
adata_loom <- connect(filename = "mouse_cochlea_P25_P28.loom",
                      mode = "r+",skip.validate = TRUE)
matrix=adata_loom[["matrix"]][,]
matrix=t(matrix)
dim(matrix)
gene = adata_loom$row.attrs$var_names[]
barcode = adata_loom$col.attrs$obs_names[]

meta_data = read.csv('mouse_cochlea_P25_P28_obs.csv',row.names = 1) # as form as dataframe format
meta_feature = read.csv('mouse_cochlea_P25_P28_var.csv',row.names = 1)

colnames(matrix)= barcode
row.names(matrix)= gene
#x_scvi = adata_loom$col.attrs$X_scVI[,]
x_umap=read.csv("mouse_cochlea_P25_P28_umap.csv",row.names = 1)


seurat_object_mouse= CreateSeuratObject(counts = matrix,meta.data = meta_data,
                                  project = 'mouse_cochlea_loom',
                                  min.cells = 0, 
                                  min.features = 0)
seurat_object_mouse@assays[["RNA"]]@meta.features <- meta_feature

rownames(x_umap) = barcode
colnames(x_umap) = c('UMAP_1','UMAP_2')


seurat_object_mouse$cellid <- rownames(seurat_object_mouse@meta.data)
head(seurat_object_mouse$cellid)
rownames(x_umap)=seurat_object_mouse$cellid


seurat_object_mouse <-FindVariableFeatures(seurat_object_mouse)
seurat_object_mouse <- ScaleData(seurat_object_mouse)
seurat_object_mouse <- RunPCA(seurat_object_mouse)
seurat_object_mouse<- RunUMAP(seurat_object_mouse,dims = 1:30)
seurat_object_mouse <- FindNeighbors(seurat_object_mouse,dims = 1:30)
seurat_object_mouse <- FindClusters(seurat_object_mouse,resolution = 0.8)
DimPlot(seurat_object_mouse,reduction = "umap",group.by = "cell_type")

seurat_object_mouse@reductions[["umap"]]@cell.embeddings=as.matrix(x_umap)
DimPlot(seurat_object_mouse,reduction = "umap",group.by = "cell_type")

sce_mouse = as.SingleCellExperiment(seurat_object_mouse)
colnames(x_umap) = c('umap_1','umap_2')
sce_mouse@int_colData@listData[["reducedDims"]]@listData[["umap"]] = x_umap
#sce1@int_colData@listData[["reducedDims"]]@listData[["scANVI"]] = x_scANVI
seurat_object_mouse <- as.Seurat(sce_mouse)
#seurat_object <- RunUMAP(seurat_object,reduction = "scANVI",dims = 1:20)
DimPlot(seurat_object_mouse,reduction = "umap",group.by = "cell_type")
DimPlot(seurat_object_mouse,reduction = "umap",group.by = "batch")
dev.off()
#################################################################################
########################################################################
#####################################################################################
###########cumulative plot and dendrogram based on TFs##############################
#####################################################################################
##for mouse species

library(Seurat)
library(matrixStats)
library(ggplot2)

Idents(seurat_object_mouse) <- "cell_type"
DimPlot(seurat_object_mouse, reduction = "umap",label = TRUE)
seurat_object_mouse <- SCTransform(seurat_object_mouse, 
            vars.to.regress = c("nCount_RNA","pct_counts_mt"),
            verbose = TRUE, method="glmGamPoi")
DefaultAssay(seurat_object_mouse) <- "SCT"
object <- seurat_object_mouse
table(Idents(seurat_object_mouse))
clusters <-levels(Idents(seurat_object_mouse))
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
  cluster_count <- length(gene_row[gene_row > 0.1])
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
  cluster_count <- length(gene_row[gene_row > 0.2])
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
  cluster_count <- length(gene_row[gene_row > 0.1])
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
  cluster_count <- length(gene_row[gene_row > 0.2])
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
     ylab = "Cumulative fraction",xlim=c(0,25),ylim=c(0,1), col="black")
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
       col=c("black", "magenta","darkgreen","blue","red","orange","darkblue","darkred"), lty=2, cex=0.7)




dev.off()
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
##for human species

library(Seurat)
library(matrixStats)
library(ggplot2)

Idents(seurat_object_human) <- "cell_type"
DimPlot(seurat_object_human, reduction = "umap",label = TRUE)
seurat_object_human <- SCTransform(seurat_object_human, 
                                   vars.to.regress = c("nCount_RNA","pct_counts_mt"),
                                   verbose = TRUE, method="glmGamPoi")


DefaultAssay(seurat_object_human) <- "SCT"
object <- seurat_object_human
table(Idents(seurat_object_human))
clusters <-levels(Idents(seurat_object_human))
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
  cluster_count <- length(gene_row[gene_row > 0.1])
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
  cluster_count <- length(gene_row[gene_row > 0.2])
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
  cluster_count <- length(gene_row[gene_row > 0.1])
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
  cluster_count <- length(gene_row[gene_row > 0.2])
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
     ylab = "Cumulative fraction",xlim=c(0,32),ylim=c(0,1), col="black")
lines(xC2F,(cumsum(yC2F)/max(cumsum(yC2F))), cex=0.8, type="o", pch=19, col="magenta")
lines(xHD,(cumsum(yHD)/max(cumsum(yHD))), cex=0.8, type="o", pch=19, col="darkgreen")
lines(xGPCRs,(cumsum(yGPCRs)/max(cumsum(yGPCRs))), cex=0.8, type="o", pch=19, col="blue")
lines(xICs,(cumsum(yICs)/max(cumsum(yICs))), cex=0.8, type="o", pch=19, col="red")
lines(xCAMs,(cumsum(yCAMs)/max(cumsum(yCAMs))), cex=0.8, type="o", pch=19, col="orange")
lines(xPGs,(cumsum(yPGs)/max(cumsum(yPGs))), cex=0.8, type="o", pch=19, col="darkblue")
lines(xRibos,(cumsum(yRibos)/max(cumsum(yRibos))), cex=0.8, type="o", pch=19, col="darkred")

legend("bottomright", legend=c("all TFs", "zinc finger-C2H2 ","homeodomain","GPCRs ","ion channels ","cell adhesion molecules ","proteoglycans ","ribosomal genes "),
       col=c("black", "magenta","darkgreen","blue","red","orange","darkblue","darkred"), lty=1, cex=0.4)


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
sce=seurat_object_mouse
sce <- SCTransform(sce, method = "glmGamPoi")
Idents(sce) <- "cell_type"
table(Idents(sce))

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

# First, make sure your cell types are properly set
Idents(sce) <- "cell_type"
table(Idents(sce))

# Set the default assay
DefaultAssay(sce) <- "SCT"

# Get average expression - make sure to return in correct format
neurons_avg <- AverageExpression(sce, assays = "SCT", group.by = "cell_type", return.seurat = FALSE)

# Extract the SCT matrix
neurons_avg_data <- neurons_avg[["SCT"]]
write.csv(neurons_avg_data,file="mouse_average_data.csv")
# Your NoisyGenes function remains the same
noisy.liz <- NoisyGenes(sce, 0.5, levels(Idents(sce)))

# Filter the data
neurons_avg_data_filtered <- neurons_avg_data[!rownames(neurons_avg_data) %in% noisy.liz,]

# Filter for TFs
neurons_avg_TFs <- neurons_avg_data_filtered[rownames(neurons_avg_data_filtered) %in% Human_TFs$HGNC.symbol,]

# Convert to data frame and ensure proper cell type names
tfs <- as.data.frame(neurons_avg_TFs)

# If column names are still numbers, you can explicitly set them:
if(all(colnames(tfs) == as.character(1:ncol(tfs)))) {
  colnames(tfs) <- levels(Idents(sce))
}
write.csv(tfs,file="mouse_tfs_for_hierarchical_clustering.csv")
############################################ 131 TFs remain##################################################
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

tfs$gene <- rownames(tfs)
DefaultAssay(sce) <- "SCT"

DotPlot(sce, features = tfs$gene,cols = c("white","#284aaa")
)+theme(axis.text.x = element_text(angle = 90, hjust = 1))
pdf(file="131_TFS_mouse.pdf", width=50, height=12)
DotPlot(sce, features = tfs$gene,
        dot.scale = 8,cols  =c("white","#284aaa")
) + RotatedAxis()
dev.off(
  
)
genes <- c("TBX2","ZNF385A","CAMTA1","RORB",  #ihc
           "YBX1","GFI1","MLXIP","GTF2I","SALL1","SIX1", "POU4F3","SIX2","ISL1","IKZF2","NR2F2", #OHC
          "PURB","NFE2L1","FLYWCH1","RUNX1",'RXRG',"KLF7","MAFB","CUX2","TSHZ3","MAF","ZNF207",
          "LMX1A","ESRRG","NR6A1","MEIS2","PRDM16","YBX2",
          "MEIS1","BHLHE41","NFIX","PROX1",
          "GATA3","NME2","NR2F1",
          "PBX1","ARHGAP35","EBF1","DPF3","THRB","TSC22D1",
          "PBX3","ZHX3","SOX5",
          "MYC","MEF2C","POU3F3","JUND","TBX18","AEBP1","ZEB2",
          "PAX3","ELF2","HIF1A","BNC2", "FOXP2",  "GLI2","ZFHX4","SMYD3","ZBTB16","ZNF521",
          "SATB2","RUNX2","EPAS1"
)

DotPlot(sce, features = genes,cols = c("white","#284aaa"))+theme(axis.text.x = element_text(angle = 90, hjust = 1))
pdf(file="MOUSE_TFS_COCH.pdf", width=20, height=8)
DotPlot(sce, features = genes,
        dot.scale = 8,cols  =c("white","#284aaa")
)+ RotatedAxis()
dev.off()
##############################################################################
##############################################################################
##############################################################################
#FOR human
sce=seurat_object_human
sce <- SCTransform(sce, method = "glmGamPoi")
Idents(sce) <- "cell_type"
table(Idents(sce))

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

noisy.liz <- NoisyGenes(sce, 0.5, levels(Idents(sce)))#0.65
neurons_avg_data_filtered <- neurons_avg_data[!rownames(neurons_avg_data) %in% noisy.liz,]
neurons_avg_TFs <- neurons_avg_data_filtered[rownames(neurons_avg_data_filtered) %in% Human_TFs$HGNC.symbol,]
tfs <- as.data.frame(neurons_avg_TFs)

# First, make sure your cell types are properly set
Idents(sce) <- "cell_type"
table(Idents(sce))

# Set the default assay
DefaultAssay(sce) <- "SCT"

# Get average expression - make sure to return in correct format
neurons_avg <- AverageExpression(sce, assays = "SCT", group.by = "cell_type", return.seurat = FALSE)

# Extract the SCT matrix
neurons_avg_data <- neurons_avg[["SCT"]]
write.csv(neurons_avg_data,file="human_average_data.csv")
# Your NoisyGenes function remains the same
noisy.liz <- NoisyGenes(sce, 0.2, levels(Idents(sce)))

# Filter the data
neurons_avg_data_filtered <- neurons_avg_data[!rownames(neurons_avg_data) %in% noisy.liz,]

# Filter for TFs
neurons_avg_TFs <- neurons_avg_data_filtered[rownames(neurons_avg_data_filtered) %in% Human_TFs$HGNC.symbol,]

# Convert to data frame and ensure proper cell type names
tfs <- as.data.frame(neurons_avg_TFs)

# If column names are still numbers, you can explicitly set them:
if(all(colnames(tfs) == as.character(1:ncol(tfs)))) {
  colnames(tfs) <- levels(Idents(sce))
}
write.csv(tfs,file="human_tfs_for_hierarchical_clustering.csv")
############################################ 491 TFs remain##################################################
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

tfs$gene <- rownames(tfs)
DefaultAssay(sce) <- "SCT"

DotPlot(sce, features = tfs$gene,cols = c("white","#284aaa")
)+theme(axis.text.x = element_text(angle = 90, hjust = 1))
pdf(file="491_TFS_human.pdf", width=125, height=12)
DotPlot(sce, features = tfs$gene,
        dot.scale = 8,cols  =c("white","#284aaa")
) + RotatedAxis()
dev.off(
  
)
genes <- c("ARHGAP35","ESR2","GATA3","ZNF10","ISL1","ZHX2","IKZF2","SMAD9",
           "MXD1","NACC2","NR6A1","PBX1", "PBX3",  "RFX3","TBX2","ARNTL",
           "RORB","THRB","ZFP64","NFATC3","NFATC2",
           "DACH1","ESRRB","ESRRG","AKAP8L","CASZ1","CUX2","EMX2","LMX1A",
           "MEIS1","MEIS2","AEBP2","BHLHE41","BNC2","ETS1","ETV5","FLYWCH1",
           "IRF4","ARID3A","ATF6","CREM","FOSL2","KLF3","NFE2",
           "EBF2","HEYL","EPAS1","ETS2","CUX1","POU3F1","ZEB2",
           "AHR","CDC5L","CREB5","FOXP2","MEF2A","MEF2C","RUNX2",
           "CAMTA1","MYT1L","PROX1","TSHZ3","ZNF804A",
           "ARNT2","KLF4","MAF","POU3F3","TFDP1","ERG","GLI2","NR2F1",
           "TBX18"
           
)

DotPlot(sce, features = genes,cols = c("white","#284aaa"))+theme(axis.text.x = element_text(angle = 90, hjust = 1))
pdf(file="HUMAN_TFS_COCH.pdf", width=20, height=8)
DotPlot(sce, features = genes,
        dot.scale = 8,cols  =c("white","#284aaa")
)+ RotatedAxis()
dev.off()
#################################################################################
#################################################################################
#######################################################################################
#######################################################################################
#load cross-species integration dataset
adata_loom <- connect(filename = "cross_species_integration_final_v2.loom",
                      mode = "r+",skip.validate = TRUE)
matrix=adata_loom[["matrix"]][,]
matrix=t(matrix)
dim(matrix)
gene = adata_loom$row.attrs$var_names[]
barcode = adata_loom$col.attrs$obs_names[]

meta_data = read.csv('cross_species_integration_final_obs.csv',row.names = 1) # as form as dataframe format
meta_feature = read.csv('cross_species_integration_final_var.csv',row.names = 1)

colnames(matrix)= barcode
row.names(matrix)= gene
#x_scvi = adata_loom$col.attrs$X_scVI[,]
x_umap=read.csv("cross_species_integration_final_umap.csv",row.names = 1)


seurat_object_integration= CreateSeuratObject(counts = matrix,meta.data = meta_data,
                                        project = 'cochlea_loom',
                                        min.cells = 0, 
                                        min.features = 0)
seurat_object_integration@assays[["RNA"]]@meta.features <- meta_feature

rownames(x_umap) = barcode
colnames(x_umap) = c('UMAP_1','UMAP_2')


seurat_object_integration$cellid <- rownames(seurat_object_integration@meta.data)
head(seurat_object_integration$cellid)
rownames(x_umap)=seurat_object_integration$cellid


seurat_object_integration <-FindVariableFeatures(seurat_object_integration)
seurat_object_integration <- ScaleData(seurat_object_integration)
seurat_object_integration <- RunPCA(seurat_object_integration)
seurat_object_integration<- RunUMAP(seurat_object_integration,dims = 1:30)
seurat_object_integration <- FindNeighbors(seurat_object_integration,dims = 1:30)
seurat_object_integration <- FindClusters(seurat_object_integration,resolution = 0.8)
DimPlot(seurat_object_integration,reduction = "umap",group.by = "cell_type")

seurat_object_integration@reductions[["umap"]]@cell.embeddings=as.matrix(x_umap)
DimPlot(seurat_object_integration,reduction = "umap",group.by = "cell_type_final")

sce_integration = as.SingleCellExperiment(seurat_object_integration)
colnames(x_umap) = c('umap_1','umap_2')
sce_integration@int_colData@listData[["reducedDims"]]@listData[["umap"]] = x_umap
#sce1@int_colData@listData[["reducedDims"]]@listData[["scANVI"]] = x_scANVI
seurat_object_integration <- as.Seurat(sce_integration)
#seurat_object <- RunUMAP(seurat_object,reduction = "scANVI",dims = 1:20)
DimPlot(seurat_object_integration,reduction = "umap",group.by = "cell_type_final")
DimPlot(seurat_object_integration,reduction = "umap",group.by = "region")
dev.off()
saveRDS(seurat_object_human,file="seurat_object_human.rds")
saveRDS(seurat_object_mouse,file="seurat_object_mouse.rds")
saveRDS(seurat_object_integration,file="seurat_object_integration.rds")
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
sce=seurat_object_integration
table(sce$leiden1.2)
Idents(sce) <- "leiden1.2"
DimPlot(sce, reduction = "umap",label = TRUE,repel = TRUE,group.by = "leiden1.2")

# 定义 leiden1.2 到 cell_type 的精确映射
cluster_mapping <- c(
  "0" = "Osteocytes/Osteoblasts",
  "1" = "Intermediate_stria",
  "2" = "Macrophages",
  "3" = "GCs",
  "4" = "Fibrocytes",
  "5" = "Claudius/Inner-Outer_sulcus_cells",
  "6" = "Interdental_cells",
  "7" = "Reissner_membrane",
  "8" = "Deiters_cells/Pillar_cells",
  "9" = "Endothelial_cellss",  
  "10" = "Inner_border-phalangeal/Hensen_cell",
  "11" = "Tympanic_border_cells",
  "12" = "Marginal_stria",
  "13" = "Basal_stria",
  "14" = "Deiters_cells/Pillar_cells", 
  "15" = "SS",
  "16" = "SGNs",
  "17" = "SGNs",
  "18" = "Root_cells",
  "19" = "HCs",
  "20" = "Inner_border-phalangeal/Hensen_cell",
  "21" = "Marginal_stria",
  "22" = "Tympanic_border_cells",
  "23" = "Spindle_cells",
  "24" = "Basal_stria",
  "25" = "HCs",
  "26" = "Pericytes",
  "27" = "SS",
  "28" = "SS",
  "29" = "SGNs"
)

# 将映射应用到 sce 对象
sce$integrated_cell_type <- cluster_mapping[as.character(sce$leiden1.2)]
Idents(sce) <- sce$integrated_cell_type  # 设置 cell_type 为新的 Ident

# 验证映射
table(Idents(sce), sce$leiden1.2)
DimPlot(sce, reduction = "umap",label = TRUE,repel = TRUE)
#remove cell_type_final 'Pericytes','Root_cells','Spindle_cells'
sce <- subset(sce, 
           subset = cell_type_final %in% c('Pericytes', 'Root_cells', 'Spindle_cells',
                                           "Endothelial_cells"),
              invert = TRUE)

sce <- SCTransform(sce, 
                                   vars.to.regress = c("nCount_RNA","pct_counts_mt"),
                                   verbose = TRUE, method="glmGamPoi")
DefaultAssay(sce) <- "SCT"
Idents(sce) <- "cell_type_final"
table(sce$cell_type_final)
brain_integrated_markers_conc <- matrix(,0,16)
table(sce$leiden1.2)
cluster_list <- c("Basal_stria","Claudius/Inner-Outer_sulcus_cells",
                  "Deiters_cells/Pillar_cells",
                  "Fibrocytes","GCs","HCs",
                  "Inner_border-phalangeal/Hensen_cell",
                  "Interdental_cells","Intermediate_stria","Macrophages",
                  "Marginal_stria","Osteocytes/Osteoblasts","Reissner_membrane",
                  "SGNs","Tympanic_border_cells",'SS'
  
  
                  )# 


Idents(sce) <- "integrated_cell_type"
brain_integrated_markers_conc <- matrix(,0,16)
cluster_list <- c("Osteocytes/Osteoblasts","Intermediate_stria","Macrophages",
                  "GCs","Fibrocytes","Claudius/Inner-Outer_sulcus_cells",
                  "Interdental_cells","Reissner_membrane","Deiters_cells/Pillar_cells",
                  "Inner_border-phalangeal/Hensen_cell",
                  "Marginal_stria","Basal_stria",
                 'SS',"SGNs",
                  "HCs",
                  "Tympanic_border_cells"
                
                  
                  
)# 
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
                organism = "hsapiens",
                ordered_query = FALSE,
                multi_query = FALSE, 
                significant = TRUE, 
                exclude_iea = FALSE,
                measure_underrepresentation = FALSE, 
                evcodes = FALSE,
                user_threshold = 0.01, 
                correction_method = "g_SCS",
                domain_scope = "annotated",
                custom_bg = NULL,
                numeric_ns = "",  
                as_short_link = FALSE,
                sources = c("GO:MF", "GO:BP","GO:CC"))


tf_gene_list <- read.csv("human_tf_gene_list.csv",row.names = 1)
brain_integrated_markers_conc$HGNC.symbol <- brain_integrated_markers_conc$gene
conserved_tf <- merge(x=brain_integrated_markers_conc,y=tf_gene_list, by="HGNC.symbol", all=FALSE)

write.csv(conserved_tf,file = "conserved_tf_across_human_and_mouse_new.csv")

####################################################################################################
#METANEIGHBOR ANALYSIS
cluster_mapping <- c(
  "0" = "Lateral Wall & SS",
  "1" = "Lateral Wall & SS",
  "2" = "Circulating & Others",
  "3" = "Modiolus",
  "4" = "Lateral Wall & SS",
  "5" = "Cochlear_epithelium",
  "6" = "Cochlear_epithelium",
  "7" = "Lateral Wall & SS",
  "8" = "Cochlear_epithelium",
  "9" = "Lateral Wall & SS",  
  "10" = "Cochlear_epithelium",
  "11" = "Lateral Wall & SS",
  "12" = "Lateral Wall & SS",
  "13" = "Lateral Wall & SS",
  "14" = "Cochlear_epithelium", 
  "15" = "Lateral Wall & SS",
  "16" = "Modiolus",
  "17" = "Modiolus",
  "18" = "Lateral Wall & SS",
  "19" = "Cochlear_epithelium",
  "20" = "Cochlear_epithelium",
  "21" = "Lateral Wall & SS",
  "22" = "Lateral Wall & SS",
  "23" = "Lateral Wall & SS",
  "24" = "Lateral Wall & SS",
  "25" = "Cochlear_epithelium",
  "26" = "Lateral Wall & SS",
  "27" = "Lateral Wall & SS",
  "28" = "Lateral Wall & SS",
  "29" = "Modiolus"
)

# 将映射应用到 sce 对象
sce$class_label <- cluster_mapping[as.character(sce$leiden1.2)]
Idents(sce) <- sce$class_label  # 设置 cell_type 为新的 Ident

# 验证映射
table(Idents(sce), sce$leiden1.2)
DimPlot(sce, reduction = "umap",label = TRUE,repel = TRUE)


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



################## Prep SingleCellExperiment object for analysis  #############################################
library(SingleCellExperiment)
library(Matrix)
#sce$class_label = "CoE"
table(sce$class_label)
Idents(sce) <- "class_label"
sce$class.species <- paste(sce$class_label, sce$species, sep = "_")
table(sce$class.species)
Idents(sce) <- "species"
table(sce$species)
hs <- subset(sce, idents ="human")
ms <- subset(sce, idents ="mouse")

hs <- GetAssayData(object = hs[["RNA"]], slot = "counts")
ms <- GetAssayData(object = ms[["RNA"]], slot = "counts")
gs=read.csv("orthologTable_human_macaque_mouse.csv",header=T,stringsAsFactors=F)  ### read in ortholog mapping table


hs <- as.matrix(hs)
ms <- as.matrix(ms)

p2 <- subset(sce, idents ="human")
p3<- subset(sce, idents ="mouse")

p2$sample_id <- colnames(p2)
p3$sample_id <- colnames(p3)

p2$study_id="human"
p3$study_id="mouse"
p2<- as.matrix(p2@meta.data)
p3<- as.matrix(p3@meta.data)

m<-match(as.character(rownames(hs)),as.character(gs$human_symbol))  ### match gene IDs from expression data to those in the ortholog table
f.a=!is.na(m)
f.b=m[f.a]
hs=hs[f.a,]
rownames(hs)=gs[f.b,"human_symbol"]  ## convert to human IDs (14.9K genes)
m<-match(as.character(rownames(ms)),as.character(gs$human_symbol))
f.a=!is.na(m)
f.b=m[f.a]
ms=ms[f.a,]
rownames(ms)=as.character(gs[f.b,"human_symbol"]) ## 14350 genes
rownames(p3)=colnames(ms)
m<-match(as.character(rownames(hs)),as.character(rownames(ms)))
f.a=!is.na(m)
f.b=m[f.a]

dat=cbind(hs[f.a,],ms[f.b,])
x=as.vector(rownames(dat))
rownames(dat)=NULL
rownames(dat)=x
p_pri2=rbind(p2,p3)
rownames(p_pri2)=colnames(dat)
sce_all=SingleCellExperiment(assays=list(counts=dat),colData=p_pri2)  
sce_all

saveRDS(sce_all,file="hs_ms_sce_V3.rds")

sce_all$sample_id_append <- sce_all$sample_id
head(sce_all$sample_id_append)
sce$sample_id <-colnames(sce)
head(sce$sample_id)
m<-match(sce_all$sample_id_append,sce$sample_id)
sum(!is.na(m))
f.a=!is.na(m)
f.b=m[f.a]

sce_all$cell_type_final[f.a]=as.character(sce$cell_type_final[f.b])
#sce_all$final_integrated_cluster_color[f.a]=cochlea_epithelium$[f.b]
head(sce_all$study_id)
head(sce_all$cell_type_final[f.a])
sce_all=sce_all[,f.a]
head(sce_all$study_id)
saveRDS(sce_all,file="hs_ms_sce_v3.rds")
sce_all

########## Within- and cross-species classification with highly variable genes ###################

library(SingleCellExperiment)
library(Matrix)
source("1v1_analysis.R")
classes=unique(sce_all$class_label)
#classes <- classes[-1]
vgs=vector("list",length=length(classes))
names(vgs)=names(classes)
for(i in seq_along(classes)){
  f=sce_all$class_label==classes[i]
  vgs[[i]]=get_variable_genes(sce_all[,f.b])
}

mn_1v1_subclass=vector("list",length=length(classes))
for(i in seq_along(classes)){
  f=sce_all$class_label==classes[i]
  
  mn_1v1_subclass[[i]]=compute_best_hits(sce_all[vgs[[i]],f.b],
                                         sce_all$cell_type_final[f.b],one_vs_one=TRUE)
}
saveRDS(mn_1v1_subclass,file="coch_1v1_celltype_level_v3.Rdata")

b<- as.data.frame(mn_1v1_subclass)
write.csv(b,file = "mn_1v1_coch_celltype_v3.csv")
options(stringsAsFactors = FALSE)

d<- read.csv("mn_1v1_coch_celltype_v3.csv",row.names = 1)
library(pheatmap)
library(viridis)
pheatmap(d,cluster_rows = FALSE,cluster_cols = FALSE,border=FALSE)
pheatmap::pheatmap(d,cluster_rows = FALSE,cluster_cols = FALSE,border=FALSE,color = viridis(50))
#color = colorRampPalette(c("navy", "white", "firebrick3"))(100)
pheatmap(d,cluster_rows = FALSE,cluster_cols = FALSE,border=FALSE,color = colorRampPalette(c("navy", "white", "firebrick3"))(100))
pheatmap::pheatmap(d,cluster_rows = FALSE,cluster_cols = FALSE,border=FALSE,color = cividis(50))
########## Within- and cross-species classification using HGNC and SynGO gene sets ###################
library(SingleCellExperiment)
library(Matrix)
options(stringsAsFactors = FALSE)
source("metaneighbor.R")
gs=readRDS("hgnc_syngo.rds")
################ within mouse
head(sce_all$study_id)
spe_f=sce_all$study_id=="mouse"
sce_ms=sce_all[,spe_f]

rand3cv=sample(c(1:3),size=sum(spe_f),replace=T)
sce_ms$study_id=rand3cv
sce_ms$class_label='cochlea'
classes=unique(sce_ms$class_label)
res=vector("list",length=length(classes))
names(res)=classes
head(sce_ms$cell_type_final)
for(i in seq_along(classes)){
  sub_f=sce_ms$class_label==classes[i]
  x=model.matrix(~(as.character(sce_ms$cell_type_final[sub_f]))+0)
  colnames(x)=names(table(sce_ms$cell_type_final[sub_f]))
  res[[i]]=call_my_metaneighbor(sce_ms[,sub_f],gs,x)
}

saveRDS(res,file="within_mouse_MN_cluster_level_v3.rds")
res_mouse=res
rm(sce_ms)
################ within human
spe_f=sce_all$study_id=="human"
sce_hs=sce_all[,spe_f]
rand3cv=sample(c(1:3),size=sum(spe_f),replace=T)
sce_hs$study_id=rand3cv
head(sce_hs$study_id)
sce_hs$class_label='cochlea'
classes=unique(sce_hs$class_label)
res=vector("list",length=length(classes))
names(res)=classes
head(sce_hs$cell_type_final,n=20)
for(i in seq_along(classes)){
  sub_f=sce_hs$class_label==classes[i]
  x=model.matrix(~(as.character(sce_hs$cell_type_final[sub_f]))+0)
  colnames(x)=names(table(sce_hs$cell_type_final[sub_f]))
  res[[i]]=call_my_metaneighbor(sce_hs[,sub_f],gs,x)
}
saveRDS(res,file="within_human_MN_cluster_level_v3.rds")
res_human=res
rm(sce_mac)
################ cross macaque and mouse
sce_all2=sce_all
sce_all2$study_original=sce_all2$study_id
sce_all2$class_label="cochlea"
classes=unique(sce_all2$class_label)
res=vector("list",length=length(classes))
names(res)=classes

for(i in seq_along(classes)){
  sub_f=sce_all2$class_label==classes[i]
  x=model.matrix(~(as.character(sce_all2$cell_type_final[sub_f]))+0)
  colnames(x)=names(table(sce_all2$cell_type_final[sub_f]))
  res[[i]]=call_my_metaneighbor(sce_all2[,sub_f],gs,x)
}


save(res,file="cross_human_mouse_cluster_level2_v3.rds")
res_human_mouse=res


write.csv(res_human,file = "1within_human_mn_cluster_1.csv")
write.csv(res_mouse,file = "1within_mouse_mn_cluster_level_v3.csv")
write.csv(res_human_mouse,file = "1cross_human_mouse_mn_cluster_level_v3.csv")



b <- read.csv("1within_human_mn_cluster_1.csv")
b
b<- aggregate(b$cochlea.auroc,by=list(type=b$cochlea.gene_set),mean)
c <- read.csv("1within_mouse_mn_cluster_level_v3.csv")
c
c<- aggregate(c$cochlea.auroc,by=list(type=c$cochlea.gene_set),mean)

f <- read.csv("1cross_human_mouse_mn_cluster_level_v3.csv")
f
f<- aggregate(f$cochlea.auroc,by=list(type=f$cochlea.gene_set),mean)

g <- merge(x=b,y=c,by="type")
g <- merge(x=g,y=f,by="type")
write.csv(g, file="MN_CLUSTER_meanROC_coe_v3.csv")
g<-read.csv("MN_CLUSTER_meanROC_coe_v3.csv")
# Library
library(ggplot2)
library(hrbrthemes)
p3 <- ggplot(g, aes(g$within_species_meanROC,g$human_meanROC)) +
  geom_point(color="#69b3a2") +
  geom_smooth(method=lm , color="black", se=FALSE) +
  theme_ipsum()+
  xlim(0.4,1)+
  ylim(0.4,1)
p3
p4 <- ggplot(g, aes(g$within_species_meanROC,g$mouse_meanROC)) +
  geom_point(color="#69b3a2") +
  geom_smooth(method=lm , color="black", se=FALSE) +
  theme_ipsum()+
  xlim(0.4,1)+
  ylim(0.4,1)
p4
p7 <- ggplot(g, aes(g$within_species_meanROC,g$human_mouse_meanROC)) +
  geom_point(color="#69b3a2") +
  geom_smooth(method=lm , color="black", se=FALSE) +
  theme_ipsum()+xlim(0.4,1)+
  ylim(0.4,1)
p7
coef(lm(g$human_mouse_meanROC~g$within_species_meanROC))[2]
g<-read.csv("MN_CLUSTER_meanROC_coe_v3.csv")
h<-read.csv("Supplementary Table 9_2.csv")
i<- merge(x=g,y=h,by="gene_set_label")
j<- i[!duplicated(i$gene_set_label), ]
k<- merge(x=j,y=g, all = TRUE)
write.csv(k,file = "coch_mn_cluster_meanROC_v3.csv")

library(tidyverse)
library(tidyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
mn.df <- read.csv("coch_mn_cluster_meanROC_v3.csv")
mn.df$gene_set_type <- factor(mn.df$gene_set_type, 
                              levels = c("Other", "Signaling", 
                                         "Cell Adhesion", "Ion Channel"))
roc.tests <- c("human_meanROC", "mouse_meanROC", 
               "human_mouse_meanROC")
mn.l <- mn.df %>% 
  subset(cell_class == "cochlea") %>%
  arrange(gene_set_type) %>% 
  gather(roc.tests, key = "comp", value = "ROC")
mn.l$comp <- factor(mn.l$comp, levels = roc.tests)
levels(mn.l$comp) <- c("Human", "Mouse", 
                       
                       "Human vs. Mouse")
mn.l$comp_type <- ifelse(mn.l$comp %in% c("Human", "Mouse"), 
                         "Within-species", "Cross-species")
geneset.pal <- c("#f70ff1", "#2626d3", "#bbce32", "#4ce03d","#998a82")
g.scatter <- ggplot(mn.l, aes(x = within_species_meanROC, y = ROC)) +
  facet_grid(. ~ comp) +
  geom_hline(yintercept = 0.5, color = "grey", size = 0.25) +
  geom_vline(xintercept = 0.5, color = "grey", size = 0.25) +
  geom_abline(slope = 1, intercept = 0, color = "grey", size = 0.25) +
  geom_point(alpha = 0.5, aes(color = gene_set_type)) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  scale_colour_manual(values=geneset.pal) +
  coord_fixed(xlim = c(0.3, 1), ylim = c(0.3, 1)) +
  scale_x_continuous(breaks = c(0.3, 0.5,0.7, 0.9)) +
  scale_y_continuous(breaks = c(0.3, 0.5,0.7, 0.9)) +
  xlab("Within-species mean AUROC") +
  ylab("AUROC") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

plot(g.scatter)
mn.l %>% 
  group_by(comp) %>%
  group_modify(~ broom::tidy(lm(ROC ~ within_species_meanROC, data = .x)))
