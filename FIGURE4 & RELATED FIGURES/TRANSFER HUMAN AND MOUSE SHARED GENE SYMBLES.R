rm(list=ls())
#remotes::install_github("zh542370159/SCP")
#remotes::install_version("ggplot2", version = "3.5.0")
#devtools::install_github("zhanghao-njmu/SCP")
#library(SCP)
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
#remotes::install_version("ggplot2", version = "3.4.2")
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
adata_loom <- connect(filename = "human_adata_sgn.loom",
                      mode = "r+",skip.validate = TRUE)
matrix=adata_loom[["matrix"]][,]
matrix=t(matrix)
dim(matrix)
gene = adata_loom$row.attrs$var_names[]
barcode = adata_loom$col.attrs$obs_names[]

meta_data = read.csv('human_adata_sgnn_obs.csv',row.names = 1) # as form as dataframe format
meta_feature = read.csv('human_adata_sgn_var.csv',row.names = 1)

colnames(matrix)= barcode
row.names(matrix)= gene
#x_scvi = adata_loom$col.attrs$X_scVI[,]
#x_umap=read.csv("human_HC_annotation_umap.csv",row.names = 1)
x_scANVI=read.csv("human_adata_sgn_scANVI.csv",row.names = 1)

seurat_object= CreateSeuratObject(counts = matrix,meta.data = meta_data,
                                  project = 'human_sgn_loom',
                                  min.cells = 0, 
                                  min.features = 0)
seurat_object@assays[["RNA"]]@meta.features <- meta_feature


#rownames(x_umap) = barcode
#colnames(x_umap) = c('UMAP_1','UMAP_2')
rownames(x_scANVI) = barcode
colnames(x_scANVI) = c('scANVI_1','scANVI_2','scANVI_3','scANVI_4','scANVI_5','scANVI_6',
                       'scANVI_7','scANVI_8','scANVI_9','scANVI_10','scANVI_11','scANVI_12',
                       'scANVI_13','scANVI_14','scANVI_15','scANVI_16','scANVI_17','scANVI_18',
                       'scANVI_19','scANVI_20','scANVI_21','scANVI_22','scANVI_23','scANVI_24',
                       'scANVI_25','scANVI_26','scANVI_27','scANVI_28','scANVI_29','scANVI_30')

seurat_object$cellid <- rownames(seurat_object@meta.data)
head(seurat_object$cellid)
#rownames(x_umap)=seurat_object$cellid
rownames(x_scANVI)=seurat_object$cellid
seurat_object <-FindVariableFeatures(seurat_object)
seurat_object <- ScaleData(seurat_object)
seurat_object <- RunPCA(seurat_object)
seurat_object<- RunUMAP(seurat_object,dims = 1:30)
seurat_object <- FindNeighbors(seurat_object,dims = 1:30)
seurat_object <- FindClusters(seurat_object,resolution = 0.8)
DimPlot(seurat_object,reduction = "umap",group.by = "gw")

#seurat_object@reductions[["umap"]]@cell.embeddings=as.matrix(x_umap)
DimPlot(seurat_object,group.by = "gw")
saveRDS(seurat_object,file="human_sgn.rds")
human_sgn <- readRDS("F:/PROJECTS/PROJECT_HUMAN_FETAL_COCHLEAE/WORKPLACE/R/sc_pipeline/SGN analysis/human_sgn.rds")
#human_hc <- readRDS("F:/PROJECTS/PROJECT_HUMAN_FETAL_COCHLEAE/WORKPLACE/R/sc_pipeline/HC_DEGs_analysis_GSEA_analysis/human_hc.rds")
library(tidyverse)
library(DT)
library(biomaRt)
# 读取文件，MOUSE ID TO HUMAN ID
#genesV2 <- read.table("mouse_to_human_genes.txt", sep="\t", header=T)
#genesV2 <- read.csv("orthologTable_human_macaque_mouse.csv",row.names = 1)
genesV2 <- read.csv("mouse_to_human_genes.csv")
#write.csv(genesV2,file = "mouse_to_human_genes.csv")
## Extract Expression Data
human_sgn_counts <- as.matrix(GetAssayData(human_sgn, slot = "counts"))
human_sgn_counts <- data.frame(gene=rownames(human_sgn_counts), human_sgn_counts, check.names = F)
dim(human_sgn_counts)

# 转人类基因名为人类基因名
human_sgn_counts$Gene <- genesV2[match(human_sgn_counts$gene, genesV2[,2]),2]
human_sgn_counts <- subset(human_sgn_counts, Gene!='NA')
human_sgn_counts <- dplyr::select(human_sgn_counts, Gene, everything())
human_sgn_counts <- human_sgn_counts[, !(colnames(human_sgn_counts) %in% 'gene')]
dim(human_sgn_counts)
rownames(human_sgn_counts)<-human_sgn_counts[,1]
human_sgn_counts<-human_sgn_counts[,-1]
head(rownames(human_sgn_counts))
meta_data=human_sgn_counts@meta.data
human_sgn_counts_v2 <- CreateSeuratObject(counts = human_sgn_counts,
                                         meta.data = meta_data,
                                      project = 'human_sgn_loom',
                                      min.cells = 0, 
                                      min.features = 0)

#rownames(x_umap) = barcode
#colnames(x_umap) = c('UMAP_1','UMAP_2')
rownames(x_scANVI) = barcode
colnames(x_scANVI) = c('scANVI_1','scANVI_2','scANVI_3','scANVI_4','scANVI_5','scANVI_6',
                       'scANVI_7','scANVI_8','scANVI_9','scANVI_10','scANVI_11','scANVI_12',
                       'scANVI_13','scANVI_14','scANVI_15','scANVI_16','scANVI_17','scANVI_18',
                       'scANVI_19','scANVI_20','scANVI_21','scANVI_22','scANVI_23','scANVI_24',
                       'scANVI_25','scANVI_26','scANVI_27','scANVI_28','scANVI_29','scANVI_30')

human_sgn_counts_v2$cellid <- rownames(human_sgn_counts_v2@meta.data)
head(human_sgn_counts_v2$cellid)
#rownames(x_umap)=human_hc_counts_v2$cellid
rownames(x_scANVI)=human_sgn_counts_v2$cellid

human_sgn_counts_v2 <-FindVariableFeatures(human_sgn_counts_v2)
human_sgn_counts_v2 <- ScaleData(human_sgn_counts_v2)
human_sgn_counts_v2 <- RunPCA(human_sgn_counts_v2)
human_sgn_counts_v2<- RunUMAP(human_sgn_counts_v2,dims = 1:30)
human_sgn_counts_v2 <- FindNeighbors(human_sgn_counts_v2,dims = 1:30)
human_sgn_counts_v2 <- FindClusters(human_sgn_counts_v2,resolution = 0.8)
DimPlot(human_sgn_counts_v2,reduction = "umap",group.by = "gw")

#human_sgn_counts_v2@reductions[["umap"]]@cell.embeddings=as.matrix(x_umap)
DimPlot(human_sgn_counts_v2,reduction = "umap",group.by = "gw")
saveRDS(human_sgn_counts_v2,file="human_sgn_filtered.rds")

FeaturePlot(human_sgn_counts_v2,features = c("DLX3","EPHA5","TUBB3","PRPH",
                                             "EPHA4","MAP2"
),pt.size = 0.8)
human_sgn_filtered <- readRDS("F:/PROJECTS/PROJECT_HUMAN_FETAL_COCHLEAE/WORKPLACE/R/sc_pipeline/SGN analysis/human_sgn_filtered.rds")

human_sgn_counts_v2=human_sgn_filtered
Idents(human_sgn_counts_v2) <- "seurat_clusters"
FeaturePlot(human_sgn_counts_v2,features = c("DLX3","EPHA5","TUBB3","PRPH",
                                             "EPHA4","MAP2"
),label = "seurat_clusters")

DimPlot(human_sgn_counts_v2,reduction = "umap",group.by = "seurat_clusters")
#remove cluster 0,5,9
# 加载必要的库
library(Seurat)
library(dplyr)

# 假设你的 Seurat 对象是 seurat_obj
# 假设你想去除的 cluster 是 cluster_to_remove

# 查看当前的 cluster 分布
table(human_sgn_counts_v2@meta.data$seurat_clusters)
cluster_to_remove=c("0","5","9")
# 选择要保留的细胞
cells_to_keep <- WhichCells(human_sgn_counts_v2, 
                            idents = setdiff(levels(human_sgn_counts_v2@meta.data$seurat_clusters), cluster_to_remove))

# 创建一个新的 Seurat 对象，只包含要保留的细胞
human_sgn_counts_v2_filtered <- subset(human_sgn_counts_v2, cells = cells_to_keep)

# 查看过滤后的 cluster 分布
table(human_sgn_counts_v2_filtered@meta.data$seurat_clusters)

FeaturePlot(human_sgn_counts_v2_filtered,features = c("DLX3","EPHA5","TUBB3","PRPH",
                                             "EPHA4","MAP2"
),pt.size = 0.9)
#harmony integration
library(harmony)
table(human_sgn_counts_v2_filtered$orig.ident)
sce_all=human_sgn_counts_v2_filtered
DefaultAssay(sce_all) <- "RNA"
sce_all <- sce_all %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData()
sce_all <- RunPCA(sce_all, features = VariableFeatures(sce_all), npcs = 50)

sce_all_harmony <- RunHarmony(sce_all,group.by.vars = 'orig.ident',reduction = "pca",
                              dims.use = 1:8,assay.use = "RNA")#lambda = 3,theta = 0.6
#sce_all_harmony <- RunHarmony(sce_all,group.by.vars = 'orig.ident',reduction = "pca",
#dims.use = 1:20,assay.use = "RNA")
sce_all[["harmony"]] <- sce_all_harmony[["harmony"]]
sce_all <- RunUMAP(sce_all,dims = 1:6,
                   reduction = "harmony",reduction.name = "umap_harmony")

p1 <- DimPlot(sce_all, reduction = "umap_harmony", group.by = "orig.ident",pt.size = 1) + 
  ggtitle("UMAP Harmony")
p1
p2 <- DimPlot(sce_all, reduction = "umap_harmony",pt.size = 1,
              group.by = "seurat_clusters",repel = TRUE) + 
  ggtitle("UMAP Harmony")
p2




