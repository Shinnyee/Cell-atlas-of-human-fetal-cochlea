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
adata_loom <- connect(filename = "adata_mouse_v2_OHC.loom",
                      mode = "r+",skip.validate = TRUE)#chang file for respective dataset
matrix=adata_loom[["matrix"]][,]
matrix=t(matrix)
dim(matrix)
gene = adata_loom$row.attrs$var_names[]
barcode = adata_loom$col.attrs$obs_names[]

meta_data = read.csv('adata_mouse_v2_OHC_obs.csv',row.names = 1) # chang file for respective dataset
meta_feature = read.csv('adata_mouse_v2_OHC_var.csv',row.names = 1)#chang file for respective dataset

colnames(matrix)= barcode
row.names(matrix)= gene


seurat_object= CreateSeuratObject(counts = matrix,meta.data = meta_data,
                                  project = 'mouse_OHC_loom',
                                  min.cells = 0, 
                                  min.features = 0)
seurat_object@assays[["RNA"]]@meta.features <- meta_feature
seurat_object <-FindVariableFeatures(seurat_object)
seurat_object <- ScaleData(seurat_object)
seurat_object <- RunPCA(seurat_object)
seurat_object<- RunUMAP(seurat_object,dims = 1:30)
seurat_object <- FindNeighbors(seurat_object,dims = 1:30)
seurat_object <- FindClusters(seurat_object,resolution = 0.8)
DimPlot(seurat_object,reduction = "umap",group.by = "age_bins2")

mouse_OHC=seurat_object
mouse_IHC=seurat_object
human_IHC=seurat_object
human_OHC=seurat_object
human_IHC2=human_IHC
human_OHC2=human_OHC
mouse_OHC2=mouse_OHC
mouse_IHC2=mouse_IHC
################################################################################
################################################################################
DefaultAssay(human_IHC2) <- "RNA"
DefaultAssay(human_OHC2) <- "RNA"
DefaultAssay(mouse_OHC2) <- "RNA"
DefaultAssay(mouse_IHC2) <- "RNA"
human_IHC2 <- human_IHC2 %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData()
human_OHC2 <- human_OHC2 %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData()
mouse_OHC2 <- mouse_OHC2 %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData()
mouse_IHC2 <- mouse_IHC2 %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData()
################################################################################
Idents(mouse_IHC2) <- "age_bins2"
Idents(mouse_OHC2) <- "age_bins2"

Idents(human_OHC2) <- "age_bins2"
Idents(human_IHC2) <- "age_bins2"
table(human_IHC2$age_bins2)
deg_hu_ihc <- FindAllMarkers(human_IHC2, assay = "RNA", slot = "data",
                               test.use = "roc")
deg_mu_ihc <- FindAllMarkers(mouse_IHC2, assay = "RNA", slot = "data",
                             test.use = "roc")
deg_hu_ihc <- deg_hu_ihc[which(deg_hu_ihc$avg_diff > 0), ]
deg_mu_ihc <- deg_mu_ihc[which(deg_mu_ihc$avg_diff > 0), ]

##################################################################################
library(tidyverse)
library(readxl)
library(data.table)
library (ggplot2)
#install.packages("ggpubr")
library(ggpubr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(biomaRt)

cluster1=names(table(deg_hu_ihc$cluster))
cluster2=names(table(deg_mu_ihc$cluster))
## jaccard
df_ja=c()

for (i in cluster1) {
  ja=c()
  for (j in cluster2) {
    a=deg_hu_ihc[deg_hu_ihc$cluster==i,]
    a=a$gene
    b=deg_mu_ihc[deg_mu_ihc$cluster==j,]
    b=b$gene
    jaccard=length(intersect(a,b))/length(union(a,b))
    
    ja=c(ja,jaccard)
  }
  df_ja=rbind(df_ja,ja)
}

rownames(df_ja)=paste0(cluster1)
colnames(df_ja)=paste0(cluster2)
library(pheatmap)

pheatmap(df_ja,cluster_cols = F,
         cluster_rows = F,show_rownames = T, show_colnames =T,
         color = colorRampPalette(c(rep("blue",1), "white", rep("red",1)))(100),
         breaks = seq(0,0.25,length.out = 100),border_color = NA)
library(RColorBrewer)
## Set the common color scale
breaksList = seq(0, 0.2, by = 0.00001)


heatmap_up_noleg = pheatmap(mat = df_ja, 
                            cluster_rows = T,
                            cluster_cols = T,
                            display_numbers = F,
                            number_format = "%.1f",
                            scale = "none",
                            na_col = "#DDDDDD", 
                            border_color = NA, 
                            color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), 
                            breaks = breaksList, legend = FALSE)


heatmap_up_noleg = pheatmap(mat = df_ja, 
                            cluster_rows = F,
                            cluster_cols = F,
                            display_numbers = F,
                            number_format = "%.1f",
                            scale = "none",
                            na_col = "#DDDDDD", 
                            border_color = NA, 
                            color = colorRampPalette(rev(brewer.pal(n = 7,
                                                                    name = "RdYlBu")))(length(breaksList)), 
                            breaks = breaksList, legend = FALSE)
##########FOR OUTER HAIR CELLS
################################################################################

deg_hu_ohc <- FindAllMarkers(human_OHC2, assay = "RNA", slot = "data",
                             test.use = "roc")
deg_mu_ohc <- FindAllMarkers(mouse_OHC2, assay = "RNA", slot = "data",
                             test.use = "roc")
deg_hu_ohc <- deg_hu_ohc[which(deg_hu_ohc$avg_diff > 0), ]
deg_mu_ohc <- deg_mu_ohc[which(deg_mu_ohc$avg_diff > 0), ]

##################################################################################

cluster1=names(table(deg_hu_ohc$cluster))
cluster2=names(table(deg_mu_ohc$cluster))
## jaccard
df_ja=c()

for (i in cluster1) {
  ja=c()
  for (j in cluster2) {
    a=deg_hu_ohc[deg_hu_ohc$cluster==i,]
    a=a$gene
    b=deg_mu_ohc[deg_mu_ohc$cluster==j,]
    b=b$gene
    jaccard=length(intersect(a,b))/length(union(a,b))
    
    ja=c(ja,jaccard)
  }
  df_ja=rbind(df_ja,ja)
}

rownames(df_ja)=paste0(cluster1)
colnames(df_ja)=paste0(cluster2)
library(pheatmap)

pheatmap(df_ja,cluster_cols = F,
         cluster_rows = F,show_rownames = T, show_colnames =T,
         color = colorRampPalette(c(rep("blue",1), "white", rep("red",1)))(100),
         breaks = seq(0,0.25,length.out = 100),border_color = NA)
library(RColorBrewer)
## Set the common color scale
breaksList = seq(0, 0.2, by = 0.00001)


heatmap_up_noleg = pheatmap(mat = df_ja, 
                            cluster_rows = T,
                            cluster_cols = T,
                            display_numbers = F,
                            number_format = "%.1f",
                            scale = "none",
                            na_col = "#DDDDDD", 
                            border_color = NA, 
                            color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), 
                            breaks = breaksList, legend = FALSE)


heatmap_up_noleg = pheatmap(mat = df_ja, 
                            cluster_rows = F,
                            cluster_cols = F,
                            display_numbers = F,
                            number_format = "%.1f",
                            scale = "none",
                            na_col = "#DDDDDD", 
                            border_color = NA, 
                            color = colorRampPalette(rev(brewer.pal(n = 7,
                                                                    name = "RdYlBu")))(length(breaksList)), 
                            breaks = breaksList, legend = FALSE)
