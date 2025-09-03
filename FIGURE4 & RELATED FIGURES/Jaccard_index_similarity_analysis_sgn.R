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
adata_loom <- connect(filename = "adata_mouse_sgn_jaccard_index.loom",
                      mode = "r+",skip.validate = TRUE)#chang file for respective dataset
matrix=adata_loom[["matrix"]][,]
matrix=t(matrix)
dim(matrix)
gene = adata_loom$row.attrs$var_names[]
barcode = adata_loom$col.attrs$obs_names[]

meta_data = read.csv('adata_mouse_sgn_jaccard_index_obs.csv',row.names = 1) # chang file for respective dataset
meta_feature = read.csv('adata_mouse_sgn_jaccard_index_var.csv',row.names = 1)#chang file for respective dataset

colnames(matrix)= barcode
row.names(matrix)= gene


seurat_object= CreateSeuratObject(counts = matrix,meta.data = meta_data,
                                  project = 'human_sgn_loom',
                                  min.cells = 0, 
                                  min.features = 0)
seurat_object@assays[["RNA"]]@meta.features <- meta_feature
seurat_object <-FindVariableFeatures(seurat_object)
seurat_object <- ScaleData(seurat_object)
seurat_object <- RunPCA(seurat_object)
seurat_object<- RunUMAP(seurat_object,dims = 1:30)
seurat_object <- FindNeighbors(seurat_object,dims = 1:30)
seurat_object <- FindClusters(seurat_object,resolution = 0.8)
DimPlot(seurat_object,reduction = "umap",group.by = "celltype")
human_sgn=seurat_object
mouse_sgn=seurat_object



################################################################################
################################################################################
DefaultAssay(human_sgn) <- "RNA"
DefaultAssay(mouse_sgn) <- "RNA"

human_sgn <- human_sgn %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData()
mouse_sgn <- mouse_sgn %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData()
################################################################################
Idents(human_sgn) <- "celltype"
Idents(mouse_sgn) <- "celltype"


deg_hu_sgn <- FindAllMarkers(human_sgn, assay = "RNA", slot = "data",
                             test.use = "roc")
deg_mu_sgn <- FindAllMarkers(mouse_sgn, assay = "RNA", slot = "data",
                             test.use = "roc")
deg_hu_sgn <- deg_hu_sgn[which(deg_hu_sgn$avg_diff > 0), ]
deg_mu_sgn <- deg_mu_sgn[which(deg_mu_sgn$avg_diff > 0), ]

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

cluster1=names(table(deg_hu_sgn$cluster))
cluster2=names(table(deg_mu_sgn$cluster))
## jaccard
df_ja=c()

for (i in cluster1) {
  ja=c()
  for (j in cluster2) {
    a=deg_hu_sgn[deg_hu_sgn$cluster==i,]
    a=a$gene
    b=deg_mu_sgn[deg_mu_sgn$cluster==j,]
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
