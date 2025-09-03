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
adata_loom <- connect(filename = "human_cochlea_annotation.loom",
                      mode = "r+",skip.validate = TRUE)
matrix=adata_loom[["matrix"]][,]
matrix=t(matrix)
dim(matrix)
gene = adata_loom$row.attrs$var_names[]
barcode = adata_loom$col.attrs$obs_names[]

meta_data = read.csv('human_cochlea_annotation_obs.csv',row.names = 1) # as form as dataframe format
meta_feature = read.csv('human_cochlea_annotation_var.csv',row.names = 1)

colnames(matrix)= barcode
row.names(matrix)= gene
#x_scvi = adata_loom$col.attrs$X_scVI[,]
x_umap=read.csv("human_cochlea_annotation_umap.csv",row.names = 1)
x_scANVI=read.csv("human_cochlea_annotation_scANVI.csv",row.names = 1)

seurat_object= CreateSeuratObject(counts = matrix,meta.data = meta_data,
                                      project = 'human_cochlea_loom',
                                      min.cells = 0, 
                                      min.features = 0)
seurat_object@assays[["RNA"]]@meta.features <- meta_feature


rownames(x_umap) = barcode
colnames(x_umap) = c('UMAP_1','UMAP_2')
rownames(x_scANVI) = barcode
colnames(x_scANVI) = c('scANVI_1','scANVI_2')

seurat_object$cellid <- rownames(seurat_object@meta.data)
head(seurat_object$cellid)
rownames(x_umap)=seurat_object$cellid
rownames(x_scANVI)=seurat_object$cellid

seurat_object <-FindVariableFeatures(seurat_object)
seurat_object <- ScaleData(seurat_object)
seurat_object <- RunPCA(seurat_object)
seurat_object<- RunUMAP(seurat_object,dims = 1:30)
seurat_object <- FindNeighbors(seurat_object,dims = 1:30)
seurat_object <- FindClusters(seurat_object,resolution = 0.8)
DimPlot(seurat_object,reduction = "umap",group.by = "gw")

#devtools::install_version("ggplot2", version = "3.4.2")
#devtools::install_version("Matrix", version = "1.5.4")

sce1 = as.SingleCellExperiment(seurat_object)
head(as.matrix(x_umap))



colnames(x_umap) = c('umap_1','umap_2')
colnames(x_scANVI) = c('scanvi_1','scanvi_2')
#sce1@int_colData@listData[["reducedDims"]]@listData[["scVI"]] = x_scvi
sce1@int_colData@listData[["reducedDims"]]@listData[["umap"]] = x_umap
sce1@int_colData@listData[["reducedDims"]]@listData[["scanvi"]] = x_umap
#human_cochlea_seurat <- as.Seurat(sce1)
#DimPlot(human_cochlea_seurat,reduction = "umap",group.by = "age_bins")

#sce1 = as.SingleCellExperiment(human_cochlea_seurat)
#sce1@int_colData@listData[["reducedDims"]]@listData[["scVI"]] = x_scvi
#sce1@int_colData@listData[["reducedDims"]]@listData[["umap"]] = x_umap
##############################################################################
############################################################################
#MILOR TO TEST DIFERENTIAL ABUNDANCE BETWEEN CONTROL AND NEOMYCIN CONDITIONS
library(miloR)

#############################################################################
#back-up counterpart from original NBT code
#without batch removal


#dimensionality reduction

#dimensionality reduction
# BY SCANVI 

set.seed(44)
remove(seurat_object)
sce2=sce1
table(sce2$cell_type_final)
#sce2 <- runUMAP(sce2, dimred="scANVI", ncomponents=2)
sce2$condition=sce2$age_bins
sce2$group=sce2$gw
scater::plotUMAP(sce2, colour_by="condition", point_alpha=1,  point_size=0.4)
scater::plotUMAP(sce2, colour_by="group", point_alpha=0.4,  point_size=0.4)
#sce2$cell_type=sce2$cell_type2
scater::plotUMAP(sce2, colour_by="cell_type_final", point_alpha=0.4, 
                 point_size=0.4, text_by='cell_type')
dev.off()

#Differential Abundance analysis with Milo
#We test for differential abundance between healthy and neomycin cochlea. 
#We start by defining neighbourhoods with refined sampling on the KNN graph.
#We inspect the size of neighbourhoods.
# defined by age bins <12; <28
sce2_milo <- Milo(sce2)
sce2_milo
## Build KNN graph
sce2_milo <- buildGraph(sce2_milo, d = 10, k=30)

## Compute neighbourhoods with refined sampling
sce2_milo <- makeNhoods(sce2_milo, k=30, d=10, prop = 0.2, refined=TRUE)


plotNhoodSizeHist(sce2_milo)
sce2_milo<-countCells(sce2_milo,
                      meta.data=as.data.frame(colData(sce2_milo)),
                      sample="sample")

head(nhoodCounts(sce2_milo))

scRNA_design<-data.frame(colData(sce2))[,c("sample","condition","group")]


scRNA_design$group<-as.factor(scRNA_design$group)
scRNA_design<-distinct(scRNA_design)
rownames(scRNA_design)<-scRNA_design$sample
scRNA_design
sce2_milo<-calcNhoodDistance(sce2_milo,d=10,reduced.dim = "scanvi")
#把批次和分组信息加到design中去
da_results <-testNhoods(sce2_milo,design = ~condition, design.df = scRNA_design)
head(da_results)

da_results%>%
  arrange(SpatialFDR)%>%
  head()


ggplot(da_results,aes(PValue)) + geom_histogram(bins=50)
ggplot(da_results,aes(logFC,-log10(SpatialFDR)))+
  geom_point()+
  geom_hline(yintercept=1)

scRNA <- buildNhoodGraph(sce2_milo)

##Plotsingle-cellUMAP
umap_pl<-plotReducedDim(scRNA,dimred="umap",
                          colour_by="group",text_by="cell_type_final",
                          text_size=3,point_size=0.5)+guides(fill="none")
umap_pl


##Plotneighbourhoodgraph
nh_graph_pl<-plotNhoodGraphDA(scRNA,da_results,
                                layout="umap",alpha=0.75)#alpha默认0.1

umap_pl+nh_graph_pl+plot_layout(guides="collect")


da_results<-annotateNhoods(scRNA,
                             da_results,
                             coldata_col="cell_type_final")
head(da_results)

ggplot(da_results,aes(cell_type_final_fraction))+geom_histogram(bins=50)

str(da_results)
table(da_results$cell_type_final)
range(da_results$SpatialFDR)
plotDAbeeswarm(da_results,group.by="cell_type_final",alpha=0.9)#alpha默认0.1


###################################################################################
###################################################################################
###################################################################################
# defined by age bins <12; <16;<20;<24;<28;
adata_loom <- connect(filename = "human_cochlea_annotation_v2.loom",
                      mode = "r+",skip.validate = TRUE)
matrix=adata_loom[["matrix"]][,]
matrix=t(matrix)
dim(matrix)
gene = adata_loom$row.attrs$var_names[]
barcode = adata_loom$col.attrs$obs_names[]

meta_data = read.csv('human_cochlea_annotation_obs_v2.csv',row.names = 1) # as form as dataframe format
meta_feature = read.csv('human_cochlea_annotation_var_v2.csv',row.names = 1)

colnames(matrix)= barcode
row.names(matrix)= gene
#x_scvi = adata_loom$col.attrs$X_scVI[,]
x_umap=read.csv("human_cochlea_annotation_umap.csv",row.names = 1)
x_scANVI=read.csv("human_cochlea_annotation_scANVI.csv",row.names = 1)

seurat_object= CreateSeuratObject(counts = matrix,meta.data = meta_data,
                                  project = 'human_cochlea_loom',
                                  min.cells = 0, 
                                  min.features = 0)
seurat_object@assays[["RNA"]]@meta.features <- meta_feature


rownames(x_umap) = barcode
colnames(x_umap) = c('UMAP_1','UMAP_2')
rownames(x_scANVI) = barcode
colnames(x_scANVI) = c('scANVI_1','scANVI_2')

seurat_object$cellid <- rownames(seurat_object@meta.data)
head(seurat_object$cellid)
rownames(x_umap)=seurat_object$cellid
rownames(x_scANVI)=seurat_object$cellid

seurat_object <-FindVariableFeatures(seurat_object)
seurat_object <- ScaleData(seurat_object)
seurat_object <- RunPCA(seurat_object)
seurat_object<- RunUMAP(seurat_object,dims = 1:30)
seurat_object <- FindNeighbors(seurat_object,dims = 1:30)
seurat_object <- FindClusters(seurat_object,resolution = 0.8)
DimPlot(seurat_object,reduction = "umap",group.by = "gw")

#devtools::install_version("ggplot2", version = "3.4.2")
#devtools::install_version("Matrix", version = "1.5.4")

sce1 = as.SingleCellExperiment(seurat_object)
head(as.matrix(x_umap))



colnames(x_umap) = c('umap_1','umap_2')
colnames(x_scANVI) = c('scanvi_1','scanvi_2')
#sce1@int_colData@listData[["reducedDims"]]@listData[["scVI"]] = x_scvi
sce1@int_colData@listData[["reducedDims"]]@listData[["umap"]] = x_umap
sce1@int_colData@listData[["reducedDims"]]@listData[["scanvi"]] = x_umap
set.seed(66)
remove(seurat_object)
sce2=sce1
table(sce2$cell_type_final)
#sce2 <- runUMAP(sce2, dimred="scANVI", ncomponents=2)
sce2$condition=sce2$age_bins
sce2$group=sce2$gw
scater::plotUMAP(sce2, colour_by="condition", point_alpha=1,  point_size=0.4)
scater::plotUMAP(sce2, colour_by="group", point_alpha=0.4,  point_size=0.4)
#sce2$cell_type=sce2$cell_type2
scater::plotUMAP(sce2, colour_by="cell_type_final", point_alpha=0.4, 
                 point_size=0.4, text_by='cell_type')
dev.off()

sce2_milo <- Milo(sce2)
sce2_milo
## Build KNN graph
sce2_milo <- buildGraph(sce2_milo, d = 10, k=30)

## Compute neighbourhoods with refined sampling
sce2_milo <- makeNhoods(sce2_milo, k=30, d=10, prop = 0.1, refined=TRUE)


plotNhoodSizeHist(sce2_milo)
sce2_milo<-countCells(sce2_milo,
                      meta.data=as.data.frame(colData(sce2_milo)),
                      sample="sample")

head(nhoodCounts(sce2_milo))

scRNA_design<-data.frame(colData(sce2))[,c("sample","condition","group")]


scRNA_design$group<-as.factor(scRNA_design$group)
scRNA_design<-distinct(scRNA_design)
rownames(scRNA_design)<-scRNA_design$sample
scRNA_design
sce2_milo<-calcNhoodDistance(sce2_milo,d=10,reduced.dim = "scanvi")



#把批次和分组信息加到design中去
da_results <-testNhoods(sce2_milo,design = ~condition, design.df = scRNA_design)
head(da_results)

da_results%>%
  arrange(SpatialFDR)%>%
  head()


ggplot(da_results,aes(PValue)) + geom_histogram(bins=50)
ggplot(da_results,aes(logFC,-log10(SpatialFDR)))+
  geom_point()+
  geom_hline(yintercept=1)

scRNA <- buildNhoodGraph(sce2_milo)

##Plotsingle-cellUMAP
umap_pl<-plotReducedDim(scRNA,dimred="umap",
                        colour_by="group",text_by="cell_type_final",
                        text_size=3,point_size=0.5)+guides(fill="none")
umap_pl


##Plotneighbourhoodgraph
nh_graph_pl<-plotNhoodGraphDA(scRNA,da_results,
                              layout="umap",alpha=0.75)#alpha默认0.1

umap_pl+nh_graph_pl+plot_layout(guides="collect")


da_results<-annotateNhoods(scRNA,
                           da_results,
                           coldata_col="cell_type_final")
head(da_results)

ggplot(da_results,aes(cell_type_final_fraction))+geom_histogram(bins=50)

str(da_results)
table(da_results$cell_type_final)
range(da_results$SpatialFDR)
plotDAbeeswarm(da_results,group.by="cell_type_final",alpha=0.9)#alpha默认0.1
remove(sce1)
write.csv(da_results,file="miloR_CELL_ABUNDANCE.csv")
saveRDS(sce2_milo,file = "miloR.rds")






