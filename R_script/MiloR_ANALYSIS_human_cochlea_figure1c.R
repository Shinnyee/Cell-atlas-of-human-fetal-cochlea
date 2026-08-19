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
adata_loom <- connect(filename = "human_cochlea_annotation_REVISION.loom",
                      mode = "r+",skip.validate = TRUE)
matrix=adata_loom[["matrix"]][,]
matrix=t(matrix)
dim(matrix)
gene = adata_loom$row.attrs$var_names[]
barcode = adata_loom$col.attrs$obs_names[]

meta_data = read.csv('human_cochlea_annotation_obs_REVISION.csv',row.names = 1) # as form as dataframe format
meta_feature = read.csv('human_cochlea_annotation_var_REVISION.csv',row.names = 1)

colnames(matrix)= barcode
row.names(matrix)= gene
#x_scvi = adata_loom$col.attrs$X_scVI[,]
x_umap=read.csv("human_cochlea_annotation_umap_REVISION.csv",row.names = 1)
x_scANVI=read.csv("human_cochlea_annotation_scANVI_REVISION.csv",row.names = 1)

seurat_object= CreateSeuratObject(counts = matrix,meta.data = meta_data,
                                      project = 'human_cochlea_loom',
                                      min.cells = 0, 
                                      min.features = 0)
seurat_object@assays[["RNA"]]@meta.features <- meta_feature


rownames(x_umap) = barcode
colnames(x_umap) = c('UMAP_1','UMAP_2')
rownames(x_scANVI) = barcode
colnames(x_scANVI) = c('scANVI_1','scANVI_2','scANVI_3','scANVI_4','scANVI_5','scANVI_6',
                       'scANVI_7','scANVI_8','scANVI_9','scANVI_10','scANVI_11','scANVI_12',
                       'scANVI_13','scANVI_14','scANVI_15','scANVI_16','scANVI_17','scANVI_18',
                       'scANVI_19','scANVI_20','scANVI_21','scANVI_22','scANVI_23','scANVI_24',
                       'scANVI_25','scANVI_26','scANVI_27','scANVI_28','scANVI_29','scANVI_30')

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
DimPlot(seurat_object,reduction = "umap",group.by = "cell_type_final")
dev.off()
#devtools::install_version("ggplot2", version = "3.4.2")
#devtools::install_version("Matrix", version = "1.5.4")

sce1 = as.SingleCellExperiment(seurat_object)
head(as.matrix(x_umap))



#colnames(x_umap) = c('umap_1','umap_2')
#colnames(x_scANVI) = c('scanvi_1','scanvi_2')
#sce1@int_colData@listData[["reducedDims"]]@listData[["scVI"]] = x_scvi
sce1@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] = x_umap
sce1@int_colData@listData[["reducedDims"]]@listData[["scanvi"]] = x_scANVI
plotUMAP(sce1, colour_by="cell_type_final")
#human_cochlea_seurat <- as.Seurat(sce1)
#DimPlot(human_cochlea_seurat,reduction = "umap",group.by = "age_bins")
human_cochlea_seurat <- as.Seurat(sce1)
DimPlot(human_cochlea_seurat,reduction = "UMAP",group.by = "cell_type_final",raster=FALSE)

sce = as.SingleCellExperiment(human_cochlea_seurat)
plotUMAP(sce, colour_by="cell_type_final")
dev.off()
sce@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] = x_umap
plotUMAP(sce, colour_by="age_bins")

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
sce2=sce
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
sce2_milo <- buildGraph(sce2_milo, 
                        k = 30, 
                        d = 30, 
                        reduced.dim = "PCA")

## Compute neighbourhoods with refined sampling
sce2_milo <- makeNhoods(sce2_milo, k=30, d=30, prop = 0.2, refined=TRUE)


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
sce2_milo<-calcNhoodDistance(sce2_milo,d=10,reduced.dim = "pca")
#?????κͷ?????Ϣ?ӵ?design??ȥ
da_results <-testNhoods(sce2_milo,design = ~condition, design.df = scRNA_design)
head(da_results)

da_results%>%
  arrange(SpatialFDR)%>%
  head()
cat("\nSignificant neighbourhoods (FDR < 0.1):", 
    sum(da_results$SpatialFDR < 0.1, na.rm = TRUE), "\n")
cat("Total neighbourhoods:", nrow(da_results), "\n")

ggplot(da_results,aes(PValue)) + geom_histogram(bins=50)
ggplot(da_results,aes(logFC,-log10(SpatialFDR)))+
  geom_point()+
  geom_hline(yintercept=1)

scRNA <- buildNhoodGraph(sce2_milo)

##Plotsingle-cellUMAP
umap_pl<-plotReducedDim(scRNA,dimred="UMAP",
                          colour_by="group",text_by="cell_type_final",
                          text_size=3,point_size=0.5)+guides(fill="none")
umap_pl


##Plotneighbourhoodgraph
nh_graph_pl<-plotNhoodGraphDA(scRNA,da_results,
                                layout="UMAP",alpha=0.3)#alphaĬ??0.1

umap_pl+nh_graph_pl+plot_layout(guides="collect")

dev.off()
da_results<-annotateNhoods(scRNA,
                             da_results,
                             coldata_col="cell_type_final")


head(da_results)

ggplot(da_results,aes(cell_type_final_fraction))+geom_histogram(bins=50)

str(da_results)
table(da_results$cell_type_final)
range(da_results$SpatialFDR)
plotDAbeeswarm(da_results,group.by="cell_type_final",alpha=0.3)#alphaĬ??0.1

#####################################################################################
######################################################################################
#7. 
write.csv(da_results, "milo_da_results.csv", row.names = TRUE)
######################################################################################
# ==============================================================
# 
# ==============================================================

library(miloR)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)



# ==============================================================
# 1. 
# ==============================================================
write.csv(da_results, "milo_da_results.csv", row.names = TRUE)
print("✅ DA results exported: milo_da_results.csv")

# ==============================================================
# 2. 
# ==============================================================



if (!is.null(scRNA@nhoodGraph)) {
  
  nhood_coords <- as.data.frame(scRNA@nhoodGraph@nhoodGraph)
  colnames(nhood_coords) <- c("X", "Y")
  nhood_coords$nhood <- rownames(nhood_coords)
  
 
  write.csv(nhood_coords, "milo_nhood_coords.csv", row.names = FALSE)
  print("✅ Nhood coordinates exported: milo_nhood_coords.csv")
} else {
  print("⚠️ nhoodGraph not found, using alternative method...")
}


if (!exists("nhood_coords")) {
  
  nhood_idx <- which(colData(scRNA)$nhood_ixs_refined == 1)
  nhood_cells <- rownames(colData(scRNA))[nhood_idx]
  
 
  umap_coords <- reducedDim(scRNA, "UMAP")
  rownames(umap_coords) <- rownames(colData(scRNA))
  
  
  nhood_coords <- as.data.frame(umap_coords[nhood_cells, ])
  colnames(nhood_coords) <- c("X", "Y")
  nhood_coords$nhood <- rownames(nhood_coords)
  
  write.csv(nhood_coords, "milo_nhood_coords.csv", row.names = FALSE)
  print("✅ Nhood coordinates exported (from UMAP): milo_nhood_coords.csv")
}

# ==============================================================
# 3. 
# ==============================================================

# 
nhood_sizes <- colSums(nhoods(scRNA))
nhood_size_df <- data.frame(
  nhood = names(nhood_sizes),
  Nhood_size = as.numeric(nhood_sizes)
)
write.csv(nhood_size_df, "milo_nhood_sizes.csv", row.names = FALSE)
print("✅ Nhood sizes exported: milo_nhood_sizes.csv")

# ==============================================================
# 4.
# ==============================================================


nhood_anno <- da_results[, c("cell_type_final", "cell_type_final_fraction")]
nhood_anno$nhood <- rownames(nhood_anno)
write.csv(nhood_anno, "milo_nhood_annotations.csv", row.names = FALSE)
print("✅ Nhood annotations exported: milo_nhood_annotations.csv")

# ==============================================================
# 5. 
# ==============================================================


if ("cell_type_final_colors" %in% names(metadata(scRNA))) {
  cell_colors <- metadata(scRNA)$adult_covid19_pbmc_annot_colors
  color_df <- data.frame(
    cell_type = names(cell_colors),
    color = cell_colors
  )
  write.csv(color_df, "milo_cell_colors.csv", row.names = FALSE)
  print("✅ Cell colors exported: milo_cell_colors.csv")
}

# ==============================================================
# 6. 
# ==============================================================

umap_all <- as.data.frame(reducedDim(scRNA, "UMAP"))
colnames(umap_all) <- c("UMAP1", "UMAP2")
umap_all$cell <- rownames(umap_all)
umap_all$condition <- colData(scRNA)$condition
umap_all$cell_type <- colData(scRNA)$cell_type_final

write.csv(umap_all, "milo_umap_coords_all.csv", row.names = FALSE)
print("✅ All cell UMAP coords exported: milo_umap_coords_all.csv")


saveRDS(sce2_milo,file = "miloR.rds")

################################################################################
# ==============================================================
#
# ==============================================================

# 1. 
nhood_mat <- nhoods(scRNA)

# 2.
umap_coords <- reducedDim(scRNA, "UMAP")

# 3. 
#
umap_matrix <- as.matrix(umap_coords)
nhood_sums <- t(nhood_mat) %*% umap_matrix

nhood_counts <- colSums(nhood_mat)

# 4. 
nhood_centers <- nhood_sums / nhood_counts

# 5. 
nhood_coords <- as.data.frame(nhood_centers)
colnames(nhood_coords) <- c("X", "Y")
nhood_coords$nhood_id <- 1:nrow(nhood_coords)
nhood_coords$nhood <- paste0("nhood_", 1:nrow(nhood_coords))

# 6. 
head(nhood_coords)
write.csv(nhood_coords, "milo_nhood_coords.csv", row.names = FALSE)
cat("✅ Saved: milo_nhood_coords.csv\n")

# 7. 
da_results$nhood_id <- as.numeric(rownames(da_results))
da_results <- da_results %>%
  left_join(nhood_coords %>% select(nhood_id, X, Y), by = "nhood_id")

write.csv(da_results, "milo_da_results_with_coords.csv", row.names = TRUE)
cat("✅ Saved: milo_da_results_with_coords.csv\n")
save.image(file = 'miloR_ANALYSIS_MD_pbmc.RData')















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



#?????κͷ?????Ϣ?ӵ?design??ȥ
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
                              layout="umap",alpha=0.75)#alphaĬ??0.1

umap_pl+nh_graph_pl+plot_layout(guides="collect")


da_results<-annotateNhoods(scRNA,
                           da_results,
                           coldata_col="cell_type_final")
head(da_results)

ggplot(da_results,aes(cell_type_final_fraction))+geom_histogram(bins=50)

str(da_results)
table(da_results$cell_type_final)
range(da_results$SpatialFDR)
plotDAbeeswarm(da_results,group.by="cell_type_final",alpha=0.9)#alphaĬ??0.1
remove(sce1)
write.csv(da_results,file="miloR_CELL_ABUNDANCE.csv")
saveRDS(sce2_milo,file = "miloR.rds")






