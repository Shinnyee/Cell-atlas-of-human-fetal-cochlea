# load packages
rm(list = ls())
library(Seurat)
library(dplyr)
library(reticulate)
library(sctransform)
library(cowplot)
library(ggplot2)
library(viridis)
library(tidyr)
library(magrittr)
library(reshape2)
library(readxl)
library(progeny)
library(readr)
library(stringr)
library(tidyverse)
options(stringsAsFactors = FALSE)
library(future)
plan("multicore", workers = 8)
options(future.globals.maxSize= 62914560000) # 60GB (60000*1024^2)
library(SingleCellExperiment)
library(Matrix)
#设定阈值
nFeature_lower <- 250
nFeature_upper <- 15000
nCount_lower <- 1000
nCount_upper <- 100000
pMT_lower <- 0
pMT_upper <- 30
pHB_lower <- 0
pHB_upper <- 5
HU_23W <- Read10X(data.dir = "C:/Users/Dell/Desktop/WORKPLACE/R/MATRICES/filtered_feature_bc_matrix_23W")
HU_25W <- Read10X(data.dir = "C:/Users/Dell/Desktop/WORKPLACE/R/MATRICES/filtered_feature_bc_matrix_25W")
HU_26W <- Read10X(data.dir = "C:/Users/Dell/Desktop/WORKPLACE/R/MATRICES/filtered_feature_bc_matrix_26W")
HU_23W = CreateSeuratObject(counts =HU_23W, project = "HU_23W")
HU_25W = CreateSeuratObject(counts =HU_25W, project = "HU_25W")
HU_26W = CreateSeuratObject(counts =HU_26W, project = "HU_26W")
mito_genes1 <- rownames(HU_23W)[grep("^MT-",rownames(HU_23W))]

HB_genes1 <- rownames(HU_23W)[grep("^HBA|^HBB",rownames(HU_23W))]

RP_genes1 <- rownames(HU_23W)[grep("^RPS|^RPL",rownames(HU_23W))]

HU_23W <- PercentageFeatureSet(HU_23W, pattern = "^HBA|^HBB", col.name = "pHB")
HU_25W <- PercentageFeatureSet(HU_25W, pattern = "^HBA|^HBB", col.name = "pHB")
HU_26W <- PercentageFeatureSet(HU_26W, pattern = "^HBA|^HBB", col.name = "pHB")
HU_23W <- PercentageFeatureSet(HU_23W, pattern = "^RPS|^RPL", col.name = "pRP")
HU_25W <- PercentageFeatureSet(HU_25W, pattern = "^RPS|^RPL", col.name = "pRP")
HU_26W <- PercentageFeatureSet(HU_26W, pattern = "^RPS|^RPL", col.name = "pRP")
HU_23W <- PercentageFeatureSet(HU_23W, pattern = "^MT-", col.name = "pMT")
HU_25W <- PercentageFeatureSet(HU_25W, pattern = "^MT-", col.name = "pMT")
HU_26W <- PercentageFeatureSet(HU_26W, pattern = "^MT-", col.name = "pMT")

qcparams <- c("nFeature_RNA", "nCount_RNA", "pHB", "pRP","pMT")
VlnPlot(object = HU_23W, features = qcparams, group.by = "orig.ident",
        pt.size = 0,ncol = 5)
VlnPlot(object = HU_25W, features = qcparams, group.by = "orig.ident", 
        pt.size = 0,ncol = 5)
VlnPlot(object = HU_26W, features = qcparams, group.by = "orig.ident", 
        pt.size = 0,ncol = 5)

## 过滤
HU_23W_filtered <- subset(HU_23W, subset = nFeature_RNA > nFeature_lower 
                          & nFeature_RNA < nFeature_upper  
                          & pMT < pMT_upper)

HU_23W
HU_23W_filtered
HU_25W_filtered <- subset(HU_25W, subset = nFeature_RNA > nFeature_lower 
                          & nFeature_RNA < nFeature_upper  
                          & pMT < pMT_upper)

HU_25W
HU_25W_filtered
HU_26W_filtered <- subset(HU_26W, subset = nFeature_RNA > nFeature_lower 
                          & nFeature_RNA < nFeature_upper 
                          & pMT < pMT_upper)

HU_26W
HU_26W_filtered
VlnPlot(object = HU_23W_filtered, features = qcparams, group.by = "orig.ident",
        pt.size = 0,ncol = 5)
VlnPlot(object = HU_25W_filtered, features = qcparams, group.by = "orig.ident",
        pt.size = 0,ncol = 5)
VlnPlot(object = HU_26W_filtered, features = qcparams, group.by = "orig.ident", 
        pt.size = 0,ncol = 5)
####################SEURAT CCA PIPELINE##################################
### CHOOSE SCTRNASFORMED DATA FOR INTEGRATION##############################
#merge data
HU_DATA <- merge(x = HU_23W_filtered, y = list(HU_25W_filtered,HU_26W_filtered),
                 add.cell.ids = c("Human_23W", "Human_25W","Human_26W"))
saveRDS(HU_DATA,file="HU_DATA_FILTERED.rds")
combined.list <- SplitObject(HU_DATA, split.by = "orig.ident")
rm(HU_DATA_FILTERED)
#combined.list<-lapply(X=combined.list,FUN = SCTransform, method="glmGamPoi")
for (i in 1:length(combined.list)) {
  combined.list[[i]] <- SCTransform(combined.list[[i]], 
                                    vars.to.regress = c("nCount_RNA","pMT"),
                                    verbose = TRUE, method="glmGamPoi")
}
features <- SelectIntegrationFeatures(object.list =combined.list, nfeatures = 2000)

#########################################data integration###################################################
combined.list <- PrepSCTIntegration(object.list =combined.list,
                                    anchor.features = features)
#combined.list <- lapply(X = combined.list, FUN = RunPCA, features = features)
cochlea.anchors <- FindIntegrationAnchors(object.list = combined.list, 
                                          normalization.method = "SCT", 
                                          anchor.features = features,
                                          dims = 1:40)
#cochlea.anchors <- FindIntegrationAnchors(object.list = combined.list, normalization.method = "SCT", anchor.features = features, reduction = "rpca", k.anchor = 20)
saveRDS(cochlea.anchors,file="cochlea_anchors_hu.rds")
cochlea.combined.sct <- IntegrateData(anchorset = cochlea.anchors,
                                      normalization.method = "SCT",dims = 1:40)
#cochlea.combined.sct <- IntegrateData(anchorset = cochlea.anchors, normalization.method = "SCT"， dims=1:50)
cochlea.combined.sct <- RunPCA(cochlea.combined.sct,npcs=200)
#cochlea.combined.sct <- RunPCA(cochlea.combined.sct, features = features, npcs=30)
ElbowPlot(cochlea.combined.sct , ndims = 200)
DefaultAssay(cochlea.combined.sct) <- "integrated"
cochlea.combined.sct <- FindNeighbors(cochlea.combined.sct, reduction = "pca",
                                      dims = 1:100)
cochlea.combined.sct <- FindClusters(cochlea.combined.sct)
cochlea.combined.sct <- RunUMAP(cochlea.combined.sct, reduction = "pca", 
                                dims = 1:100)
#cochlea.combined.sct <- RunTSNE(cochlea.combined.sct, reduction = "pca", dims = 1:30)
# Visualization
p1 <- DimPlot(cochlea.combined.sct, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(cochlea.combined.sct, reduction = "umap", label = TRUE)
p1+p2
table(cochlea.combined.sct$seurat_clusters)
DimPlot(cochlea.combined.sct, reduction = "umap", split.by = "orig.ident",label = TRUE)
DimPlot(cochlea.combined.sct, reduction = "umap",label = TRUE)
DefaultAssay(cochlea.combined.sct) <- "SCT"
FeaturePlot(cochlea.combined.sct,features = "TMC1",label = TRUE)
VlnPlot(cochlea.combined.sct,features = c("TMC1","SLC26A5","SLC17A8","OTOF","TBX2","OCM"),ncol = 1)
VlnPlot(cochlea.combined.sct,features = c("SNAP25","TUBB3","NEFL","NEFM","NEFH","MBP"),ncol = 1)
VlnPlot(cochlea.combined.sct,features = c("NEFH","CALB2","CALB1","RUNX1","LYPD1","ANO2"),ncol = 1)
VlnPlot(cochlea.combined.sct,features =c("NEFL","NEFM","SNAP25","SYT2","EPHA4","TSC22D3","RYR2",
                                         "GATA3","SMAD6","ATP2B4","ANO2","TH",
                                         "CALB2","MDGA1","CD24","SEMA3E","NTNG1","CALB1",
                                         "RUNX1","GRM8","POU4F1","LYPD1",
                                         "MBP","MPZ","MOBP"
),stack = TRUE,
fill.by = "ident",group.by = "seurat_clusters",
sort = FALSE,flip = TRUE) + theme(legend.position = "none")
saveRDS(cochlea.combined.sct,file = "HUMAN_cochlea_combined.rds")

cochlea.combined.sct=HUMAN_cochlea_combined
#find markers for every cluster compared to all remaining cells, report only the positive ones.
cochlea.combined.sct2 <- PrepSCTFindMarkers(cochlea.combined.sct)
all.markers <- FindAllMarkers(cochlea.combined.sct, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(all.markers,file = "HUMAN_all_markers_unidentified.csv")

# TO EXCLUDE AMBIGUOUS EXPRESSION PATTERNS FOR FURTHER ANNOTATION
a<-PrepSCTFindMarkers(cochlea.combined.sct)
all.markers2 <- FindAllMarkers(a, assay = "SCT", slot = "data", test.use = "roc")
write.csv(all.markers2,file = "HUMAN_all_markers_unidentified_roctest.csv")
top10 <- all.markers2 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_diff)
DoHeatmap(cochlea.combined.sct, features = top10$gene) + NoLegend()
# EXCLUDE CLUSTER 0#,16#
cochlea.combined.sct <- subset(cochlea.combined.sct,idents = c("0","16"),inver=TRUE)
dev.off()


DimPlot(cochlea.combined.sct, reduction = "umap",label = TRUE)
DefaultAssay(cochlea.combined.sct) <- "SCT"
markers.to.plot<-c("TMC1","SLC26A5","OTOF","MYO7A","SLC17A8","POU4F3",
                   "TSC22D1","EPYC","FBXO2","TMOD1","FGFR3","GATA3","GJB2",'DPYSL2',
                   "DGKB","CEACAM16","SLC1A3","S100A6","OTOG","OTOGL","FGF10",
                   "EMILIN2",	"AXIN2",	"TJP1",	"ATP1A2","OTOR","RARRES1","NOTUM",
                   "ENAH","CLU","NUDT4","PTGDS",
                   "ANXA1","KCNQ5","MEIS2","ADGRA3",
                   "NEFH","NEFL","SNAP25","PVALB","NEFM","TUBB3","CALB2",
                   "MBP","MOBP",	"MOG",	"TUBB4A",	"PMP22",	"PLP1",
                   "MPZ",
                   "DCT","TYR","MET",	"KCNJ13",	"KCNJ10",
                   "DCLK1","ESRRB","KCNQ1",	"KCNE1","DCN",	"ITIH5",
                   "EMCN","ARHGAP6","CLDN11","ATP6V0A4",
                   "VWF","MECOM",	"GBP7",
                   "CREB5","ITGA8",	"IFIT3",
                   "VEPH1","CHST9",	"POSTN",	"LUM","PDGFRB","MAP1B","NEBL",
                   "COCH","SLC4A10","SLC4A11",	"COL9A1",	"COL9A2",
                   "SLC7A11","SLC13A3","KCNK2",
                   "ATP2B2","SLC26A7","MEIS1","DKK2","SLC26A4",
                   "CEMIP","THSD4",
                   "RUNX2","COL1A2","RGS5","COL5A2","COL11A1","COL1A1","COL11A2",
                   "CD163","MEF2C","CX3CR1",	"CD68",	"C1QA","DIO2","SP7","TMEM119","IBSP",
                   "SLC4A4","GPC6","VMO1",
                   "CD3G","CD4","CD8A","SLC26A2"
)
DotPlot(cochlea.combined.sct, features = markers.to.plot, dot.scale = 8,
        cols  =c("white", "#ad9300")) + RotatedAxis()
### by using canonical markers that were from mouse species, we still hardly cannot annotate 3 clusters in which
# these clusters still showed robust distinct expressional profiles.
DotPlot(cochlea.combined.sct, features = markers.to.plot, cols = c("royalblue1", "maroon4", "sienna2"),
        dot.scale = 6, split.by = "orig.ident") +
  RotatedAxis()
#levels = c("PC1","PC2","TBC", "DC_OPC1","DC_OPC2",
#"DC_OPC3","IPh_IBC1","IPh_IBC2","HeC","OHC",
#"IHC","MC1","MC2","BC","IC",
#"CEC","PVM/M","Fibro1", "Fibro2", "Fibro3", 
#"Fibro4","Fib","SMC1","SMC2","SC1",
#"SC2","SC3","SC4","SC5","OL",
#"SGN1","SGN2","Chond/RMC","OB/M"))
#"HC","DC_OPC","IPh_IBC","TBC","PC",
#"HeC","SGN","OL","SC","IC",
#"MC","BC","CEC","PVM/M","Fib",
#"Fibro1","Fibro2","Fibro3","Fibro4","SMC",
#"OB/M","RMC"
table(Idents(cochlea.combined.sct))
####################################rename clusters#########################################################
Idents(cochlea.combined.sct) <- "seurat_clusters"
table(cochlea.combined.sct$seurat_clusters)
new.cluster.ids <- c("CEC", "SMC","BC", "GC", "Fibro1", "TBC", "GC", "Fibro2",
                     "Fibro1","Fibro1",
                     "Fibro1","IC","SMC","IC","GC","DC_PC_HeC","CEC","CEC",
                     "BC","SGN","Undefined_1","PVM/M","Undefined_2",
                     "IPh_BC","Immun","Fibro4","Fibro3",
                     "IC","Undefined_3","MC","CEC","Fib","SMC","HC")
names(new.cluster.ids) <- levels(cochlea.combined.sct)
cochlea.combined.sct<- RenameIdents(cochlea.combined.sct, new.cluster.ids)
DimPlot(cochlea.combined.sct, reduction = "umap",label = TRUE,repel = TRUE)
cochlea.combined.sct$subclasses <- Idents(cochlea.combined.sct)
Idents(cochlea.combined.sct) <- "subclasses"
table(Idents(cochlea.combined.sct))
DefaultAssay(cochlea.combined.sct) <- "integrated"
cochlea.combined.sct2 <- RunUMAP(cochlea.combined.sct, reduction = "pca", dims = 1:120)
DimPlot(cochlea.combined.sct2, reduction = "umap",label = TRUE,repel = TRUE)

cluster_cols <- c("#0868AC","#5AB2A8","#91D4CA","#F09150","#1A9850","#C7EAE5","#E47440",
                  "#D85730","#BF1D10","#A6D96A","#00C3FF","#6E4DA2","#3CA0C9","#288A82",
                  "#FDAE61","#2B8DBF","#085584","#4EB3D3","#b680a9","#60B85D","#01665E",
                  "#197AB5")#CB3A20
DimPlot(cochlea.combined.sct2,reduction = "umap",label = TRUE,repel = TRUE,
        group.by = "subclasses",cols = cluster_cols)
cluster_cols2 <- c( "#cc2d86","#6070b3","#7bc9a7")
p1 <- DimPlot(cochlea.combined.sct2,reduction = "umap",label = TRUE,repel = TRUE, 
              group.by = "orig.ident",cols = cluster_cols2,pt.size = 0.5)
p2 <- DimPlot(cochlea.combined.sct2,reduction = "umap",label = TRUE,repel = TRUE,
              group.by = "subclasses",cols = cluster_cols,pt.size = 0.5)
p2+p1
DimPlot(cochlea.combined.sct2, reduction = "umap", split.by = "orig.ident",
        label = TRUE,repel= TRUE,cols = cluster_cols)
saveRDS(cochlea.combined.sct2,file = "HUMAN_cochlea_combined_filtered.rds")
cochlea.combined.sct2=HUMAN_cochlea_combined_filtered
#########################BY USING METASCAPE WE WERE ABLE TO ASSIGN THOSE UNDEFINED CLUSTERS
######WE FOUND THAT UNDEFINED_1 AS ENDOTHELIAL CELLS (ECS)
######UNDEFINED_2 AS ERYTHROBLASTS (EBS)
######UNDEFINED_3 AS CYCLING CELLS (CCS)

Idents(cochlea.combined.sct2) <- "subclasses"
table(cochlea.combined.sct2$subclasses)
new.cluster.ids <- c("CEC", "SMC","BC", "GC", "Fibro1", "TBC", "Fibro2", 
                     "IC","DC_PC_HeC",
                     "SGN","ECs","PVM/M","EBs","IPh_BC","Immun","Fibro4","Fibro3",
                     "CCs","MC","Fib","HC")
names(new.cluster.ids) <- levels(cochlea.combined.sct2)
cochlea.combined.sct2<- RenameIdents(cochlea.combined.sct2, new.cluster.ids)
DimPlot(cochlea.combined.sct2, reduction = "umap",label = TRUE,repel = TRUE)
cochlea.combined.sct2$subclasses <- Idents(cochlea.combined.sct2)
Idents(cochlea.combined.sct2) <- "subclasses"
table(Idents(cochlea.combined.sct2))
DefaultAssay(cochlea.combined.sct2) <- "integrated"


cluster_cols <- c("#0868AC","#5AB2A8","#91D4CA","#F09150","#1A9850","#C7EAE5","#E47440",
                  "#D85730","#BF1D10","#A6D96A","#00C3FF","#6E4DA2","#3CA0C9","#288A82",
                  "#FDAE61","#2B8DBF","#085584","#4EB3D3","#b680a9","#60B85D","#01665E",
                  "#197AB5")#CB3A20
DimPlot(cochlea.combined.sct2,reduction = "umap",label = TRUE,repel = TRUE,
        group.by = "subclasses",cols = cluster_cols)
cluster_cols2 <- c( "#cc2d86","#6070b3","#7bc9a7")
p1 <- DimPlot(cochlea.combined.sct2,reduction = "umap",label = TRUE,repel = TRUE, 
              group.by = "orig.ident",cols = cluster_cols2,pt.size = 0.5)
p2 <- DimPlot(cochlea.combined.sct2,reduction = "umap",label = TRUE,repel = TRUE,
              group.by = "subclasses",cols = cluster_cols,pt.size = 0.5)
p2+p1
DimPlot(cochlea.combined.sct2, reduction = "umap", split.by = "orig.ident",
        label = TRUE,repel= TRUE,cols = cluster_cols)
saveRDS(cochlea.combined.sct2,file = "HUMAN_cochlea_combined_filtered.rds")




###### Dot plot for making cell-type specific markers across samples#########################################
DefaultAssay(cochlea.combined.sct2) <- "SCT"
table(Idents(cochlea.combined.sct2))
Idents(cochlea.combined.sct2) <- factor(Idents(cochlea.combined.sct2), levels = c("HC","DC_PC_HeC","IPh_IBC","TBC","SGN","GC","IC","MC",
                                                                                  "BC","PVM/M","CEC","Fib","Fibro1","Fibro2","Fibro3","Fibro4","SMC",
                                                                                  "Immun","ECs","EBs","CCs"))

cochlea.combined.sct2$subclasses2=Idents(cochlea.combined.sct2)

#find markers for every cluster compared to all remaining cells, report only the positive ones.
cochlea.combined.sct2 <- PrepSCTFindMarkers(cochlea.combined.sct2)
all.markers <- FindAllMarkers(cochlea.combined.sct2, only.pos = TRUE, 
                              min.pct = 0.25, logfc.threshold = 0.25)
write.csv(all.markers,file = "DEGs_human_cell_types_by_wilcox_test.csv")
all.markers2 <- FindAllMarkers(cochlea.combined.sct2, assay = "SCT", 
                               slot = "data", test.use = "roc")
write.csv(all.markers2,file = "DEGs_human_cell_types_by_roc_test.csv")

DefaultAssay(cochlea.combined.sct2) <- "RNA"
cochlea.combined.sct2 <-SCTransform(cochlea.combined.sct2,
                                    vars.to.regress = c("nCount_RNA","pMT"),
                                    verbose = TRUE, method="glmGamPoi")
Idents(cochlea.combined.sct2) <- factor(Idents(cochlea.combined.sct2), levels = c("HC","DC_PC_HeC","IPh_IBC","TBC","SGN","GC","IC","MC",
                                                                                  "BC","PVM/M","CEC","Fib","Fibro1","Fibro2","Fibro3","Fibro4","SMC",
                                                                                  "Immun","ECs","EBs","CCs"))

cochlea.combined.sct2$subclasses2=Idents(cochlea.combined.sct2)

markers.to.plot<-c("TMC1","PCP4","GJB2","GATA3","TSC22D1","EPYC","OTOG","OTOGL","EMILIN2","RARRES1","NEFL","SNAP25","MBP","PLP1",
                   "MPZ","PMP22","DCT","TYR","DCLK1","ESRRB","TJP1","CLDN11", "ITGA8",	"CREB5",   "VWF","MECOM",
                   "VEPH1","CHST9","COCH","SLC4A10","SLC7A11","SLC13A3","ATP2B2","SLC26A7",
                   "CEMIP","THSD4","RUNX2","COL1A2","CD163","MEF2C","IGFBP7","COL4A1",
                   "HBG2","HBM","TOP2A","HMGB2")
DotPlot(cochlea.combined.sct2, features = markers.to.plot, dot.scale = 8,
        group.by = "subclasses2",
        cols  =c("white", "#ad9300")) + RotatedAxis()
#cluster_cols2 <- c( "#cc2d86","#6070b3","#7bc9a7")
DotPlot(cochlea.combined.sct2, features = markers.to.plot, cols = c("#cc2d86","#6070b3","#7bc9a7"),
        dot.scale = 8, split.by = "orig.ident") +
  RotatedAxis()
DotPlot(cochlea.combined.sct, features = markers.to.plot, cols = c("royalblue1", "maroon4"),
        dot.scale = 6) +
  RotatedAxis()


#cochlea.combined.sct=HUMAN_cochlea_combined
DefaultAssay(cochlea.combined.sct2) <- "SCT"
markers.to.plot<-c("TMC1","SLC26A5","OTOF","MYO7A","SLC17A8","POU4F3",
                   "TSC22D1","EPYC","FBXO2","TMOD1","FGFR3","GATA3","GJB2",'DPYSL2',
                   "DGKB","CEACAM16","SLC1A3","S100A6","OTOG","OTOGL",
                   "EMILIN2",	"AXIN2",	"TJP1",	"ATP1A2","OTOR",
                   "ENAH","CLU","NUDT4","PTGDS",
                   "ANXA1","KCNQ5",
                   "NEFH","NEFL","SNAP25","PVALB","NEFM","TUBB3","CALB2",
                   "MBP","MOBP",	"MOG",	"TUBB4A",	"PMP22",	"PLP1",
                   "MPZ",
                   "DCT","TYR","MET",	"KCNJ13",	"KCNJ10",
                   "DCLK1","ESRRB","KCNQ1",	"KCNE1","DCN",	"ITIH5",
                   "EMCN","ARHGAP6","CLDN11",
                   "VWF","MECOM",	"GBP7",
                   "CREB5","ITGA8",	"IFIT3",
                   "VEPH1","CHST9",	"POSTN",	"LUM","PDGFRB","MAP1B","NEBL",
                   "COCH","SLC4A10","SLC4A11",	"COL9A1",	"COL9A2",
                   "SLC7A11","SLC13A3","KCNK2",
                   "ATP2B2","SLC26A7","MEIS1","DKK2",
                   "CEMIP","THSD4",
                   "RUNX2","COL1A2","RGS5","COL5A2",
                   "CD163","MEF2C","CX3CR1",	"CD68",	"C1QA",
                   "SLC4A4","GPC6","VMO1"
)
DotPlot(cochlea.combined.sct2, features = markers.to.plot, dot.scale = 8,
        cols  =c("white", "#ad9300")) + RotatedAxis()
FeaturePlot(cochlea.combined.sct2,features = c("TMC1","SLC26A5","SLC17A8","OTOF",
                                               "TBX2","KDM5B","IKZF2","PBX3"),
            label = TRUE)
FeaturePlot(cochlea.combined.sct2,features = c("TBX2","OTOF"),
            label = F,pt.size = 1)
VlnPlot(cochlea.combined.sct2,features = c("TMC1","SLC26A5","SLC17A8","OTOF",
                                           "TBX2","KDM5B","IKZF2","PBX3"
),stack = TRUE)#split.by = "orig.ident",
VlnPlot(cochlea.combined.sct2,features = c("LGR5","LGR6","PBX3","OTOGL","EGFR",'EGF'
),stack = TRUE,fill.by = "ident",flip = TRUE) + RotatedAxis()

table(cochlea.combined.sct2$orig.ident)



#############################################################################################################
#############################reproducibility of cell-type accuracy by MetaNeighbor analysis######################
#METANEIGHBOR ANALYSIS
dev.off()
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
HUMAN_cochlea_combined_filtered <- readRDS("C:/Users/Dell/Desktop/WORKPLACE/R/HUMAN_cochlea_combined_filtered.rds")
sce=cochlea.combined.sct2
sce=HUMAN_cochlea_combined_filtered
DefaultAssay(sce) <- "SCT"
sce$class_label = "CoE"
table(sce$class_label)
sce$sample_id <- colnames(sce)
sce$study_id <- sce$orig.ident# for each replicate
sce$study_id <- "human"# as uniform
sce$joint_subclass_label <- sce$subclasses
library(scater)
sce <- as.SingleCellExperiment(sce)
sce
library(MetaNeighbor)
dim(sce)
head(colData(sce))
table(sce$subclasses,sce$study_id)
global_hvgs <- variableGenes(dat=sce,exp_labels=sce$study_id)
length(global_hvgs)
global_hvgs=global_hvgs[1:500]
aurocs <- MetaNeighborUS(var_genes = global_hvgs,
                         
                         dat = sce,
                         
                         study_id = sce$study_id,
                         
                         cell_type = sce$subclasses,
                         
                         fast_version=TRUE)

plotHeatmap(aurocs,cex=0.5)
dev.off()
top_hits = topHits(aurocs,
                   
                   dat = sce,
                   
                   study_id = sce$study_id,
                   
                   cell_type = sce$subclasses,
                   
                   threshold = 0.9)
####################THESE PROCESSES COULD IGNORE##################
level1_split = splitClusters(aurocs, k = 2)
level1_split
first_split = level1_split[[2]]
full_labels = makeClusterName(sce$study_id, sce$subclasses)
subdata = sce[, full_labels %in% first_split]
dim(subdata)

var_genes = variableGenes(dat = subdata, exp_labels = subdata$study_id)
var_genes=var_genes[1:500]
aurocs = MetaNeighborUS(var_genes = var_genes,
                        dat = subdata, fast_version = TRUE,
                        study_id = subdata$study_id,
                        cell_type = subdata$subclasses)
plotHeatmap(aurocs, cex = 0.7)
dev.off()

level2_split = splitClusters(aurocs, k = 3)
my_split = level2_split[[3]]
subdata = sce[, full_labels %in% my_split]
var_genes = variableGenes(dat = subdata, exp_labels = subdata$study_id)
length(var_genes)
var_genes=var_genes[1:500]
aurocs = MetaNeighborUS(var_genes = var_genes,
                        dat = subdata, fast_version = TRUE,
                        study_id = subdata$study_id,
                        cell_type = subdata$subclasses)
plotHeatmap(aurocs, cex = 0.5)

################################################################################

best_hits = MetaNeighborUS(var_genes = global_hvgs,
                           dat = sce,
                           study_id = sce$study_id,
                           cell_type = sce$subclasses,
                           fast_version = TRUE,
                           one_vs_best = TRUE, symmetric_output = FALSE)
plotHeatmap(best_hits, cex = 0.5)
dev.off()
write.csv(best_hits,file = "best_hits_subclasses_human.csv")
d<- read.csv("best_hits_subclasses_human.csv",row.names = 1)
dev.off()
library(pheatmap)
library(viridis)
pheatmap(d,cluster_rows = FALSE,cluster_cols = FALSE,border=FALSE)
pheatmap::pheatmap(d,cluster_rows = FALSE,
                   cluster_cols = FALSE,border=FALSE,color = viridis(50))
#color = colorRampPalette(c("navy", "white", "firebrick3"))(100)
pheatmap(d,cluster_rows = FALSE,cluster_cols = FALSE,border=FALSE,color = colorRampPalette(c("navy", "white", "firebrick3"))(100))
pheatmap::pheatmap(d,cluster_rows = FALSE,cluster_cols = FALSE,border=FALSE,color = cividis(50))
dev.off()




mclusters <- extractMetaClusters(best_hits,threshold=0.7)

mcsummary <- scoreMetaClusters(mclusters,best_hits)

plotUpset(mclusters)

plotMetaClusters(mclusters,best_hits)
dev.off()
cluster_graph = makeClusterGraph(best_hits, low_threshold = 0.3)
plotClusterGraph(cluster_graph, sce$study_id,
                 sce$subclasses, size_factor=3)


####################################stacked barplot for average proportion##################################

library(dittoSeq)
cochlea.combined.sct2$subclasses.species <- paste(Idents(cochlea.combined.sct2), cochlea.combined.sct2$orig.ident, sep = "_")
cochlea.combined.sct2$subclasses <- Idents(cochlea.combined.sct2)

dittoBarPlot(object = cochlea.combined.sct2,var = "orig.ident",group.by = "subclasses")
dittoBarPlot(object = cochlea.combined.sct2,var = "subclasses",group.by = "orig.ident")


####################################correlation of each clusters############################################
av <-AverageExpression(cochlea.combined.sct2,
                       group.by = "subclasses",
                       assays = "SCT")
av=av[[1]]
#av <- ScaleData(av)
av <- log1p(av)
head(av)
range(av)
cg=names(tail(sort(apply(av, 1, sd)),500))
pheatmap::pheatmap(cor(av[cg,]),color = colorRampPalette(c("#008080", "white", "#ad9300"))(100))
dev.off()
qcparams <- c("nFeature_RNA", "nCount_RNA", "pHB", "pRP")
sce=cochlea.combined.sct2
sce <- subset(sce, subset = nFeature_RNA > nFeature_lower 
              & nFeature_RNA < nFeature_upper 
              & pMT < pMT_upper)
VlnPlot(object = sce, features = "nFeature_RNA", group.by = "subclasses", pt.size = 0,split.by = "orig.ident")
VlnPlot(object = sce, features = "nCount_RNA", group.by = "subclasses", pt.size = 0,split.by = "orig.ident")
VlnPlot(object = sce, features = "pMT", group.by = "subclasses", pt.size = 0,split.by = "orig.ident")
VlnPlot(object = sce, features = "pRP", group.by = "subclasses", pt.size = 0,split.by = "orig.ident")


##############heatmap showing top 4 to 6 genes and GO terms
table(cochlea.combined.sct2$subclasses.species)
Idents(cochlea.combined.sct2) <- "subclasses"

table(cochlea.combined.sct2$subclasses)
Idents(cochlea.combined.sct2) <- "subclasses"
Idents(cochlea.combined.sct2) <- factor(cochlea.combined.sct2$subclasses,levels =
                                          c("HC","DC_PC_HeC","IPh_BC","TBC","SGN","GC","IC","MC",
                                            "BC","PVM/M","CEC","Fib","Fibro1","Fibro2","Fibro3","Fibro4","SMC",
                                            "Immun","CCs","ECs","EBs"))
table(Idents(cochlea.combined.sct2))
cochlea.combined.sct2$celltype = Idents(cochlea.combined.sct2)
table(cochlea.combined.sct2$celltype)
cochlea.combined.sct2<- PrepSCTFindMarkers(cochlea.combined.sct2,assay = "SCT",verbose = TRUE)

cochlea_genes <- FindAllMarkers(cochlea.combined.sct2, assay = "SCT", slot = "data", test.use = "roc")
cochlea_genes <- cochlea_genes[which(cochlea_genes$avg_diff > 0), ]
#write.csv(cochlea_genes,file = "degs_HUMAN_cochlea.csv")


library(ClusterGVis)
library(org.Hs.eg.db)
library(clusterProfiler)
group <- data.frame(gene=cochlea_genes$gene,
                    group= cochlea_genes$cluster)
Gene_ID <- bitr(cochlea_genes$gene,fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = "org.Hs.eg.db")
data  <- merge(Gene_ID,group,by.x='SYMBOL',by.y='gene')
enrich2 <- compareCluster(
  ENTREZID~group, 
  data=data, 
  fun="enrichGO", 
  OrgDb="org.Hs.eg.db",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.5
)#qvalueCutoff = 0.05
enrich2_discription <- enrich2@compareClusterResult[!duplicated(enrich2@compareClusterResult[,3]),] %>%
  dplyr::group_by(Cluster) %>%
  dplyr::top_n(n=5,wt=Count) %>% dplyr::top_n(n=5,wt=ID)
enrich2_discription <- enrich2_discription[!duplicated(enrich2_discription[,3]),]
enrich2_discription <- enrich2_discription[,c(3,2,4,7,5)]
write.csv(enrich2_discription,file = "goenrich_HUMAN_subclasses.csv")
enrich2_discription <- read.csv("goenrich_HUMAN_subclasses.csv",row.names = 1)



##### plot GO TERMS RELATED TO SUBCLASSES
library(ggplot2)
df <- read.csv("GO_TERMS_RELATED_TO_SUBCLASSES.csv", header = T)
df$LogP <- -df$Log.q.value.

df$labelx=rep(0,nrow(df))
df$labely=seq(nrow(df),1)

ggplot(data = df, 
       aes(LogP, reorder(Description,LogP))) +
  geom_bar(stat="identity",
           alpha=0.5,
           fill="#FE8D3C",
           width = 0.5) + 
  geom_text(aes(x=labelx,
                y=labely,
                label = Description),
            size=3, 
            hjust =0.001)+
  theme_classic()+
  theme(axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(colour = 'black', linewidth = 1),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.ticks.x = element_line(colour = 'black', linewidth = 1),
        axis.title.x = element_text(colour = 'black', size = 10))+
  xlab("-log10(qvalue)")+
  ggtitle("Enrichment")+
  scale_x_continuous(expand = c(0,0))


# add gene name
markGenes = unique(cochlea_genes$gene)[sample(1:length(unique(cochlea_genes$gene)),50,
                                              replace = F)]
# no average cells
cochlea.markers1 <- cochlea_genes %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 5, wt = avg_diff)
# retain duplicate diff gene in multiple clusters
sce=cochlea.combined.sct2

DefaultAssay(sce) <- "SCT"
Idents(sce) <- "subclasses"
st.data <- prepareDataFromscRNA(object = sce,
                                diffData = cochlea.markers1,
                                showAverage = FALSE)
str(st.data)
st.data1 <- prepareDataFromscRNA(object = cochlea.combined.sct2,
                                 diffData = cochlea.markers1,
                                 group.by = "subclasses",
                                 showAverage = FALSE)
str(st.data1)
enrich <- enrichCluster(object = st.data1,
                        OrgDb = org.Hs.eg.db,
                        type = "BP",
                        organism = "Hsa",
                        pvalueCutoff = 0.5,
                        topn = 5,
                        seed = 5201314,
                        readable = TRUE)
# add GO annotation
#pdf('HEATMAP_SUBCLASS_human_DEGS_GO.pdf',height = 20,width = 30,onefile = F)
library(viridis)
visCluster(object = st.data,
           plot.type = "both",
           column_title_rot = 45,
           markGenes = unique(cochlea.markers1$gene),
           markGenes.side = "left",
           annoTerm.data = enrich,
           
           go.size = 5,
           border = FALSE,
           line.side = "left",
           cluster.order = c(1:18),go.col = "black",
           add.bar = T)
#pdf('sc2.pdf',height = 10,width = 14,onefile = F)
visCluster(object = st.data1,
           plot.type = "both",
           column_title_rot = 45,
           show_row_dend = F,
           markGenes = unique(cochlea.markers1$gene),
           
           annoTerm.data = enrich,
           cluster.order = c(1:18),
           markGenes.side = "left",
           genes.gp = c('italic',fontsize = 12,col = "orange"),
           go.size = 5
           
)

dev.off()



