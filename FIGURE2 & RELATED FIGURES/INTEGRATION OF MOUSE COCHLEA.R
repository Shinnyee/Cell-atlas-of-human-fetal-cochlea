library(Seurat)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(sctransform)
library(tidyverse)
library(harmony)
options(stringsAsFactors = FALSE)
library(tidyverse)
library(DT)
library(biomaRt)
# 读取文件，MOUSE ID TO HUMAN ID
genesV2 <- read.table("mouse_to_human_genes.txt", sep="\t", header=T)
# mouse dataset
# create Seurat object
mouse_p28_1<- Read10X(data.dir = "F:/Chai Lab/PROJECT_INTEGRATION ACROSS SPECIES/Cochlea_Evo_new/MATRICES_MOUSE/P28_1")
mouse_p28_1<- CreateSeuratObject(counts =mouse_p28_1)
mouse_p28_1$orig.ident <- "mouse_p28_1"
mouse_p28_2<- Read10X(data.dir = "F:/Chai Lab/PROJECT_INTEGRATION ACROSS SPECIES/Cochlea_Evo_new/MATRICES_MOUSE/P28_2")
mouse_p28_2<- CreateSeuratObject(counts =mouse_p28_2)
mouse_p28_2$orig.ident <- "mouse_p28_2"
mouse_p25_1<- Read10X(data.dir = "F:/Chai Lab/PROJECT_INTEGRATION ACROSS SPECIES/Cochlea_Evo_new/MATRICES_MOUSE/mm1_filtered_feature_bc_matrix_P25")
mouse_p25_1<- CreateSeuratObject(counts =mouse_p25_1)
mouse_p25_1$orig.ident <- "mouse_p25_1"
mouse_p25_2<- Read10X(data.dir = "F:/Chai Lab/PROJECT_INTEGRATION ACROSS SPECIES/Cochlea_Evo_new/MATRICES_MOUSE/mm2_filtered_feature_bc_matrix_P25")
mouse_p25_2<- CreateSeuratObject(counts =mouse_p25_2)
mouse_p25_2$orig.ident <- "mouse_p25_2"
## Extract Expression Data
mouse_counts_p28_1 <- as.matrix(GetAssayData(mouse_p28_1, slot = "counts"))
mouse_counts_p28_1 <- data.frame(gene=rownames(mouse_counts_p28_1), mouse_counts_p28_1, check.names = F)
dim(mouse_counts_p28_1)
mouse_counts_p28_2 <- as.matrix(GetAssayData(mouse_p28_2, slot = "counts"))
mouse_counts_p28_2 <- data.frame(gene=rownames(mouse_counts_p28_2), mouse_counts_p28_2, check.names = F)
dim(mouse_counts_p28_2)
mouse_counts_p25_1 <- as.matrix(GetAssayData(mouse_p25_1, slot = "counts"))
mouse_counts_p25_1 <- data.frame(gene=rownames(mouse_counts_p25_1), mouse_counts_p25_1, check.names = F)
dim(mouse_counts_p25_1)
mouse_counts_p25_2 <- as.matrix(GetAssayData(mouse_p25_2, slot = "counts"))
mouse_counts_p25_2 <- data.frame(gene=rownames(mouse_counts_p25_2), mouse_counts_p25_2, check.names = F)
dim(mouse_counts_p25_2)

# gene conversion from mouse to human orthology
mouse_counts_p28_1$Gene <- genesV2[match(mouse_counts_p28_1$gene, genesV2[,1]),2]
mouse_counts_p28_1 <- subset(mouse_counts_p28_1, Gene!='NA')
mouse_counts_p28_1 <- dplyr::select(mouse_counts_p28_1, Gene, everything())
mouse_counts_p28_1 <- mouse_counts_p28_1[, !(colnames(mouse_counts_p28_1) %in% 'gene')]
dim(mouse_counts_p28_1)
write.csv(mouse_counts_p28_1,file="rawcounts_mouse_counts_p28_1_to_human_transformed.csv")
mouse_counts_p28_2$Gene <- genesV2[match(mouse_counts_p28_2$gene, genesV2[,1]),2]
mouse_counts_p28_2 <- subset(mouse_counts_p28_2, Gene!='NA')
mouse_counts_p28_2 <- dplyr::select(mouse_counts_p28_2, Gene, everything())
mouse_counts_p28_2 <- mouse_counts_p28_2[, !(colnames(mouse_counts_p28_2) %in% 'gene')]
dim(mouse_counts_p28_2)
write.csv(mouse_counts_p28_2,file="rawcounts_mouse_counts_p28_2_to_human_transformed.csv")
mouse_counts_p25_1$Gene <- genesV2[match(mouse_counts_p25_1$gene, genesV2[,1]),2]
mouse_counts_p25_1 <- subset(mouse_counts_p25_1, Gene!='NA')
mouse_counts_p25_1 <- dplyr::select(mouse_counts_p25_1, Gene, everything())
mouse_counts_p25_1 <- mouse_counts_p25_1[, !(colnames(mouse_counts_p25_1) %in% 'gene')]
dim(mouse_counts_p25_1)
write.csv(mouse_counts_p25_1,file="rawcounts_mouse_counts_p25_1_to_human_transformed.csv")
mouse_counts_p25_2$Gene <- genesV2[match(mouse_counts_p25_2$gene, genesV2[,1]),2]
mouse_counts_p25_2 <- subset(mouse_counts_p25_2, Gene!='NA')
mouse_counts_p25_2 <- dplyr::select(mouse_counts_p25_2, Gene, everything())
mouse_counts_p25_2 <- mouse_counts_p25_2[, !(colnames(mouse_counts_p25_2) %in% 'gene')]
dim(mouse_counts_p25_2)
write.csv(mouse_counts_p25_2,file="rawcounts_mouse_counts_p25_2_to_human_transformed.csv")
#create mouse_data into Seuratobject
rm(list=ls())
mouse_p28_1<-read.csv(file="rawcounts_mouse_counts_p28_1_to_human_transformed.csv")
rownames(mouse_p28_1)<-mouse_p28_1[,1]
mouse_p28_1<-mouse_p28_1[,-1]
head(rownames(mouse_p28_1))
rownames(mouse_p28_1)<-mouse_p28_1[,1]
mouse_p28_1<-mouse_p28_1[,-1]
head(rownames(mouse_p28_1))
mouse_p28_1 <- CreateSeuratObject(counts = mouse_p28_1, 
                    
                          min.cells = 10, 
                          min.features = 10)
mouse_p28_1[['percent.mt']] <- PercentageFeatureSet(mouse_p28_1, pattern = '^MT-')
mouse_p28_1 <- subset(mouse_p28_1, subset = nFeature_RNA > 500 & nFeature_RNA < 10000 & percent.mt < 15 )
mouse_p28_1$orig.ident <- "Mouse_P28_1"
mouse_p28_1
mouse_p28_2<-read.csv(file="rawcounts_mouse_counts_p28_2_to_human_transformed.csv")
rownames(mouse_p28_2)<-mouse_p28_2[,1]
mouse_p28_2<-mouse_p28_2[,-1]
head(rownames(mouse_p28_2))
rownames(mouse_p28_2)<-mouse_p28_2[,1]
mouse_p28_2<-mouse_p28_2[,-1]
head(rownames(mouse_p28_2))
mouse_p28_2 <- CreateSeuratObject(counts = mouse_p28_2, 
                                  
                                  min.cells = 10, 
                                  min.features = 10)
mouse_p28_2[['percent.mt']] <- PercentageFeatureSet(mouse_p28_2, pattern = '^MT-')
mouse_p28_2 <- subset(mouse_p28_2, subset = nFeature_RNA > 500 & nFeature_RNA < 10000 & percent.mt < 15 )
mouse_p28_2$orig.ident <- "Mouse_P28_2"
mouse_p28_2
mouse_p25_1<-read.csv(file="rawcounts_mouse_counts_p25_1_to_human_transformed.csv")
rownames(mouse_p25_1)<-mouse_p25_1[,1]
mouse_p25_1<-mouse_p25_1[,-1]
head(rownames(mouse_p25_1))
rownames(mouse_p25_1)<-mouse_p25_1[,1]
mouse_p25_1<-mouse_p25_1[,-1]
head(rownames(mouse_p25_1))
mouse_p25_1 <- CreateSeuratObject(counts = mouse_p25_1, 
                                  
                                  min.cells = 10, 
                                  min.features = 10)
mouse_p25_1[['percent.mt']] <- PercentageFeatureSet(mouse_p25_1, pattern = '^MT-')
mouse_p25_1 <- subset(mouse_p25_1, subset = nFeature_RNA > 500 & nFeature_RNA < 10000 & percent.mt < 15 )
mouse_p25_1$orig.ident <- "Mouse_P25_1"
mouse_p25_1
mouse_p25_2<-read.csv(file="rawcounts_mouse_counts_p25_2_to_human_transformed.csv")
rownames(mouse_p25_2)<-mouse_p25_2[,1]
mouse_p25_2<-mouse_p25_2[,-1]
head(rownames(mouse_p25_2))
rownames(mouse_p25_2)<-mouse_p25_2[,1]
mouse_p25_2<-mouse_p25_2[,-1]
head(rownames(mouse_p25_2))
mouse_p25_2 <- CreateSeuratObject(counts = mouse_p25_2, 
                                  
                                  min.cells = 10, 
                                  min.features = 10)
mouse_p25_2[['percent.mt']] <- PercentageFeatureSet(mouse_p25_2, pattern = '^MT-')
mouse_p25_2 <- subset(mouse_p25_2, subset = nFeature_RNA > 500 & nFeature_RNA < 10000 & percent.mt < 15 )
mouse_p25_2$orig.ident <- "Mouse_P25_2"
mouse_p25_2
#merge sample
sce_all <- merge(x = mouse_p25_1, 
                 y = c(mouse_p25_2, mouse_p28_1,mouse_p28_2),
                 add.cell.ids = c("Mouse_P25_1","Mouse_P25_2","Mouse_P28_1","Mouse_P28_2"))
sce_all
combined.list <- SplitObject(sce_all, split.by = "orig.ident")
combined.list<-lapply(X=combined.list,FUN = SCTransform, method="glmGamPoi")
features <- SelectIntegrationFeatures(object.list =combined.list, nfeatures = 3000)
#Need to increase memory for server usage
library(future)
options(future.globals.maxSize = 10000 * 1024^2)
#########################################data integration###################################################
combined.list <- PrepSCTIntegration(object.list =combined.list, anchor.features = features)
combined.list <- lapply(X = combined.list, FUN = RunPCA, features = features)
cochlea.anchors <- FindIntegrationAnchors(object.list = combined.list, normalization.method = "SCT", anchor.features = features)
#cochlea.anchors <- FindIntegrationAnchors(object.list = combined.list, normalization.method = "SCT", anchor.features = features, reduction = "rpca", k.anchor = 20)
saveRDS(cochlea.anchors,file="cochlea_anchors_MOUSE.rds")
cochlea_anchors_MOUSE <- readRDS("H:/Chai Lab/PROJECT_INTEGRATION ACROSS SPECIES/Cochlea_Evo_new/mouse_snRNA_v2/cochlea_anchors_MOUSE.rds")
cochlea.combined.sct <- IntegrateData(anchorset = cochlea_anchors_MOUSE, normalization.method = "SCT",dims = 1:40)
#cochlea.combined.sct <- IntegrateData(anchorset = cochlea.anchors, normalization.method = "SCT"， dims=1:50)
cochlea.combined.sct <- RunPCA(cochlea.combined.sct,npcs=200)
#cochlea.combined.sct <- RunPCA(cochlea.combined.sct, features = features, npcs=30)
ElbowPlot(cochlea.combined.sct , ndims = 200)
cochlea.combined.sct <- FindNeighbors(cochlea.combined.sct, reduction = "pca", dims = 1:50)
cochlea.combined.sct <- FindClusters(cochlea.combined.sct)
cochlea.combined.sct <- RunUMAP(cochlea.combined.sct, reduction = "pca", dims = 1:100)
#cochlea.combined.sct <- RunTSNE(cochlea.combined.sct, reduction = "pca", dims = 1:100)
# Visualization
p1 <- DimPlot(cochlea.combined.sct, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(cochlea.combined.sct, reduction = "umap", label = TRUE)
p1+p2

#find markers for every cluster compared to all remaining cells, report only the positive ones.
cochlea.combined.sct2 <- PrepSCTFindMarkers(cochlea.combined.sct)
all.markers2 <- FindAllMarkers(cochlea.combined.sct2, assay = "SCT", slot = "data", test.use = "roc")
write.csv(all.markers2,file = "MOUSE_all_markers_unidentified_roctest.csv")
top10 <- all.markers2 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_diff)
DoHeatmap(cochlea.combined.sct2, features = top10$gene) + NoLegend()
# EXCLUDE CLUSTER 5#,16#,23#
cochlea.combined.sct2 <- subset(cochlea.combined.sct2,idents = c("5","16","23"),inver=TRUE)
dev.off()
DimPlot(cochlea.combined.sct2, reduction = "umap",label = TRUE)
cochlea.combined.sct2=MOUSE_cochlea_combined_filtered_pnas
DefaultAssay(cochlea.combined.sct2) <- "SCT"
table(Idents(cochlea.combined.sct2))


Idents(cochlea.combined.sct2) <- factor(Idents(cochlea.combined.sct2), levels = 
                     c("HC","DC_PC_HeC","IPh_IBC","TBC","IDC",     "SGN","GC","IC","MC",
            "BC", "SpC_RC",    "PVM/M","Fib","Fibro1","Fibro2","Fibro3","Fibro4","OB/M",
                      "Macrophages"))



markers.to.plot<-c("TMC1","MYO7A","GJB2","GATA3","TSC22D1","EPYC","OTOG","OTOGL","EMILIN2","RARRES1","NEFL","SNAP25","MBP","PLP1",
                   "MPZ","PMP22","DCT","TYR","DCLK1","ESRRB","TJP1","CLDN11", "ITGA8",	"CREB5",   "VWF","MECOM",
                   "VEPH1","CHST9","COCH","SLC4A10","SLC7A11","SLC13A3","ATP2B2","SLC26A7",
                   "CEMIP","THSD4","RUNX2","COL1A2","CD163","MEF2C"
                   )
DotPlot(cochlea.combined.sct2, features = markers.to.plot, dot.scale = 8,
        cols  =c("white", "#ad9300")) + RotatedAxis()



markers.to.plot<-c("OTOA", #INTERDENTAL CELLS
                   "TMC1", #HAIR CELLS
                   "SMPX","LGR6", #PILLAR CELLS
                   "EDNRB","PRDM12",#PILLAR CELLS, DEITERS' CELLS
                   "CEACAM16","RBP7","FABP3",#DEITERS' CELLS
                   "AQP4","TECTB", #INNER BOR.-PH./HENSEN'S CELLS
                   "EPYC", #INNER BOR.-PH./HENSEN'S CELLS; CLAUD./IN.-OUT.SULCUS CELLS 
                   "OTOG","OTOGL","USH1C",    # SUPPORTING CELLS
                   "SNAP25","NEFL","PVALB",#SGN
                   "OLIG1","OLIG2", #GLIAL PRECURSOR CELLS
                   "MAG","MOG","MOBP", #SATELLITE GLIAL CELLS
                   "MPZ","PRX", #SCHWANN CELLS
                   "ESAM", #ENDO.CELLS
                   "RGS5", #PERICYTES
                   "TAGLN", #SMOOTH MUSCLE CELLS (ALDN1A2)
                   "KCNE1", "ALDN1A2",  #MARGINAL CELLS
                   "TYR", #INTER.STRIA CELLS
                   "CLDN11","ATP6V0A4",    #BASAL STRIA CELLS
                   "SLC26A4","ANXA1", #SPINDLE CELLS; ROOT CELLS
                   "OTOS","CAR3", "COL9A2","COL9A3",    #FIBROCYTES
                   "CAR2","SLC12A2","GJB2",# FIBROCYTES SUBTYPES I TO V
                   "DLX5","BGLAP","PHEX", "DMP1","IBSP",   #OSTEOBLASTS
                   "OSR2","BMP6","SLC7A11","CAVIN2",    # SURROUNDING STRUCTURES
                   "SLC26A7", #RM
                   "EMILIN2","NOTUM","RARRES1", #TYMPANIC BORDER CELLS
                   "RHD", #ERYTHROCYTES
                   "LY6G", #NEUTROPHILS
                   "CD19", # B CELLS
                   "CD3G", #T CELLS
                   "NKG7", #NK CELLS
                   "CD68", #MONOCYTES
                   "AIF1", #MACROPHAGES
                   "PRSS34" #MAST CELLS
                   
                   
)
DotPlot(cochlea.combined.sct2, features = markers.to.plot, dot.scale = 8,
        cols  =c("white", "#ad9300")) + RotatedAxis()


# EXCLUDE CLUSTER 26#,29#
cochlea.combined.sct2 <- subset(cochlea.combined.sct2,idents = c("26","29"),inver=TRUE)
dev.off()
DimPlot(cochlea.combined.sct2, reduction = "umap",label = TRUE)

##rename clusters
Idents(cochlea.combined.sct2) <- "seurat_clusters"
table(cochlea.combined.sct2$seurat_clusters)
new.cluster.ids <- c("TBC","PVM/M","DC_PC_HeC","MC","Fibro3","HC","SGN","Fibro4",
                     "IPh_IBC","Fibro1","IC","Macrophages","BC","BC","IDC","SGN",
                     "DC_PC_HeC","SpC_RC","IPh_IBC","HC","Fibro2","DC_PC_HeC","Fib",
                     "GC","IPh_IBC","OB/M","SGN")
names(new.cluster.ids) <- levels(cochlea.combined.sct2)
cochlea.combined.sct2<- RenameIdents(cochlea.combined.sct2, new.cluster.ids)
DimPlot(cochlea.combined.sct2, reduction = "umap",label = TRUE,repel = TRUE)
cochlea.combined.sct2$subclasses <- Idents(cochlea.combined.sct2)

Idents(cochlea.combined.sct2) <- "subclasses"
table(Idents(cochlea.combined.sct2))
DefaultAssay(cochlea.combined.sct2) <- "integrated"
cochlea.combined.sct2 <- RunUMAP(cochlea.combined.sct2, reduction = "pca", dims = 1:50)
DimPlot(cochlea.combined.sct2, reduction = "umap",label = TRUE,repel = TRUE)

cluster_cols <- c("#5089a0","#01665e","#00abe0","#91d4ca","#d85730","#6edbf4","#a6d96a",
                  "#cb3a20","#165f7c","#f09150","#c7eae5","#4b6662","#5ab2a8","#636370",
                  "#243835","#e47440","#fdae61","#416b3f","#6e4da2"
                  )#CB3A20
DimPlot(cochlea.combined.sct2,reduction = "umap",label = TRUE,repel = TRUE,
        group.by = "subclasses",cols = cluster_cols)
cluster_cols2 <- c( "#ad7e5b","#bfac62","#42859f","#6bc7c2")
table(cochlea.combined.sct2$orig.ident)
p1 <- DimPlot(cochlea.combined.sct2,reduction = "umap",label = TRUE,repel = TRUE, 
              group.by = "orig.ident",cols = cluster_cols2,pt.size = 0.5)

p1

DimPlot(cochlea.combined.sct2, reduction = "umap", split.by = "orig.ident",
        label = TRUE,repel= TRUE,cols = cluster_cols)
saveRDS(cochlea.combined.sct2,file = "MOUSE_cochlea_combined_filtered_pnas.rds")
MOUSE_cochlea_combined_filtered_pnas <- readRDS("F:/PROJECTS/PROJECT_INTEGRATION ACROSS SPECIES/Cochlea_Evo_new/mouse_snRNA_v2/MOUSE_cochlea_combined_filtered_pnas.rds")

####################################stacked barplot for average proportion##################################
#BiocManager::install("dittoSeq")
library(dittoSeq)
dittoBarPlot(object = cochlea.combined.sct2,var = "orig.ident",group.by = "subclasses")
dittoBarPlot(object = cochlea.combined.sct2,var = "subclasses",group.by = "orig.ident")


##rename subclasses
Idents(cochlea.combined.sct2) <- "subclasses"
table(cochlea.combined.sct2$subclasses)
new.cluster.ids <- c("TBC","PVM/M","DC_PC_HeC","MC","Fibro3","HC","SGN",
                     "Fibro4","IPh_IBC","Fibro1","IC","Immun","BC",
                     "IDC","SpC_RC","Fibro2","Fib","GC","OB"
)
names(new.cluster.ids) <- levels(cochlea.combined.sct2)
cochlea.combined.sct2<- RenameIdents(cochlea.combined.sct2, new.cluster.ids)
DimPlot(cochlea.combined.sct2, reduction = "umap",label = TRUE,repel = TRUE)
cochlea.combined.sct2$subclasses <- Idents(cochlea.combined.sct2)
saveRDS(cochlea.combined.sct2,file = "MOUSE_cochlea_combined_filtered_pnas.rds")

MOUSE_cochlea_combined_filtered_pnas <- readRDS("H:/Chai Lab/PROJECT_INTEGRATION ACROSS SPECIES/Cochlea_Evo_new/mouse_snRNA_v2/MOUSE_cochlea_combined_filtered_pnas.rds")
cochlea.combined.sct2=MOUSE_cochlea_combined_filtered_pnas
DefaultAssay(cochlea.combined.sct2) <- "SCT"
Idents(cochlea.combined.sct2) <- factor(Idents(cochlea.combined.sct2), levels = 
                                          c("HC","DC_PC_HeC","IPh_IBC","TBC","IDC",     "SGN","GC","IC","MC",
                                            "BC", "SpC_RC",    "PVM/M","Fib","Fibro1","Fibro2","Fibro3","Fibro4","OB",
                                            "Immun"))

markers.to.plot<-c("TMC1","PCP4","GJB2","GATA3","TSC22D1","EPYC","OTOG","OTOGL","EMILIN2","RARRES1","NEFL","SNAP25","MBP","PLP1",
                   "MPZ","PMP22","DCT","TYR","DCLK1","ESRRB","TJP1","CLDN11", "ITGA8",	"CREB5",   "VWF","MECOM",
                   "VEPH1","CHST9","COCH","SLC4A10","SLC7A11","SLC13A3","ATP2B2","SLC26A7",
                   "CEMIP","THSD4","RUNX2","COL1A2","CD163","MEF2C","IGFBP7","COL4A1",
                   "HBG2","HBM","TOP2A","HMGB2","GSDMD")
DotPlot(cochlea.combined.sct2, features = markers.to.plot, dot.scale = 8,
        cols  =c("white", "#ad9300")) + RotatedAxis()
markers.to.plot<-c("OTOA", #INTERDENTAL CELLS
                   "TMC1", #HAIR CELLS
                   "SMPX","LGR6", #PILLAR CELLS
                   "EDNRB","PRDM12",#PILLAR CELLS, DEITERS' CELLS
                   "CEACAM16","RBP7","FABP3",#DEITERS' CELLS
                   "AQP4","TECTB", #INNER BOR.-PH./HENSEN'S CELLS
                   "EPYC", #INNER BOR.-PH./HENSEN'S CELLS; CLAUD./IN.-OUT.SULCUS CELLS 
                   "OTOG","OTOGL","USH1C",    # SUPPORTING CELLS
                   "SNAP25","NEFL","PVALB",#SGN
                   "OLIG1","OLIG2", #GLIAL PRECURSOR CELLS
                   "MAG","MOG","MOBP", #SATELLITE GLIAL CELLS
                   "MPZ","PRX", #SCHWANN CELLS
                   "ESAM", #ENDO.CELLS
                   "RGS5", #PERICYTES
                   "TAGLN", #SMOOTH MUSCLE CELLS (ALDN1A2)
                   "KCNE1", "ALDN1A2",  #MARGINAL CELLS
                   "TYR", #INTER.STRIA CELLS
                   "CLDN11","ATP6V0A4",    #BASAL STRIA CELLS
                   "SLC26A4","ANXA1", #SPINDLE CELLS; ROOT CELLS
                   "OTOS","CAR3", "COL9A2","COL9A3",    #FIBROCYTES
                   "CAR2","SLC12A2","GJB2",# FIBROCYTES SUBTYPES I TO V
                   "DLX5","BGLAP","PHEX", "DMP1","IBSP",   #OSTEOBLASTS
                   "OSR2","BMP6","SLC7A11","CAVIN2",    # SURROUNDING STRUCTURES
                   "SLC26A7", #RM
                   "EMILIN2","NOTUM","RARRES1", #TYMPANIC BORDER CELLS
                   "RHD", #ERYTHROCYTES
                   "LY6G", #NEUTROPHILS
                   "CD19", # B CELLS
                   "CD3G", #T CELLS
                   "NKG7", #NK CELLS
                   "CD68", #MONOCYTES
                   "AIF1", #MACROPHAGES
                   "PRSS34" #MAST CELLS
                   
                   
)
DotPlot(cochlea.combined.sct2, features = markers.to.plot, dot.scale = 8,
        cols  =c("white", "#ad9300")) + RotatedAxis()
