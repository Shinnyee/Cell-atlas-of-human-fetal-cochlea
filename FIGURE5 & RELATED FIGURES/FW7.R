# 
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
hu <- Read10X(data.dir = "F:/PROJECTS/PROJECT_HUMAN FETAL COCHLEAE/WORKPLACE/source/FW7_5")
sce=hu
sce<- CreateSeuratObject(counts =sce)
sce$orig.ident <- "human_GW7"
sce
nFeature_lower <- 500
nFeature_upper <- 5000

pMT_lower <- 0
pMT_upper <- 10
pHB_lower <- 0
pHB_upper <- 5

sce
mito_genes <- rownames(sce)[grep("^MT-",rownames(sce))]
sce <- PercentageFeatureSet(sce, pattern = "^MT-", col.name = "pMT")
sce <- subset(sce, subset = nFeature_RNA > nFeature_lower & 
                nFeature_RNA < nFeature_upper & 
                
                pMT < pMT_upper
)
sce


sce <- SCTransform(sce, verbose = T, vars.to.regress = c("nCount_RNA","pMT"), 
                   conserve.memory = T)
sce <- RunPCA(sce,npcs = 50)

DimPlot(sce,reduction = "pca")
ElbowPlot(sce , ndims = 50)
sce <- FindNeighbors(sce, reduction = "pca", dims = 1:50)
sce <- FindClusters(sce)
sce <- RunUMAP(sce, reduction = "pca", dims = 1:50)

# Visualization

DimPlot(sce, reduction = "umap",label = TRUE)

saveRDS(sce, file = "cochlea_GW7.rds")
sce=cochlea_GW7
DefaultAssay(sce) <- "SCT"


VlnPlot(sce, features = c("RORB","ISL1","LGR5","SOX2","FGF20",
                          "TECTA","FGF10","JAG1","BMP4","GATA3",
                          "OTX2","FGF9","WNT4","GSC"))



markers.to.plot<-c("RORB","ISL1","LGR5","SOX2","FGF20","ATOH1","CDKN1B",   # COCHLEAR DUCT FLOOR PROSENSORY
                   "TECTA","FGF10","JAG1",#FLOOR MEDIAL
                   "GATA3","FGFR3","PROX1","BMP4",   #FLOOR LATERAL
                   "OTX2",	"FGF9",	"WNT4",	"GSC",#COCHLEAR ROOF CELLS
                   
                   
                   "MEIS2",	"ADAMTSL1",	"OTOGL",	"USH1C",#VESTIBULAR SIPPORTING CELLS/EPITHELIAL CELLS
                   
                   "NTN1","SMOC2","WNT3",#VESTIBULAR ROOF CELLS
                   "KCNE1","ATP1B2","SPP1",#DARK CELLS
                   "STRC","OTOF","USH2A","MYO15A","SLC26A5","SLC17A8","SLC7A14",    #HAIR CELLS
                   
                   "SNAP25","TUBB3","PRPH","CALB1","PBX3","ESRRG", # SGN
                   "MBP","MPZ","PLP1","PMP22", # GC
                   "EPCAM", #EPITHELIUM
                   "PRRX1",#MESENCHYMAL
                   "ACAN", #CHONDROCYTES
                   "PECAM1", #ENDOTHELIAL
                   "PTPRC", #MACROPHAGES
                   "MLANA", #melanocytes
                   "MKI67","TOP2A","HMGB2" #CYCLING
                   
                   
                   
                   
)







DotPlot(sce,features=markers.to.plot)+ RotatedAxis()

table(Idents(sce))


sce <- subset(sce,idents = c("0","1","2","3","9","13","14","16","19",
                             "20","22","26"),invert=T)



# ANNOTATION BASED ON CELL REPORTS PAPER
# RENAME clusters 'cluster_label'
Idents(sce) <- "seurat_clusters"
table(Idents(sce))
new.cluster.ids <- c("cochlear_epithelium","Chondrocytes","CCs","Chondrocytes","SGN","Mesenchymal",
                     "Vestibular_SCs","GC","Vestibular_roof_cells","Mesenchymal","Endothelial_cells",
                     "SGN","Macrophages","SGN","Vestibular_HCs","Melanocytes"
                     
)
names(new.cluster.ids) <- levels(sce)
sce<- RenameIdents(sce, new.cluster.ids)
DimPlot(sce, reduction = "umap",label = TRUE,repel = TRUE)
sce$celltype <- Idents(sce)
Idents(sce) <- "celltype"
table(Idents(sce))
DimPlot(sce, reduction = "umap",label = TRUE,repel = TRUE)


saveRDS(sce, file = "cochlea_GW7.rds")

sce=cochlea_GW7
Idents(sce) <- "celltype"
sce <- RunUMAP(sce, reduction = "pca", dims = 1:20)
DimPlot(sce, reduction = "umap",label = TRUE)


sce <- FindNeighbors(sce, reduction = "pca", dims = 1:50)
sce <- FindClusters(sce,resolution = 1.2)

p1=DimPlot(sce, reduction = "umap",label = TRUE,repel = TRUE,group.by = "celltype")
p1
p2=DimPlot(sce, reduction = "umap",label = TRUE,repel = TRUE,group.by = "seurat_clusters")
p2
p1+p2
table(sce$celltype,sce$seurat_clusters)
markers.to.plot<-c("RORB","ISL1","LGR5","SOX2","FGF20","ATOH1","CDKN1B",   # COCHLEAR DUCT FLOOR PROSENSORY
                   "TECTA","FGF10","JAG1",#FLOOR MEDIAL
                   "GATA3","FGFR3","PROX1","BMP4",   #FLOOR LATERAL
                   "OTX2",	"FGF9",	"WNT4",	"GSC",#COCHLEAR ROOF CELLS
                   
                   
                   "MEIS2",	"ADAMTSL1",	"OTOGL",	"USH1C",#VESTIBULAR SIPPORTING CELLS/EPITHELIAL CELLS
                   
                   "NTN1","SMOC2","WNT3",#VESTIBULAR ROOF CELLS
                   "KCNE1","ATP1B2","SPP1",#DARK CELLS
                   "STRC","OTOF","USH2A","MYO15A","SLC26A5","SLC17A8","SLC7A14",    #HAIR CELLS
                   
                   "SNAP25","TUBB3","PRPH","CALB1","PBX3","ESRRG", # SGN
                   "MBP","MPZ","PLP1","PMP22", # GC
                   "EPCAM", #EPITHELIUM
                   "PRRX1",#MESENCHYMAL
                   "ACAN", #CHONDROCYTES
                   "PECAM1", #ENDOTHELIAL
                   "PTPRC", #MACROPHAGES
                   "MLANA", #melanocytes
                   "MKI67","TOP2A","HMGB2" #CYCLING
                   
                   
                   
                   
)







DotPlot(sce,features=markers.to.plot)+ RotatedAxis()


# RENAME clusters "subclasses"
Idents(sce) <- "seurat_clusters"
table(Idents(sce))
new.cluster.ids <- c("Chondrocytes","CoE_prosensory","Mesenchymal","GC",
                     "CCs","Chondrocytes","CCs","Chondrocytes","SGN",
                     "CoE_foor_roof_cells","Vestibular_SCs","Mesenchymal",
                     "SGN","Endothelial","SGN","Vestibular_Roof_cells","Vestibular_Roof_cells",
                     "Vestibular_SCs","Macrophages","Chondrocytes",
                     "Mixed","Vestibular_SCs","Vestibular_HCs","Melanocytes"
                     
)
names(new.cluster.ids) <- levels(sce)
sce<- RenameIdents(sce, new.cluster.ids)
DimPlot(sce, reduction = "umap",label = TRUE,repel = TRUE)
sce$subclasses <- Idents(sce)
Idents(sce) <- "subclasses"
table(Idents(sce))
DimPlot(sce, reduction = "umap",label = TRUE,repel = TRUE)
cluster_cols <- c("#0868AC","#5AB2A8","#91D4CA","#F09150","#1A9850","#C7EAE5","#E47440",
                  "#D85730","#BF1D10","#A6D96A","#00C3FF","#6E4DA2","#3CA0C9","#288A82",
                  "#FDAE61","#2B8DBF","#085584","#4EB3D3","#b680a9","#60B85D","#01665E",
                  "#197AB5")#CB3A20
DimPlot(sce,reduction = "umap",label = TRUE,repel = TRUE,
        group.by = "subclasses",cols = cluster_cols)
Idents(sce)<- "subclasses"
table(Idents(sce))
sce <- subset(sce,idents = "Mixed",invert=T)
Idents(sce) <- factor(Idents(sce), levels = c("CoE_prosensory","CoE_foor_roof_cells",
                                              "Vestibular_HCs","Vestibular_SCs",
                                              "Vestibular_Roof_cells",
                                              "SGN","GC",
                                              "Chondrocytes","Mesenchymal","Endothelial",
                                              "Melanocytes","CCs","Macrophages"
))
sce$subclasses=Idents(sce)


markers.to.plot<-c( "GATA3",   "RORB","ISL1","LGR5","SOX2","FGF20",   # COCHLEAR DUCT FLOOR PROSENSORY
                   "TECTA","FGF10","JAG1",#FLOOR MEDIAL
                      #FLOOR LATERAL
                   #COCHLEAR ROOF CELLS
                   "STRC","OTOF","USH2A","MYO15A",    #HAIR CELLS
                   
                   "MEIS2",	"ADAMTSL1",	"OTOGL",	"USH1C",#VESTIBULAR SIPPORTING CELLS/EPITHELIAL CELLS
                   
                   "NTN1","SMOC2","WNT3",#VESTIBULAR ROOF CELLS
                  
                  
                   
                   "SNAP25","TUBB3","PRPH","CALB1","PBX3","ESRRG", # SGN
                   "MBP","MPZ","PLP1","PMP22", # GC
                   "ACAN", #CHONDROCYTES
                   "PRRX1",#MESENCHYMAL
                   "PECAM1", #ENDOTHELIAL
                   "MLANA", #melanocytes
                   "MKI67","TOP2A","HMGB2", #CYCLING
                  
                   
                   "PTPRC" #MACROPHAGES
                 
                 
                   
                   
                   
                   
)
DotPlot(sce, features = markers.to.plot, dot.scale = 8, group.by = "subclasses",
        cols  =c("white", "#ad9300")) + RotatedAxis()


saveRDS(sce, file = "cochlea_GW7.rds")
##TRANSFER TO H5AD

library(sceasy)
library(reticulate)
use_condaenv('EnvironmentName')

DefaultAssay(sce) <- "RNA"
sceasy::convertFormat(sce, from="seurat", to="anndata",
                      outFile='Hu_GW7_python.h5ad')





