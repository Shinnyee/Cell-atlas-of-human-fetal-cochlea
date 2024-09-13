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
hu <- Read10X(data.dir = "F:/PROJECTS/PROJECT_HUMAN FETAL COCHLEAE/WORKPLACE/source/rawdata_GW15")
sce=hu
sce<- CreateSeuratObject(counts =sce)
sce$orig.ident <- "human_GW15"
sce
nFeature_lower <- 500
nFeature_upper <- 10000
nCount_lower <- 1000
nCount_upper <- 100000
pMT_lower <- 0
pMT_upper <- 15
pHB_lower <- 0
pHB_upper <- 5

sce
mito_genes <- rownames(sce)[grep("^MT-",rownames(sce))]
sce <- PercentageFeatureSet(sce, pattern = "^MT-", col.name = "pMT")
sce <- subset(sce, subset = nFeature_RNA > nFeature_lower & 
                nFeature_RNA < nFeature_upper & 
                nCount_RNA > nCount_lower & 
                nCount_RNA < nCount_upper &
                pMT < pMT_upper
)
sce


sce <- SCTransform(sce, verbose = T, vars.to.regress = c("nCount_RNA","pMT"), conserve.memory = T)
sce <- RunPCA(sce,npcs = 30)

DimPlot(sce,reduction = "pca")
ElbowPlot(sce , ndims = 50)
sce <- FindNeighbors(sce, reduction = "pca", dims = 1:30)
sce <- FindClusters(sce)
sce <- RunUMAP(sce, reduction = "pca", dims = 1:30)

# Visualization

DimPlot(sce, reduction = "umap",label = TRUE)
seq <- seq(0.1, 2, by = 0.1)
for(res in seq){
  sce <- FindClusters(sce, resolution = res)
}
sce <- RunUMAP(sce, reduction = "pca", dims = 1:30)
library(clustree)
library(patchwork)
p1 <- clustree(sce, prefix = 'SCT_snn_res.') + coord_flip()
p2 <- DimPlot(sce, group.by = 'SCT_snn_res.0.8', label = T)
p1 + p2 + plot_layout(widths = c(3, 1))

sce <- FindClusters(sce,resolution = 1.4)
sce <- RunUMAP(sce, reduction = "pca", dims = 1:30)
DimPlot(sce, reduction = "umap",label = TRUE)

saveRDS(sce, file = "cochlea_GW15.rds")
sce=cochlea_GW15
DefaultAssay(sce) <- "SCT"
VlnPlot(sce, features = c("TMC1","OTOF","ATOH1","INSM1","SOX2","CABP2",
                          "FGFR3","PROX1","FGF8","MYO7A",
                          "SLC17A8","FGF20","SLC26A5","OCM"))#hc

VlnPlot(sce, features = c("LGR5","LGR6","OTOG","OTOGL","SOX2"))


markers.to.plot<-c("TMC1","INSM1","OTOF","MYO7A","CABP2","ATOH1",
                   "TSC22D1","EPYC","FBXO2","TMOD1","FGFR3","GATA3","GJB2",'DPYSL2',
                   "DGKB","SLC1A3","S100A6","OTOG","OTOGL","LGR5","LGR6",
                   "EMILIN2",	"AXIN2",	"TJP1",	"ATP1A2","OTOR",
                   "ENAH","CLU","NUDT4","PTGDS",
                   "ANXA1","KCNQ5","MEIS2",
                   "NEFH","NEFL","SNAP25","PVALB","NEFM","TUBB3","CALB2",
                   "MBP",	"TUBB4A",	"PMP22",	"PLP1",
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
                   "CD163","MEF2C",	"CD68",	"C1QA",
                   "SLC4A4","GPC6","VMO1",
                   "CD79A","CD79B",
                   "CD4","CD8A"
)
markers.to.plot<-c("TMC1","INSM1","OTOF","MYO7A","CABP2","ATOH1",# HC
                   "TSC22D1","EPYC","FBXO2","TMOD1","FGFR3","GATA3","GJB2",'DPYSL2',#DCOPC
                   "DGKB","SLC1A3","S100A6","OTOG","OTOGL","LGR5","LGR6",#IPHIBC
                   "EMILIN2",	"AXIN2",	"TJP1",	"ATP1A2","OTOR",#TBC
                   
                  
                   "MBP",	"TUBB4A",	"PMP22",	"PLP1","MPZ",#SC/OL
                  
                   "DCT","TYR","MET",	"KCNJ13",	"KCNJ10",#IC
                   "DCLK1","ESRRB","KCNQ1",	"KCNE1","DCN",	"ITIH5",#MC
                   "EMCN","ARHGAP6","CLDN11",#BC
                   "VWF","MECOM",	"GBP7",#CEC
                   
                   "ATP2B2","SLC26A7","MEIS1","DKK2",#FIBRO3 
                   "SMOC2","OTOA","PTGDS","CDKN1C", # Interdental
                   "UBE2C","KIAA0101", #Ube2c+ 
                   "OC90","TTR","VMO1","ECEL1","ENPEP" #OC90
)







DotPlot(sce,features=markers.to.plot)+ RotatedAxis()


# ANNOTATION BASED ON NC PAPER
# RENAME clusters 'cluster_label'
Idents(sce) <- "seurat_clusters"
table(Idents(sce))
new.cluster.ids <- c("TBC","BC","BC","MC","BC","Oc90",
                     "BC","TBC","Interdental","SC/OL","BC",
                     "DCOPC","SC/OL","Ube2c+","BC","IPhIBC",
                     "HC","TBC","CEC","BC","IC"
)
names(new.cluster.ids) <- levels(sce)
sce<- RenameIdents(sce, new.cluster.ids)
DimPlot(sce, reduction = "umap",label = TRUE,repel = TRUE)
sce$celltype <- Idents(sce)
Idents(sce) <- "celltype"
table(Idents(sce))
DimPlot(sce, reduction = "umap",label = TRUE,repel = TRUE)


saveRDS(sce, file = "cochlea_GW15.rds")

sce=cochlea_GW15
Idents(sce) <- "celltype"
table(Idents(sce))
DimPlot(sce, reduction = "umap",label = TRUE,repel = TRUE,pt.size = 1.5)
# RENAME clusters 'cluster_label'
Idents(sce) <- "celltype"
table(Idents(sce))
new.cluster.ids <- c("TBC","BC","MC","Oc90","Interdental","SC/OL",
                     "DC_OPC","Ube2c+","IPh_IBC","HC","CEC","IC"
)
names(new.cluster.ids) <- levels(sce)
sce<- RenameIdents(sce, new.cluster.ids)
DimPlot(sce, reduction = "umap",label = TRUE,repel = TRUE)
sce$celltype <- Idents(sce)
Idents(sce) <- "celltype"
table(Idents(sce))
DimPlot(sce, reduction = "umap",label = TRUE,repel = TRUE)


table(sce2$celltype)
saveRDS(sce2, file = "human_GW15_coe_monocle.rds")

sce=human_GW15_coe_monocle
sce2=cochlea_GW15
DimPlot(sce, reduction = "umap",label = TRUE,repel = TRUE,pt.size = 1.5)

DimPlot(sce2, reduction = "umap",label = TRUE,repel = TRUE,pt.size = 1.5)


cochlea_GW15 <- readRDS("F:/PROJECTS/PROJECT_HUMAN FETAL COCHLEAE/WORKPLACE/R/Human_GW15_GW17/cochlea_GW15.rds")
sce=cochlea_GW15
DefaultAssay(sce) <- "SCT"
markers.to.plot<-c("TMC1","INSM1","OTOF","MYO7A","CABP2","ATOH1",# HC
                   "TSC22D1","EPYC","FBXO2","TMOD1","FGFR3","GATA3","GJB2",'DPYSL2',#DCOPC
                   "DGKB","SLC1A3","S100A6","OTOG","OTOGL","LGR5","LGR6",#IPHIBC
                   "EMILIN2",	"AXIN2",	"TJP1",	"ATP1A2","OTOR",#TBC
                   
                   
                   "MBP",	"TUBB4A",	"PMP22",	"PLP1","MPZ",#SC/OL
                   
                   "DCT","TYR","MET",	"KCNJ13",	"KCNJ10",#IC
                   "DCLK1","ESRRB","KCNQ1",	"KCNE1","DCN",	"ITIH5",#MC
                   "EMCN","ARHGAP6","CLDN11",#BC
                   "VWF","MECOM",	"GBP7",#CEC
                   
                   "ATP2B2","SLC26A7","MEIS1","DKK2",#FIBRO3 
                   "SMOC2","OTOA","PTGDS","CDKN1C", # Interdental
                   "UBE2C","KIAA0101", #Ube2c+ 
                   "OC90","TTR","VMO1","ECEL1","ENPEP" #OC90
)







DotPlot(sce,features=markers.to.plot)+ RotatedAxis()
Idents(sce)<- "celltype"
Idents(sce) <- factor(Idents(sce), levels = c("HC","DCOPC","IPhIBC","Oc90",
                                              "Interdental",
                                              "Ube2c+",
                                              "TBC","SC/OL","IC","MC",
                             "BC","CEC"
                                                                          ))
sce$subclasses=Idents(sce)
markers.to.plot<-c("TMC1","INSM1","OTOF","MYO7A","CABP2","ATOH1",# HC
                   "EPYC","FBXO2","TMOD1","FGFR3","GATA3","GJB2",'DPYSL2',#DCOPC
                   
                   "DGKB","SLC1A3","S100A6","OTOG","OTOGL","LGR5","LGR6",#IPHIBC
                   "OC90","TTR","VMO1","ECEL1","ENPEP", #OC90
                   "SMOC2","OTOA","PTGDS","CDKN1C", # Interdental
                   "UBE2C","KIAA0101", #Ube2c+ 
                   
                   "EMILIN2",	"AXIN2",	"TJP1",	"ATP1A2","OTOR",#TBC
                   
                   
                   "MBP",	"TUBB4A",	"PMP22",	"PLP1","MPZ",#SC/OL
                   
                   "DCT","TYR","MET",	"KCNJ13",	"KCNJ10",#IC
                   "DCLK1","ESRRB","KCNQ1",	"KCNE1","DCN",	"ITIH5",#MC
                   "EMCN","ARHGAP6","CLDN11",#BC
                   "VWF","MECOM",	"GBP7"#CEC
                   
                 
                   
                  
                   
)
DotPlot(sce, features = markers.to.plot, dot.scale = 8, group.by = "subclasses",
        cols  =c("white", "#ad9300")) + RotatedAxis()



Idents(sce) <- "subclasses"
table(Idents(sce))
new.cluster.ids <- c("HC","IPh_IBC","DC_PC_HeC","Oc90","IDC","Ube2c+",
                     "TBC","GC","IC","MC","BC","CEC"
)
names(new.cluster.ids) <- levels(sce)
sce<- RenameIdents(sce, new.cluster.ids)
DimPlot(sce, reduction = "umap",label = TRUE,repel = TRUE)
sce$subclasses <- Idents(sce)
Idents(sce) <- "subclasses"
table(Idents(sce))
DimPlot(sce, reduction = "umap",label = TRUE,repel = TRUE)

markers.to.plot<-c("TMC1","INSM1","OTOF","MYO7A","CABP2","ATOH1",# HC
                   
                   "EPYC","SLC1A3","S100A6","OTOG","OTOGL",#IPHIBC
                   
                   "PROX1","TMOD1","FGFR3","GATA3","GJB2","LGR5","LGR6",#DCOPC
                   
                   
                   "OC90","TTR","VMO1","ECEL1","ENPEP", #OC90
                   "SMOC2","OTOA","PTGDS","CDKN1C", # Interdental
                   "UBE2C","KIAA0101", #Ube2c+ 
                   
                   "EMILIN2",	"AXIN2",	"TJP1",	"ATP1A2","OTOR",#TBC
                   
                   
                   "MBP",	"TUBB4A",	"PMP22",	"PLP1","MPZ",#SC/OL
                   
                   "DCT","TYR","MET",	"KCNJ13",	"KCNJ10",#IC
                   "DCLK1","ESRRB","KCNQ1",	"KCNE1","DCN",	"ITIH5",#MC
                   "EMCN","ARHGAP6","CLDN11",#BC
                   "VWF","MECOM"#CEC
                   
                   
                   
                   
                   
)
DotPlot(sce, features = markers.to.plot, dot.scale = 8, group.by = "subclasses",
        cols  =c("white", "#ad9300")) + RotatedAxis()

sce <- RunUMAP(sce, reduction = "pca", dims = 1:20)
DimPlot(sce, reduction = "umap",label = TRUE)

cluster_cols <- c("#0868AC","#5AB2A8","#91D4CA","#F09150","#1A9850","#C7EAE5","#E47440",
                  "#D85730","#BF1D10","#A6D96A","#00C3FF","#6E4DA2","#3CA0C9","#288A82",
                  "#FDAE61","#2B8DBF","#085584","#4EB3D3","#b680a9","#60B85D","#01665E",
                  "#197AB5")#CB3A20
DimPlot(sce,reduction = "umap",label = TRUE,repel = TRUE,
        group.by = "subclasses",cols = cluster_cols)
saveRDS(sce,file = "HUMAN_GW15.rds")
##TRANSFER TO H5AD

library(sceasy)
library(reticulate)
use_condaenv('EnvironmentName')

DefaultAssay(sce) <- "RNA"
sceasy::convertFormat(sce, from="seurat", to="anndata",
                      outFile='Hu_GW15_python.h5ad')
