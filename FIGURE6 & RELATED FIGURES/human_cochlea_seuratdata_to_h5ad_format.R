rm(list=ls())
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
HUMAN_cochlea_combined_filtered <- readRDS("F:/PROJECTS/PROJECT_HUMAN FETAL COCHLEAE/WORKPLACE/R/HUMAN_cochlea_combined_filtered.rds")
cluster_cols <- c("#0868AC","#5AB2A8","#91D4CA","#F09150","#1A9850","#C7EAE5","#E47440",
                  "#D85730","#BF1D10","#A6D96A","#00C3FF","#6E4DA2","#3CA0C9","#288A82",
                  "#FDAE61","#2B8DBF","#085584","#4EB3D3","#b680a9","#60B85D","#01665E",
                  "#197AB5")#CB3A20
DimPlot(HUMAN_cochlea_combined_filtered,reduction = "umap",label = TRUE,repel = TRUE,
        group.by = "subclasses",cols = cluster_cols)

DefaultAssay(HUMAN_cochlea_combined_filtered) <- "SCT"
table(Idents(HUMAN_cochlea_combined_filtered))
Idents(HUMAN_cochlea_combined_filtered) <- factor(Idents(HUMAN_cochlea_combined_filtered), levels = c("HC","DC_PC_HeC","IPh_IBC","TBC","SGN","GC","IC","MC",
                                                                                  "BC","PVM/M","CEC","Fib","Fibro1","Fibro2","Fibro3","Fibro4","SMC",
                                                                                  "Immun","ECs","EBs","NPCs"))
markers.to.plot<-c("TMC1","PCP4","GJB2","GATA3","TSC22D1","EPYC","OTOG","OTOGL","EMILIN2","RARRES1","NEFL","SNAP25","MBP","PLP1",
                   "MPZ","PMP22","DCT","TYR","DCLK1","ESRRB","TJP1","CLDN11", "ITGA8",	"CREB5",   "VWF","MECOM",
                   "VEPH1","CHST9","COCH","SLC4A10","SLC7A11","SLC13A3","ATP2B2","SLC26A7",
                   "CEMIP","THSD4","RUNX2","COL1A2","CD163","MEF2C","IGFBP7","COL4A1",
                   "HBG2","HBM","TOP2A","HMGB2")
DotPlot(HUMAN_cochlea_combined_filtered, features = markers.to.plot, dot.scale = 8,
        cols  =c("white", "#ad9300")) + RotatedAxis()

MOUSE_cochlea_combined_filtered_pnas <- readRDS("F:/PROJECTS/PROJECT_HUMAN FETAL COCHLEAE/WORKPLACE/R/MOUSE_cochlea_combined_filtered_pnas.rds")


DimPlot(MOUSE_cochlea_combined_filtered_pnas,reduction = "umap",label = TRUE,repel = TRUE,
        group.by = "subclasses",cols = cluster_cols)

library(sceasy)
library(reticulate)
use_condaenv('EnvironmentName')
sce_human=HUMAN_cochlea_combined_filtered
DefaultAssay(sce_human) <- "RNA"
sceasy::convertFormat(sce_human, from="seurat", to="anndata",
                      outFile='Hu_cochlea_python.h5ad')

sce_mouse=MOUSE_cochlea_combined_filtered_pnas
DefaultAssay(sce_mouse) <- "RNA"
sceasy::convertFormat(sce_mouse, from="seurat", to="anndata",
                      outFile='Mu_cochlea_python.h5ad')

markers.to.plot<-c("TMC1","SLC17A8","GJB2","GATA3","TSC22D1","EPYC","OTOG","OTOGL","EMILIN2","RARRES1","NEFL","SNAP25","MBP","PLP1",
                   "MPZ","PMP22","DCT","TYR","DCLK1","ESRRB","TJP1","CLDN11", "ITGA8",	"CREB5",   "VWF","MECOM",
                   "VEPH1","CHST9","COCH","SLC4A10","SLC7A11","SLC13A3","ATP2B2","SLC26A7",
                   "CEMIP","THSD4","RUNX2","COL1A2","CD163","MEF2C","IGFBP7","COL4A1",
                   "HBG2","HBM","TOP2A","HMGB2")
DotPlot(MOUSE_cochlea_combined_filtered_pnas, features = markers.to.plot, dot.scale = 8,
        cols  =c("white", "#ad9300")) + RotatedAxis()

