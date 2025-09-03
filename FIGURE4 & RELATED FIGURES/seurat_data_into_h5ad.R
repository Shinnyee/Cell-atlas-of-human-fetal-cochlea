##R version 4.2.3 (2023-03-15 ucrt)
#sceasy_0.0.7                reticulate_1.35.0 
#mouse hc data export into h5ad format
rm(list=ls())
library(sceasy)
#remotes::install_version("reticulate", version = "1.40.0")
library(reticulate)
library(Seurat)
library(rtracklayer)
library(tibble)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(Hmisc)
library(cowplot)
library(zoo)
use_condaenv('EnvironmentName')
human_sgn_filtered <- readRDS("F:/PROJECTS/PROJECT_HUMAN_FETAL_COCHLEAE/WORKPLACE/R/sc_pipeline/SGN analysis/human_sgn_filtered.rds")
human_hc_filtered <- readRDS("F:/PROJECTS/PROJECT_HUMAN_FETAL_COCHLEAE/WORKPLACE/R/sc_pipeline/HC_DEGs_analysis_GSEA_analysis/human_hc_filtered.rds")
mouse_hc_filtered <- readRDS("F:/PROJECTS/PROJECT_HUMAN_FETAL_COCHLEAE/WORKPLACE/R/sc_pipeline/HC_DEGs_analysis_GSEA_analysis/mouse_hc_filtered.rds")
sce_mouse=mouse_hc_filtered
sce_human=human_sgn_filtered
DefaultAssay(sce_mouse) <- "RNA"
sceasy::convertFormat(sce_mouse, from="seurat", to="anndata",
                      outFile='mouse_developmental_hc_python_human_id.h5ad')

DefaultAssay(sce_human) <- "RNA"
sceasy::convertFormat(sce_human, from="seurat", to="anndata",
                      outFile='human_developmental_sgn_python_human_id.h5ad')

sessionInfo()
#R version 4.2.3 (2023-03-15 ucrt)
#sceasy_0.0.7                reticulate_1.35.0