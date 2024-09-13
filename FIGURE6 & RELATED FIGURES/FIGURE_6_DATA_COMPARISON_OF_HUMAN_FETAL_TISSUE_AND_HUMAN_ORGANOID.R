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
HUMAN_cochlea_combined_filtered <- readRDS("F:/PROJECTS/PROJECT_HUMAN FETAL COCHLEAE/WORKPLACE/R/HUMAN_cochlea_combined_filtered.rds")
hu_gw23=HUMAN_cochlea_combined_filtered
DimPlot(hu_gw23)
cluster_cols <- c("#0868AC","#5AB2A8","#91D4CA","#F09150","#1A9850","#C7EAE5","#E47440",
                  "#D85730","#BF1D10","#A6D96A","#00C3FF","#6E4DA2","#3CA0C9","#288A82",
                  "#FDAE61","#2B8DBF","#085584","#4EB3D3","#b680a9","#60B85D","#01665E",
                  "#197AB5")#CB3A20
DimPlot(hu_gw23,reduction = "umap",label = TRUE,repel = TRUE,
        group.by = "subclasses",cols = cluster_cols)
Idents(hu_gw23) <- "subclasses"
hu_gw23_coe <- subset( hu_gw23,idents = c("HC","DC_PC_HeC","IPh_IBC"))
DimPlot(hu_gw23_coe,reduction = "umap",label = TRUE,repel = TRUE,
        group.by = "subclasses",cols = cluster_cols)

HUMAN_GW17 <- readRDS("F:/PROJECTS/PROJECT_HUMAN FETAL COCHLEAE/WORKPLACE/R/HUMAN_GW17.rds")
hu_gw17=HUMAN_GW17

DimPlot(hu_gw17,reduction = "umap",label = TRUE,repel = TRUE,
        group.by = "subclasses",cols = cluster_cols)

hu_gw17_coe <- subset( hu_gw17,idents = c("HC","IPh_IBC"))

HUMAN_GW15 <- readRDS("F:/PROJECTS/PROJECT_HUMAN FETAL COCHLEAE/WORKPLACE/R/HUMAN_GW15.rds")


hu_gw15=HUMAN_GW15

DimPlot(hu_gw15,reduction = "umap",label = TRUE,repel = TRUE,
        group.by = "subclasses",cols = cluster_cols)

hu_gw15_coe <- subset( hu_gw15,idents = c("HC","DC_PC_HeC","IPh_IBC"))
cochlea_GW9 <- readRDS("F:/PROJECTS/PROJECT_HUMAN FETAL COCHLEAE/WORKPLACE/R/Human_GW7_GW9/cochlea_GW9.rds")
cochlea_GW7 <- readRDS("F:/PROJECTS/PROJECT_HUMAN FETAL COCHLEAE/WORKPLACE/R/Human_GW7_GW9/cochlea_GW7.rds")


hu_gw9=cochlea_GW9
hu_gw7=cochlea_GW7


DimPlot(hu_gw9,reduction = "umap",label = TRUE,repel = TRUE,
        group.by = "subclasses",cols = cluster_cols)
hu_gw7_coe <- subset(hu_gw7,idents = "CoE_prosensory")
hu_gw9_coe <- subset(hu_gw9,idents = c("CoE_prosensory"))

hu_gw7_coe$subclasses_devtime <- paste(hu_gw7_coe$subclasses,
                                       hu_gw7_coe$orig.ident,sep="_")
table(hu_gw7_coe$subclasses_devtime)
hu_gw9_coe$subclasses_devtime <- paste(hu_gw9_coe$subclasses,
                                       hu_gw9_coe$orig.ident,sep="_")
table(hu_gw9_coe$subclasses_devtime)
hu_gw15_coe$subclasses_devtime <- paste(hu_gw15_coe$subclasses,
                                       hu_gw15_coe$orig.ident,sep="_")
table(hu_gw15_coe$subclasses_devtime)

hu_gw17_coe$subclasses_devtime <- paste(hu_gw17_coe$subclasses,
                                        hu_gw17_coe$orig.ident,sep="_")
table(hu_gw17_coe$subclasses_devtime)
hu_gw23_coe$devtime = "gw23_26"
hu_gw23_coe$subclasses_devtime <- paste(hu_gw23_coe$subclasses,
                                        hu_gw23_coe$devtime,sep="_")
table(hu_gw23_coe$subclasses_devtime)

hu_fetal_all <- merge(x=hu_gw23_coe,y=c(hu_gw17_coe,hu_gw15_coe,
                                        hu_gw9_coe,hu_gw7_coe))
hu_fetal_all

DefaultAssay(hu_fetal_all) <- "RNA"
hu_fetal_all <- hu_fetal_all %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData()
hu_fetal_all <- SCTransform(hu_fetal_all, vars.to.regress = c("nCount_RNA"),
                       method="glmGamPoi",verbose = T)

saveRDS(hu_fetal_all,file = "hu_fetal_all.rds")
hu_fetal_all <- readRDS("F:/PROJECTS/PROJECT_HUMAN FETAL COCHLEAE/WORKPLACE/R/hu_fetal_all.rds")
table(hu_fetal_all$subclasses_devtime)
Idents(hu_fetal_all) <- "subclasses_devtime"

all_human_orgn_new <- readRDS("F:/PROJECTS/PROJECT_HUMAN FETAL COCHLEAE/WORKPLACE/R/HUMAN_ORGANOIDS/all_human_orgn_new.rds")

table(all_human_orgn_new$subclasses)
hu_orgn_all <- subset(all_human_orgn_new,idents = c("HCs","LGR5_pos_transitional","LGR5_neg_transitional"))
table(hu_orgn_all$subclasses)
hu_orgn_all$subclasses_devtime <- paste(hu_orgn_all$subclasses,
                                        hu_orgn_all$devtime,sep = "_")
table(hu_orgn_all$subclasses_devtime)


## Calculate intersection between human fetal  vs human orgn DE genes for each cell type according to Jaccard similarity index 
Idents(hu_fetal_all) <- "subclasses_devtime"
Idents(hu_orgn_all) <- "subclasses_devtime"
deg_hu_fetal <- FindAllMarkers(hu_fetal_all, assay = "SCT", slot = "data",
                               test.use = "roc")

deg_orgn_fetal <- FindAllMarkers(hu_orgn_all, assay = "SCT", slot = "data",
                               test.use = "roc")

#TOP50 FOR SIMILARITY COMPARISON

deg_hu_fetal_TOP100 <- deg_hu_fetal %>% group_by(cluster) %>% top_n(n=100,wt=myAUC)
deg_orgn_fetal_TOP100 <- deg_orgn_fetal %>% group_by(cluster) %>% top_n(n=100,wt=myAUC)


library(tidyverse)
library(readxl)
library(data.table)
library (ggplot2)
#install.packages("ggpubr")
library(ggpubr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(biomaRt)

cluster1=names(table(deg_hu_fetal$cluster))
cluster2=names(table(deg_orgn_fetal$cluster))

## jaccard
df_ja=c()

for (i in cluster1) {
  ja=c()
  for (j in cluster2) {
    a=deg_hu_fetal[deg_hu_fetal$cluster==i,]
    a=a$gene
    b=deg_orgn_fetal[deg_orgn_fetal$cluster==j,]
    b=b$gene
    jaccard=length(intersect(a,b))/length(union(a,b))
    
    ja=c(ja,jaccard)
  }
  df_ja=rbind(df_ja,ja)
}

rownames(df_ja)=paste0(cluster1)
colnames(df_ja)=paste0(cluster2)
library(pheatmap)

pheatmap(df_ja,cluster_cols = T,show_rownames = T, show_colnames =T,
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
####scLearn 
library(devtools)
library(SingleCellExperiment)
library(M3Drop)
install_github("bm2-lab/scLearn")
library(scLearn)
#Data preprocessing:

# loading the reference dataset
# loading the reference dataset
data=hu_fetal_all
data=as.SingleCellExperiment(data)
rawcounts<-assays(data)[[1]]
refe_ann<-as.character(data$subclasses_devtime)
names(refe_ann)<-colnames(data)
# cell quality control and rare cell type filtered and feature selection
data_qc<-Cell_qc(rawcounts,refe_ann,species="Hs")
data_type_filtered<-Cell_type_filter(data_qc$expression_profile,data_qc$sample_information_cellType,min_cell_number = 10)
high_varGene_names <- Feature_selection_M3Drop(data_type_filtered$expression_profile)
#Model learning:

# training the model. To improve the accuracy for "unassigned" cell, you can increase "bootstrap_times", but it will takes longer time. The default value of "bootstrap_times" is 10.
scLearn_model_learning_result<-scLearn_model_learning(high_varGene_names,
                                                      data_type_filtered$expression_profile,
                                                      data_type_filtered$sample_information_cellType,
                                                      bootstrap_times=10)


#Cell assignment:
# loading the quary cell and performing cell quality control.
data2=hu_orgn_all
data2=as.SingleCellExperiment(data2)
rawcounts2<-assays(data2)[[1]]
### the true labels of this test datasets 
#query_ann<-as.character(data2$cell_type1)
#names(query_ann)<-colnames(data2)
#query_ann<-query_ann[query_ann %in% c("alpha","beta","delta","gamma")]
#rawcounts2<-rawcounts2[,names(query_ann)]
#data_qc_query<-Cell_qc(rawcounts2,query_ann,species="Hs")
### 
data_qc_query<-Cell_qc(rawcounts2,species="Hs",gene_low=50,umi_low=50)
# Assignment with trained model above. To get a less strict result for "unassigned" cells, you can decrease "diff" and "vote_rate". If you are sure that the cell type of query cells must be in the reference dataset, you can set "threshold_use" as FALSE. It means you don't want to use the thresholds learned by scLearn.
scLearn_predict_result<-scLearn_cell_assignment(scLearn_model_learning_result,
                                                data_qc_query$expression_profile,
                                                diff=0.05,threshold_use=TRUE,
                                                vote_rate=0.6)
results1=scLearn_predict_result
rownames(results1)=scLearn_predict_result$Query_cell_id

hu_orgn_all <- AddMetaData(hu_orgn_all, metadata = results1)
DimPlot(hu_orgn_all,group.by = "Predict_cell_type")
#######
####################################
# Convert Seurat object to Anndata #
####################################

# Load Seurat object
# Install packages not yet installed
packages <- c("Seurat", "argparse")
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

# SeuratDisk (not currently available on CRAN)
if ("SeuratDisk" %in% rownames(installed.packages())) {
  library("SeuratDisk")
} else {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }
  remotes::install_github("mojaveazure/seurat-disk")
  library("SeuratDisk")
}

# Save intermediate format (h5Seurat)
SaveH5Seurat(hu_orgn_all, filename ="hu_organ_all.h5Seurat")


# Convert and save to h5ad file 
Convert("hu_organ_all.h5Seurat", dest = "h5ad")


# Save intermediate format (h5Seurat)
SaveH5Seurat(hu_fetal_all, filename ="hu_fetal_all.h5Seurat")


# Convert and save to h5ad file 
Convert("hu_fetal_all.h5Seurat", dest = "h5ad")



#
label_transfer <- read.csv("transferred_labels.csv",row.names = 1)
label_transfer2 <- read.csv("transferred_labels_rf.csv",row.names = 1)
hu_orgn_all <- AddMetaData(hu_orgn_all, metadata = label_transfer2)
DimPlot(hu_orgn_all,group.by = "classifier")











