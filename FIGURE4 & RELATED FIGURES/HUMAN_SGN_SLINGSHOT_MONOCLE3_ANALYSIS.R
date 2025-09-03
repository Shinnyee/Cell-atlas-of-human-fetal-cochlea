#############################################################################
#############################################################################
######======================#monocle3========================##########
#############################################################################
library(Seurat)
library(dplyr)
library(patchwork)
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
library(Rcpp)
library(harmony)
library(monocle3)
library(tidyverse)
library(dittoSeq)
#pre-run slingshot analysis pipeline
DimPlot(sce,
        cols = my_cols, label=FALSE , repel=TRUE,reduction = "umap",pt.size = 2
        ,group.by = "cell_type_final_v2")#
data <- GetAssayData(sce29,assay = 'RNA',slot='counts')
cell_metadata <- sce29@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 50)
plot_pc_variance_explained(cds)
cds <- reduce_dimension(cds,preprocess_method = 'PCA',reduction_method = 'UMAP')
plot_cells(cds)
colnames(colData(cds))
p1 <- plot_cells(cds,reduction_method = 'UMAP',
                 color_cells_by = 'HCtype')+
  ggtitle('cds.tsne')
p1  
cds.embed <- cds@int_colData$reducedDims$UMAP
DimPlot(sce29, reduction = "umap_harmony",label = TRUE)
int.embed <- Embeddings(sce29, reduction = 'umap_harmony')
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP<- int.embed
p2 <- plot_cells(cds,reduction_method = 'UMAP',
                 color_cells_by = 'HCtype')+
  ggtitle('sce.tsne')
p2
p= p1|p2
p
genes_of_interest <- c("Slc26a5","Otof","Tmc1")
plot_cells(cds, genes = genes_of_interest)
cds <- cluster_cells(cds)
plot_cells(cds,color_cells_by = 'partition',reduction_method = 'UMAP')
cds <- learn_graph(cds)
p= plot_cells(cds,color_cells_by = 'cluster',label_groups_by_cluster = FALSE,
              label_leaves = FALSE,label_branch_points = FALSE)
p
plot_cells(cds,color_cells_by = 'hctype',label_groups_by_cluster = FALSE,
           label_leaves = TRUE,label_branch_points = TRUE,graph_label_size = 2,
           label_principal_points = TRUE)
cds <- order_cells(cds)

p=plot_cells(cds,color_cells_by = 'pseudotime',label_cell_groups = FALSE,
             label_leaves = FALSE,label_branch_points = FALSE,cell_size = 1.25)
p

plot_genes_in_pseudotime(cds[c("Zfp462","Kdm5b"),],color_cells_by = "pseudotime",
                         ncol=3,cell_size = 2) 
###############################dotplot######################################################
Alltype_markers <- c("Pcp4","Lmo7","Dnm3","Tmc1","Slc26a5",
                     "Slc17a8","Cabp2","Osbpl9","Nf2","Strip2","Otof","Calb2",
                     "Kcnj13","Atp2a3","Adamts13","Shtn1","Cacng2","Tbx2","Fgf8","Tnfaip1",
                     "Atoh1","Jag2","Sox2","Insm1","Cdh1"
)
Alltype_markers2 <- c("Cdh1","Insm1","Sox2","Jag2","Atoh1","Tnfaip1","Fgf8","Tbx2","Cacng2",
                      "Shtn1","Adamts13","Atp2a3","Kcnj13","Calb2","Otof","Strip2","Nf2","Osbpl9","Cabp2",
                      "Slc17a8","Slc26a5","Tmc1","Dnm3","Lmo7","Pcp4"
)
DotPlot(sce29_copy,features = Alltype_markers,group.by = "hctype",
        cluster.idents = F)+
  RotatedAxis()
dittoDotPlot(sce29_copy,assay = "SCT",vars = Alltype_markers2,group.by = "hctype",size = 6,
             max.color ="#f39200" )+
  RotatedAxis()
dittoDotPlot(sce29_copy,assay = "RNA",vars = Alltype_markers,group.by = "hctype",size = 6)+
  RotatedAxis()

Msce=sce29
table(Msce$hctype)
Idents(Msce)<- "hctype"
Msce <- SCTransform(Msce, method = "glmGamPoi")
#cellmarker <- FindAllMarkers(Msce, assay = "SCT", slot = "data", test.use = "roc" )
#write.csv(cellmarker,file="cellmarkers_for_sce29.csv")
