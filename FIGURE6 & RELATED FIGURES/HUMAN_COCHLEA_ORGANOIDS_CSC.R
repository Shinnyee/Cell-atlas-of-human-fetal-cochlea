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
orgn_d20 <- Read10X(data.dir = "F:/PROJECTS/PROJECT_HUMAN FETAL COCHLEAE/WORKPLACE/R/HUMAN_ORGANOIDS/D20")
orgn_d80 <- Read10X(data.dir = "F:/PROJECTS/PROJECT_HUMAN FETAL COCHLEAE/WORKPLACE/R/HUMAN_ORGANOIDS/D80")
orgn_d109_1 <- Read10X(data.dir = "F:/PROJECTS/PROJECT_HUMAN FETAL COCHLEAE/WORKPLACE/R/HUMAN_ORGANOIDS/D109_1")
orgn_d109_2 <- Read10X(data.dir = "F:/PROJECTS/PROJECT_HUMAN FETAL COCHLEAE/WORKPLACE/R/HUMAN_ORGANOIDS/D109_2")

orgn_d20<- CreateSeuratObject(counts =orgn_d20)
orgn_d20$devtime <- "D20"
orgn_d20[['percent.mt']] <- PercentageFeatureSet(orgn_d20, pattern = '^MT-')


orgn_d80<- CreateSeuratObject(counts =orgn_d80)
orgn_d80$devtime <- "D80"
orgn_d80[['percent.mt']] <- PercentageFeatureSet(orgn_d80, pattern = '^MT-')

orgn_d109_1<- CreateSeuratObject(counts =orgn_d109_1)
orgn_d109_1$devtime <- "D109_1"
orgn_d109_1[['percent.mt']] <- PercentageFeatureSet(orgn_d109_1, pattern = '^MT-')

orgn_d109_2<- CreateSeuratObject(counts =orgn_d109_2)
orgn_d109_2$devtime <- "D109"
orgn_d109_2[['percent.mt']] <- PercentageFeatureSet(orgn_d109_2, pattern = '^MT-')

# REMOVE OBJECT D109_1 

VlnPlot(orgn_d20,features = c("nFeature_RNA","percent.mt"))
VlnPlot(orgn_d80,features = c("nFeature_RNA","percent.mt"))
VlnPlot(orgn_d109_2,features = c("nFeature_RNA","percent.mt"))

#QC 
orgn_d20
orgn_d20 <- subset(orgn_d20, subset = nFeature_RNA > 1500& nFeature_RNA < 6000 & percent.mt < 12.5 )


orgn_d80
orgn_d80 <- subset(orgn_d80, subset = nFeature_RNA > 2000& nFeature_RNA < 6000 & percent.mt < 12.5 )

orgn_d109_2
orgn_d109_2 <- subset(orgn_d109_2, subset = nFeature_RNA > 1500& nFeature_RNA < 6000 & percent.mt < 12.5 )


sce_all <- merge(x = orgn_d20, 
                 y = c(orgn_d80, orgn_d109_2),
                 add.cell.ids = c("D20","D80","D109"))
sce_all
saveRDS(sce_all,file = "all_human_orgn.rds")

DefaultAssay(sce_all) <- "RNA"
sce_all <- sce_all %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData()



sce_all <- RunPCA(sce_all, features = VariableFeatures(sce_all), npcs = 50)
library(harmony)
sce_all_harmony <- RunHarmony(sce_all,group.by.vars = 'devtime',reduction = "pca",
                              dims.use = 1:10,assay.use = "RNA")#,lambda = 3,theta = 0.6
#sce_all_harmony <- RunHarmony(sce_all,group.by.vars = 'orig.ident',reduction = "pca",
#dims.use = 1:20,assay.use = "RNA")
sce_all[["harmony"]] <- sce_all_harmony[["harmony"]]
sce_all <- RunUMAP(sce_all,dims = 1:10,
                   reduction = "harmony",reduction.name = "umap_harmony")

p1 <- DimPlot(sce_all, reduction = "umap_harmony", group.by = "devtime",pt.size = 1) + 
  ggtitle("UMAP Harmony")
p1
p2 <- DimPlot(sce_all, reduction = "umap_harmony",pt.size = 1,group.by = "devtime",
              split.by = "devtime",repel = TRUE) + 
  ggtitle("UMAP Harmony")
p2
p1+p2

sce_all <- SCTransform(sce_all, vars.to.regress = c("nCount_RNA","percent.mt"),
                       method="glmGamPoi",verbose = T)

saveRDS(sce_all,file = "all_human_orgn.rds")

#UMAP devtime
data_to_plot <- data.frame(sce_all@reductions$umap@cell.embeddings)

sce_all$new_id <- colnames(sce_all)
data_to_plot$devtime <- sce_all$devtime[match(rownames(data_to_plot), sce_all$new_id)]
data_to_plot$plot_order <- 0
for(i in 1:nrow(data_to_plot)){
  data_to_plot$plot_order[i] <- sample(1:100000, 1)
  print(i)
}
data_to_plot <- arrange(data_to_plot, plot_order)

p1<- ggplot(data_to_plot,
            aes(x = UMAP_1, y = UMAP_2, color = devtime)) +
  geom_point(size = 1.5)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_color_manual(values = c("#9b5792",  
                                "#74c5e9",
                                "#b7de27" ))#"#575add","#62edd0","#1db342","#e8ba69"

p1

FeaturePlot(sce_all,features = c("MYO7A","GATA3","OTOF","ATOH1"))
sce_all <- FindNeighbors(sce_all, reduction = "harmony", dims = 1:10)
DefaultAssay(sce_all) <- "RNA"
sce_all <- FindClusters(sce_all)


FeaturePlot(sce_all,features = c("MYO7A","GATA3","OTOF","ATOH1","OTOG","OTOGL"),label = T)
VlnPlot(sce_all,features = c("MYO7A","GATA3","OTOF","ATOH1","OTOG","OTOGL"))
saveRDS(sce_all,file = "all_human_orgn.rds")
all_human_orgn <- readRDS("F:/PROJECTS/PROJECT_HUMAN FETAL COCHLEAE/WORKPLACE/R/HUMAN_ORGANOIDS/all_human_orgn.rds")
sce=all_human_orgn
DefaultAssay(sce) <- "SCT"
FeaturePlot(sce,features = c("MYO7A","GATA3","OTOF","ATOH1","OTOG","OTOGL"),label = T)

markers.to.plot<-c( "NR2F1","GATA3","INSM1","HES6","FGF8","GNG8",# cochlear/ventral markers
                    "NDRG1","TEKT1","PCDH20","CD164L2","MEIS2","NEUROD6",#VESTIBULAR MARKERS
                    "ATOH1","CCER2","KCNH6","GRXCR2","MYO7A","LHX3",#HAIR CELLS
                    "COL9A2","OC90","LGR5",    #TRANSITIONAL
                    "FGFR3","PROX1",#LATERAL
                    "SHH","NTN1","SPON1",#FP NEURONS
                    "COL1A2","POSTN","TWIST1",#MESENCHYMAL
                    "FABP5","FABP7","PTGDS","NTRK2",#EARLY DRG
                    "STMN2","STMN4","DCX","RTN1",#MATURE DRG
                    "TOP2A","CDK1","UBE2C"#CYCLING CELLS

)
DotPlot(sce, features = markers.to.plot, dot.scale = 8, group.by = "seurat_clusters",
        cols  =c("white", "#ad9300")) + RotatedAxis()

DimPlot(sce,label = T)

# remove cluster 0,10,18
sce <- subset(sce,idents = c("0","10","18"),invert=T)
table(sce$seurat_clusters)
sce2 <- RunUMAP(sce,dims = 1:10,
                reduction = "harmony",reduction.name = "umap_harmony")
DimPlot(sce2,label = T)

# RENAME clusters "subclasses"
Idents(sce2) <- "seurat_clusters"
table(Idents(sce2))
new.cluster.ids <- c("Early_DRG","Mature_DRG","HCs","LGR5_neg_transitional",
                     "LGR5_neg_transitional","LGR5_pos_transitional",
                     "Early_DRG","Cycling_cells","LGR5_pos_transitional",
                     "FP_Neurons","Cycling_cells","Mesenchymal",
                     "Early_ventral_cells_1","Early_ventral_cells_2",
                     "Mature_DRG","Early_DRG"
                     
)
names(new.cluster.ids) <- levels(sce2)
sce2<- RenameIdents(sce2, new.cluster.ids)
DimPlot(sce2, reduction = "umap_harmony",label = TRUE,repel = TRUE)
sce2$subclasses <- Idents(sce2)
Idents(sce2) <- "subclasses"
table(Idents(sce2))
DimPlot(sce2, reduction = "umap_harmony",label = TRUE,repel = TRUE)
cluster_cols <- c("#0868AC","#5AB2A8","#91D4CA","#F09150","#1A9850","#C7EAE5","#E47440",
                  "#D85730","#BF1D10","#A6D96A","#00C3FF","#6E4DA2","#3CA0C9","#288A82",
                  "#FDAE61","#2B8DBF","#085584","#4EB3D3","#b680a9","#60B85D","#01665E",
                  "#197AB5")#CB3A20
DimPlot(sce2,reduction = "umap_harmony",label = TRUE,repel = TRUE,
        group.by = "subclasses",cols = cluster_cols)



#UMAP devtime
data_to_plot <- data.frame(sce2@reductions$umap_harmony@cell.embeddings)

sce2$new_id <- colnames(sce2)
data_to_plot$devtime <- sce2$devtime[match(rownames(data_to_plot), sce2$new_id)]
data_to_plot$plot_order <- 0
for(i in 1:nrow(data_to_plot)){
  data_to_plot$plot_order[i] <- sample(1:100000, 1)
  print(i)
}
data_to_plot <- arrange(data_to_plot, plot_order)

p1<- ggplot(data_to_plot,
            aes(x = UMAP_1, y = UMAP_2, color = devtime)) +
  geom_point(size = 1.5)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_color_manual(values = c("#9b5792",  
                                "#74c5e9",
                                "#b7de27" ))#"#575add","#62edd0","#1db342","#e8ba69"

p1

Idents(sce2) <- factor(Idents(sce2), levels = c("HCs","LGR5_pos_transitional","LGR5_neg_transitional",
                                                
                                                "Early_ventral_cells_1","Early_ventral_cells_2",
                                                "Early_DRG","Mature_DRG","FP_Neurons",
                                                "Mesenchymal",
                                                "Cycling_cells"
  
))
sce2$subclasses=Idents(sce2)

table(sce2$subclasses)
all_human_orgn_new <- readRDS("F:/PROJECTS/PROJECT_HUMAN FETAL COCHLEAE/WORKPLACE/R/HUMAN_ORGANOIDS/all_human_orgn_new.rds")
sce2=all_human_orgn_new
DefaultAssay(sce2) <- "SCT"
markers.to.plot<-c( 
                    "ATOH1","CCER2","KCNH6","GRXCR2","MYO7A","LHX3",#HAIR CELLS
                    "RORB",
                    "COL9A2","OC90","LGR5",    #TRANSITIONAL
                    "FGFR3","PROX1",#LATERAL
                    "NR2F1","GATA3","INSM1","HES6","FGF8","GNG8",# cochlear/ventral markers
                    "TEKT1","PCDH20","CD164L2","MEIS2","NEUROD6",#VESTIBULAR MARKERS
                    "FABP5","FABP7","PTGDS","NTRK2",#EARLY DRG
                    "STMN2","STMN4","DCX","RTN1",#MATURE DRG
                    
                    "SHH","NTN1","SPON1",#FP NEURONS
                    "TUBB3","PRPH","CALB1","ESRRG","PBX3",
                    "COL1A2","POSTN","TWIST1",#MESENCHYMAL
                   
                    "TOP2A","CDK1","UBE2C"#CYCLING CELLS
                    
)
DotPlot(sce2, features = markers.to.plot, dot.scale = 8, group.by = "subclasses",
        cols  =c("white", "#ad9300")) + RotatedAxis()





#MONOCLE3
library(monocle3)
sce_monocle3=sce2

Idents(sce_monocle3) <- "devtime"
table(Idents(sce_monocle3))
data <- GetAssayData(sce_monocle3,assay = 'SCT',slot='counts')
cell_metadata <- sce_monocle3@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 30)
plot_pc_variance_explained(cds)
cds
cds <- reduce_dimension(cds,preprocess_method = 'PCA',reduction_method = 'UMAP')
plot_cells(cds)
colnames(colData(cds))
p1 <- plot_cells(cds,reduction_method = 'UMAP',
                 color_cells_by = 'devtime')+
  ggtitle('cds.umap')
p1  
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(sce_monocle3, reduction = 'umap_harmony')
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP<- int.embed
p2 <- plot_cells(cds,reduction_method = 'UMAP',
                 color_cells_by = 'devtime')+
  ggtitle('sce.umap')
p2
p= p1|p2
p

cds <- cluster_cells(cds)
plot_cells(cds,color_cells_by = 'partition',reduction_method = 'UMAP')
cds <- learn_graph(cds)
p= plot_cells(cds,color_cells_by = 'subclasses',label_groups_by_cluster = FALSE,
              label_leaves = FALSE,label_branch_points = FALSE)
p
plot_cells(cds,color_cells_by = 'devtime',label_groups_by_cluster = FALSE,
           label_leaves = TRUE,label_branch_points = TRUE,graph_label_size = 2,
           label_principal_points = TRUE)
cds <- order_cells(cds)
#cds_sub <- choose_graph_segments(cds)
p=plot_cells(cds,color_cells_by = 'pseudotime',label_cell_groups = FALSE,
             label_leaves = FALSE,label_branch_points = FALSE,cell_size = 0.8)
p




p3=plot_cells(cds,
              color_cells_by = "pseudotime",
              label_cell_groups=FALSE,
              label_leaves=FALSE,
              label_branch_points=FALSE,
              graph_label_size=1.5,cell_size = 0.8)
p3
saveRDS(cds,file = "human_orgn_MONOCLE3.rds")
saveRDS(sce2,file = "all_human_orgn_new.rds")


###TRANSFER TO H5AD

library(sceasy)
library(reticulate)
use_condaenv('EnvironmentName')

DefaultAssay(sce2) <- "RNA"
sceasy::convertFormat(sce2, from="seurat", to="anndata",
                      outFile='Hu_organoids_python.h5ad')




