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
#########################################################################################
#substract sgn from HUMAN FETAL cochlea
seurat_human_new <- readRDS("C:/Users/Dell/Desktop/WORKPLACE/R/seurat_human_new.rds")
Idents(seurat_human_new) <- "subclasses_label"
table(seurat_human_new$subclasses_label)
DimPlot(seurat_human_new)
sgn_human <- subset(seurat_human_new, idents="SGN")
sgn_human
rm(seurat_human_new)
table(sgn_human$orig.ident, sgn_human$subclasses_label)
# re-cluster sgn by joint CCA-embedding
DefaultAssay(sgn_human) <- "RNA"
sgn_human$species <- "human"
combined.list <- SplitObject(sgn_human, split.by = "orig.ident")
for (i in 1:length(combined.list)) {
  combined.list[[i]] <- SCTransform(combined.list[[i]], vars.to.regress = "nCount_RNA",
                                    verbose = TRUE, method="glmGamPoi")
}
features <- SelectIntegrationFeatures(object.list =combined.list, nfeatures = 2000)
combined.list <- PrepSCTIntegration(object.list =combined.list, 
                                    anchor.features = features)
#combined.list <- lapply(X = combined.list, FUN = RunPCA, features = features)
sgn_human.anchors <- FindIntegrationAnchors(object.list = combined.list,
                                              normalization.method = "SCT", 
                                              anchor.features = features,dims = 1:10,
                                              reduction="cca", k.filter =20,
                                            k.anchor = 5
                                            ,k.score = 5)

#cochlea.anchors <- FindIntegrationAnchors(object.list = combined.list, normalization.method = "SCT", anchor.features = features, reduction = "rpca", k.anchor = 20)
sgn_human.combined.sct <- IntegrateData(anchorset = sgn_human.anchors, 
                                          normalization.method = "SCT",
                                        dims = 1:10,k.weight = 20)
sgn_human.combined.sct <- RunPCA(sgn_human.combined.sct, 
                                 features = features, npcs=30)
ElbowPlot(sgn_human.combined.sct , ndims = 30)

sgn_human.combined.sct <- FindNeighbors(sgn_human.combined.sct, 
                                          reduction = "pca", dims = 1:30,nn.eps = 0)


sgn_human.combined.sct <- FindClusters(sgn_human.combined.sct)
saveRDS(sgn_human.combined.sct,file="sgn_human_anchors_integrated_new.rds")
#sgn_human.combined.sct=sgn_human_anchors_integrated
sgn_human.combined.sct <- FindNeighbors(sgn_human.combined.sct, reduction = "pca",
                                          dims = 1:30)
sgn_human.combined.sct <- FindClusters(sgn_human.combined.sct,resolution = 1.2)
sgn_human.combined.sct <- RunUMAP(sgn_human.combined.sct, reduction = "pca", 
                                    dims = 1:30)

sgn_human.combined.sct$seurat_clusters.new <- as.integer(sgn_human.combined.sct$seurat_clusters)


# Visualization
Idents(sgn_human.combined.sct) <- "seurat_clusters"
p1 <- DimPlot(sgn_human.combined.sct, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(sgn_human.combined.sct, reduction = "umap", label = TRUE)
p3 <- DimPlot(sgn_human.combined.sct, reduction = "umap", label = TRUE,group.by = "cluster_label")
p1+p2+p3
table(sgn_human.combined.sct$orig.ident,sgn_human.combined.sct$seurat_clusters.new)
table(sgn_human.combined.sct$cluster_label,sgn_human.combined.sct$seurat_clusters.new)

#########################################################################################
DefaultAssay(sgn_human.combined.sct) <- "SCT"
markers.to.plot<-c("NEFL","SNAP25","TUBB3","NEFH","NEFM", "MAP2", #TYPE I SGN 
                   "CALB2","MDGA1","B3GAT1",  #IA
                   "SEMA3E", "CALB1","NTNG1",    #IB
                   "RUNX1", "POU4F1","LYPD1","GRM8",            #IC
                   "GATA3","SMAD6","ATP2B4", "ANO2" ,"PRPH",    #TYPE II
                   "MBP","MPZ","PMP22","SOX10","SOX9"
                   
                   
                   
                   
                   
)

VlnPlot(sgn_human.combined.sct, features = markers.to.plot)
VlnPlot(sgn_human.combined.sct, features = "PBX3")
DotPlot(sgn_human.combined.sct, features = markers.to.plot, dot.scale = 8,
        
        cols  =c("white", "#ad9300")) + RotatedAxis()
FeaturePlot(sgn_human.combined.sct,pt.size = 2,
            features = c("TUBB3","SNAP25","NEFH","NEFM","NEFL","MBP"),label = TRUE)
dev.off()
###############check more type1 AND 2 marker, based on the assumption that these stages of SGN
###appear likely to P14 and more mature

markers.to.plot<-c("NEFL","SNAP25","TUBB3","NEFH","NEFM", "MAP2", 
                   
                   #pan-SGN 
                  
                  
                   "SYT2","RYR2","EPHA4","PROX1","KCNIP1","SCN4B",
                   
                   
                   #TYPE I SGN
                   "GATA3","ATP2B4" ,"PDE1C","SORCS3",
                   "MAFB","NGFR","SNCA","FXYD6","EFNA5",
                  "RPH3A","TLE4","PRPH",

                   
                   #TYPE II
                  "MBP","MPZ","PMP22","SOX10","SOX9","PLP1"
                  #GLIAL CELLS
                   
)

VlnPlot(sgn_human.combined.sct, features = markers.to.plot)
DotPlot(sgn_human.combined.sct, features = markers.to.plot, dot.scale = 8,
        
        cols  =c("white", "#ad9300")) + RotatedAxis()

##############remove cluster 0,1,2,4 as containing glial cells



Idents(sgn_human.combined.sct) <- "seurat_clusters"
table(sgn_human.combined.sct$seurat_clusters)
sgn_human_new <- subset(sgn_human.combined.sct, idents=c("3","6","5"))
sgn_human_new
table(sgn_human_new$orig.ident, sgn_human_new$seurat_clusters)
# re-cluster sgn by joint CCA-embedding
DefaultAssay(sgn_human_new) <- "RNA"
sgn_human_new$species <- "human"
combined.list <- SplitObject(sgn_human_new, split.by = "orig.ident")
for (i in 1:length(combined.list)) {
  combined.list[[i]] <- SCTransform(combined.list[[i]], vars.to.regress = "nCount_RNA",
                                    verbose = TRUE, method="glmGamPoi")
}
features <- SelectIntegrationFeatures(object.list =combined.list, nfeatures = 2000)
combined.list <- PrepSCTIntegration(object.list =combined.list, 
                                    anchor.features = features)
#combined.list <- lapply(X = combined.list, FUN = RunPCA, features = features)
sgn_human.anchors_new <- FindIntegrationAnchors(object.list = combined.list,
                                            normalization.method = "SCT", 
                                            anchor.features = features,dims = 1:8,
                                            reduction="cca", k.filter =8,
                                            k.anchor = 5
                                            ,k.score = 5)

#cochlea.anchors <- FindIntegrationAnchors(object.list = combined.list, normalization.method = "SCT", anchor.features = features, reduction = "rpca", k.anchor = 20)
sgn_human.combined.sct <- IntegrateData(anchorset = sgn_human.anchors_new, 
                                        normalization.method = "SCT",
                                        dims = 1:8,k.weight = 5)
sgn_human.combined.sct <- RunPCA(sgn_human.combined.sct, 
                                 features = features, npcs=30)
ElbowPlot(sgn_human.combined.sct , ndims = 30)

sgn_human.combined.sct <- FindNeighbors(sgn_human.combined.sct, 
                                        reduction = "pca", dims = 1:30,nn.eps = 0)


sgn_human.combined.sct <- FindClusters(sgn_human.combined.sct)
saveRDS(sgn_human.combined.sct,file="sgn_human_anchors_integrated_new.rds")
sgn_human.combined.sct=sgn_human_anchors_integrated_new
DefaultAssay(sgn_human.combined.sct) <- "integrated"
sgn_human.combined.sct <- FindNeighbors(sgn_human.combined.sct, reduction = "pca",
                                        dims = 1:30)
sgn_human.combined.sct <- FindClusters(sgn_human.combined.sct)
sgn_human.combined.sct <- RunUMAP(sgn_human.combined.sct, reduction = "pca", 
                                  dims = 1:10)

sgn_human.combined.sct$seurat_clusters.new <- as.integer(sgn_human.combined.sct$seurat_clusters)


# Visualization
Idents(sgn_human.combined.sct) <- "seurat_clusters"
p1 <- DimPlot(sgn_human.combined.sct, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(sgn_human.combined.sct, reduction = "umap", label = TRUE)
p3 <- DimPlot(sgn_human.combined.sct, reduction = "umap", label = TRUE,group.by = "cluster_label")
p1+p2+p3
table(sgn_human.combined.sct$orig.ident,sgn_human.combined.sct$seurat_clusters.new)
table(sgn_human.combined.sct$cluster_label,sgn_human.combined.sct$seurat_clusters.new)
saveRDS(sgn_human.combined.sct,file="sgn_human_new.rds")
sgn_human_new <- readRDS("F:/PROJECTS/PROJECT_HUMAN FETAL COCHLEAE/WORKPLACE/R/sgn_human_new.rds")
sgn_human.combined.sct=sgn_human_new
DefaultAssay(sgn_human.combined.sct) <- "SCT"
markers.to.plot<-c("NEFL","SNAP25","TUBB3","NEFH","NEFM", "MAP2", 
                   
                   #pan-SGN 
                   
                   
                   "SYT2","RYR2","EPHA4","PROX1","KCNIP1","SCN4B",
                   "CALB1","LYPD1","RUNX1",
                   
                   
                   #TYPE I SGN
                   "GATA3","ATP2B4" ,"PDE1C","SORCS3",
                   "MAFB","NGFR","SNCA","FXYD6","EFNA5",
                   "RPH3A","TLE4","PRPH"
                   
                   
                   #TYPE II
                   
                   
)

VlnPlot(sgn_human.combined.sct, features = markers.to.plot)
dev.off()
DotPlot(sgn_human.combined.sct, features = markers.to.plot, dot.scale = 8,
        
        cols  =c("white", "#ad9300")) + RotatedAxis()

VlnPlot(sgn_human.combined.sct, features = markers.to.plot,stack = TRUE,flip = TRUE)



markers.to.plot<-c("NEFL","SNAP25","TUBB3","NEFH","NEFM", "MAP2", 
                   
                   #pan-SGN 
                   
                   "SYT2","RYR2","EPHA4","PROX1","KCNIP1","SCN4B","ESRRG", "KCNA1", # PAN-TYPE1
                   "MDGA1","B3GAT1", "CACNA1B","PALMD", "CD24", #IA "CALB2",
                   "SEMA3E", "CALB1","NTNG1","MGAT4C","PCP4L1",        #IB
                   "RUNX1", "POU4F1","LYPD1","GRM8",         #IC
                   
                   #TYPE I SGN
                   "GATA3","ATP2B4" ,"PDE1C","SORCS3",
                   "MAFB","NGFR","SNCA","FXYD6","EFNA5",
                   "RPH3A","TLE4","PRPH","SMAD6","ANO2"
                   
                   
                   #TYPE II "TH"
)
DotPlot(sgn_human.combined.sct, features = markers.to.plot, dot.scale = 8,
        
        cols  =c("white", "#ad9300")) + RotatedAxis()
DefaultAssay(sgn_human.combined.sct) <- "integrated"
sgn_human.combined.sct <- FindNeighbors(sgn_human.combined.sct, reduction = "pca",
                                        dims = 1:30)
sgn_human.combined.sct <- FindClusters(sgn_human.combined.sct,resolution = 1.5,
                                       future.seed=TRUE)
sgn_human.combined.sct <- RunUMAP(sgn_human.combined.sct, reduction = "pca", 
                                  dims = 1:10)
library(paletteer)
pal <- paletteer_d("ggsci::nrc_npg")[c(5,3,1,4,8)] 
Idents(sgn_human.combined.sct) <- "seurat_clusters"
DimPlot(sgn_human.combined.sct, reduction = "umap",label = TRUE,repel = TRUE,
        pt.size = 3,cols = pal)
pal1 <- paletteer_d("ggsci::nrc_npg")[c(9,7,6)]
pal2 <- paletteer_d("ggsci::nrc_npg")[c(5,3,4,2,1)]
p1 <- DimPlot(sgn_human.combined.sct, reduction = "umap", cols = pal1,
              group.by = "orig.ident",pt.size = 2)
p2 <- DimPlot(sgn_human.combined.sct, reduction = "umap", cols = pal2,
              label = TRUE,pt.size = 2)
p1+p2

DefaultAssay(sgn_human.combined.sct) <- "SCT"
markers.to.plot<-c("NEFL","SNAP25","TUBB3","NEFH","NEFM", "MAP2", #pan-SGN 

                   
                   #TYPE I SGN
                  
                   "SYT2","RYR2","EPHA4","PROX1","KCNIP1","SCN4B","ESRRG", "KCNA1", # PAN-TYPE1
                   "MDGA1","B3GAT1", "CACNA1B","PALMD", "CD24", #IA "CALB2",
                   "SEMA3E", "CALB1","NTNG1","MGAT4C","PCP4L1",        #IB
                   "RUNX1", "POU4F1","LYPD1","GRM8" ,       #IC
                   
                   "GATA3","ATP2B4" ,"PDE1C","SORCS3",
                   "MAFB","NGFR","SNCA","FXYD6","EFNA5",
                   "RPH3A","TLE4","PRPH","SMAD6","ANO2"
                   
                   
                   #TYPE II "TH"
)

DotPlot(sgn_human.combined.sct, features = markers.to.plot, dot.scale = 8,
        
        cols  =c("white", "#2141B5")) + RotatedAxis()

Idents(sgn_human.combined.sct) <- factor(Idents(sgn_human.combined.sct), 
                                         levels = c("0","3","4","1","2"
                                              ))

table(Idents(sgn_human.combined.sct))
sgn_human.combined.sct$seurat_clusters_2=Idents(sgn_human.combined.sct)
DotPlot(sgn_human.combined.sct, features = markers.to.plot, dot.scale = 8,
        group.by = "seurat_clusters_2",
        
        cols  =c("white", "#2141B5")) + RotatedAxis()

DimPlot(sgn_human.combined.sct)
##name subclass


Idents(sgn_human.combined.sct) <- "seurat_clusters"
table(Idents(sgn_human.combined.sct))
new.cluster.ids <- c("Type_I", "Type_II","Type_II", "Type_I", "Type_I")
names(new.cluster.ids) <- levels(sgn_human.combined.sct)
sgn_human.combined.sct<- RenameIdents(sgn_human.combined.sct, new.cluster.ids)
DimPlot(sgn_human.combined.sct, reduction = "umap",label = TRUE,repel = TRUE)
sgn_human.combined.sct$celltype <- Idents(sgn_human.combined.sct)
Idents(sgn_human.combined.sct) <- "celltype"
table(Idents(sgn_human.combined.sct))

##NAME SUBTYPE

Idents(sgn_human.combined.sct) <- "seurat_clusters"
table(Idents(sgn_human.combined.sct))
new.cluster.ids <- c("Type_IA/B", "Type_II_intermediate","Type_II", "Type_IA/B", "Type_IC")
names(new.cluster.ids) <- levels(sgn_human.combined.sct)
sgn_human.combined.sct<- RenameIdents(sgn_human.combined.sct, new.cluster.ids)
DimPlot(sgn_human.combined.sct, reduction = "umap",label = TRUE,repel = TRUE)
sgn_human.combined.sct$subtype <- Idents(sgn_human.combined.sct)
Idents(sgn_human.combined.sct) <- "subtype"
table(Idents(sgn_human.combined.sct))

Idents(sgn_human.combined.sct) <- factor(Idents(sgn_human.combined.sct), 
                                         levels = c("Type_IA/B","Type_IC","4",
                                                    "Type_II_intermediate","Type_II"
                                         ))

table(Idents(sgn_human.combined.sct))
sgn_human.combined.sct$subtype_2=Idents(sgn_human.combined.sct)
DotPlot(sgn_human.combined.sct, features = markers.to.plot, dot.scale = 8,
        group.by = "subtype_2",
        
        cols  =c("white", "#2141B5")) + RotatedAxis()

DimPlot(sgn_human.combined.sct)


pal1 <- paletteer_d("ggsci::nrc_npg")[c(9,5,3)]
pal2 <- paletteer_d("ggsci::nrc_npg")[c(3,5)]
pal3 <- paletteer_d("ggsci::nrc_npg")[c(5,3,1,2)]
p1 <- DimPlot(sgn_human.combined.sct, reduction = "umap", cols = pal1,
              group.by = "orig.ident",pt.size = 2)
p2 <- DimPlot(sgn_human.combined.sct, reduction = "umap", cols = pal2,
              label = TRUE,pt.size = 2,group.by = "celltype")
p3 <- DimPlot(sgn_human.combined.sct, reduction = "umap", cols = pal3,
              label = TRUE,pt.size = 2,group.by = "subtype")
p1+p2+p3






library(dittoSeq)
dittoDimPlot(sgn_human.combined.sct,"PRPH",size = 4,
             min.color = "#E0FFFF",max.color = "#C71585",order = "increasing",assay = "SCT")

dittoDimPlot(sgn_human.combined.sct,"PBX3",size = 4,
             min.color = "#E0FFFF",max.color = "#C71585",assay = "SCT")



p1=multi_dittoDimPlot(sgn_human.combined.sct,c("SYT2","PROX1","PRPH","GATA3","PALMD","CACNA1B",
                                               "RUNX1","POU4F1","LYPD1","CALB1"),size = 2,
                      min.color = "#E0FFFF",max.color = "#C71585",order = "increasing",
                      nrow = 2,ncol = 5
)
p1

saveRDS(sgn_human.combined.sct,file = "human_sgn_annotation.rds")
#####################################################################################
######################################################################################
#########################################################################################
human_sgn_annotation <- readRDS("F:/PROJECTS/PROJECT_HUMAN FETAL COCHLEAE/WORKPLACE/R/human_sgn_annotation.rds")
sgn_human.combined.sct=human_sgn_annotation
library(dittoSeq)
p1=multi_dittoDimPlot(sgn_human.combined.sct,c("PBX3","GATA3","PRPH","CALB1"),size = 2,
                      min.color = "#E0FFFF",max.color = "#C71585",order = "increasing",
                      nrow = 2,ncol = 2
)
p1
p2=multi_dittoDimPlot(sgn_human.combined.sct,c("PBX3","PROX1","PRPH","TUBB3"),size = 2,
                      min.color = "#E0FFFF",max.color = "#C71585",order = "increasing",
                      nrow = 2,ncol = 2
)
p2
################################################################################################
################################################################################################
################################################################################################
#########################
##plot CALB+ PRPH+ CELLS

p1=DimPlot(sgn_human.combined.sct,group.by = "subtype",
           cols = c("grey","grey",
                    "grey","grey"),
           pt.size = 2.5)+NoLegend()

p1


gene1 <- "PRPH"
gene2 <- "CALB1"


pos_gene=sgn_human.combined.sct@reductions$umap@cell.embeddings
pos_gene1=pos_gene[sgn_human.combined.sct@assays$SCT@data[gene1,]>0
                   &sgn_human.combined.sct@assays$SCT@data[gene2,]>0,]
pos_gene2=pos_gene[sgn_human.combined.sct@assays$SCT@data[gene1,]>0,]
pos_gene3=pos_gene[sgn_human.combined.sct@assays$SCT@data[gene2,]>0,]


p1+geom_point(data = as.data.frame(pos_gene1),
              aes(x=UMAP_1,y=UMAP_2), 
              shape = 16,size = 2, color="blue")
  
p1+geom_point(data = as.data.frame(pos_gene2),
              aes(x=UMAP_1,y=UMAP_2), 
              shape = 16,size = 2, color="red")
p1+geom_point(data = as.data.frame(pos_gene3),
              aes(x=UMAP_1,y=UMAP_2), 
              shape = 16,size = 2, color="green")



select_cells1 <- WhichCells(sgn_human.combined.sct,slot = "counts",
                            expression = CALB1 >0 & PRPH >0)
select_cells2 <- WhichCells(sgn_human.combined.sct,slot = "counts",
                            expression = CALB1 >0 )
select_cells3 <- WhichCells(sgn_human.combined.sct,slot = "counts",
                            expression = PRPH >0 )
head(select_cells1)
head(list(select_cells1))


select_cells1 <- as.matrix(select_cells1)
select_cells1 <- as.data.frame(select_cells1)

select_cells2 <- as.matrix(select_cells2)
select_cells2 <- as.data.frame(select_cells2)

select_cells3 <- as.matrix(select_cells3)
select_cells3 <- as.data.frame(select_cells3)


all.cells <- data.frame(cell = unique(c(select_cells1$V1, select_cells2$V1, 
                                        select_cells3$V1)))


all.cells$calb1_prph_POS <- as.character(match(all.cells$cell, select_cells1$V1))
all.cells$calb1_POS <- as.character(match(all.cells$cell, select_cells2$V1))
all.cells$prph_POS <- as.character(match(all.cells$cell, select_cells3$V1))


all.cells$calb1_prph_POS[which(is.na(all.cells$calb1_prph_POS))] <- FALSE
all.cells$calb1_prph_POS[which(all.cells$calb1_prph_POS != FALSE)] <- TRUE
all.cells$calb1_POS[which(is.na(all.cells$calb1_POS))] <- FALSE
all.cells$calb1_POS[which(all.cells$calb1_POS != FALSE)] <- TRUE
all.cells$prph_POS[which(is.na(all.cells$prph_POS))] <- FALSE
all.cells$prph_POS[which(all.cells$prph_POS != FALSE)] <- TRUE

all.cells$calb1_prph_POS <- as.logical(all.cells$calb1_prph_POS)
all.cells$calb1_POS <- as.logical(all.cells$calb1_POS)
all.cells$prph_POS <- as.logical(all.cells$prph_POS)

all.cells2=all.cells
library(eulerr)
plot(euler(
  all.cells2[ ,2:4]),
  quantities = list(cex = 3),
  labels = TRUE,
  main = "human sgn",
  fills = c("royalblue1", "maroon4", "sienna2")
)
cells_calb1_prph_pos <- all.cells[which(all.cells$calb1_prph_POS==TRUE),]

DimPlot(sgn_human.combined.sct,label = T, cells.highlight = cells_calb1_prph_pos$cell,
        cols.highlight = c("blue"),cols = "gray",pt.size = 1)
cells_calb1_pos_only <- all.cells[which(all.cells$calb1_prph_POS==FALSE&all.cells$calb1_POS==TRUE),]
DimPlot(sgn_human.combined.sct,label = T, cells.highlight = cells_calb1_pos_only$cell,
        cols.highlight = c("blue"),cols = "gray",pt.size = 1)
cells_prph_pos_only <- all.cells[which(all.cells$calb1_prph_POS==FALSE&all.cells$prph_POS==TRUE),]
DimPlot(sgn_human.combined.sct,label = T, cells.highlight = cells_prph_pos_only$cell,
        cols.highlight = c("blue"),cols = "gray",pt.size = 1)


cells_highlight <- list(cells_calb1_prph_pos$cell,cells_calb1_pos_only$cell,cells_prph_pos_only$cell)
DimPlot(sgn_human.combined.sct,label = T, cells.highlight = cells_highlight,
        cols.highlight = c("pink","green","blue"),cols = "gray",pt.size = 1)
dev.off()
p1=DimPlot(sgn_human.combined.sct,group.by = "subtype",
           cols = c("#fbee9c","#fbee9c",
                    "#fbee9c","#fbee9c"),
           pt.size = 2.5)+NoLegend()

p1


gene1 <- "PRPH"
gene2 <- "CALB1"


pos_gene=sgn_human.combined.sct@reductions$umap@cell.embeddings
pos_gene1=pos_gene[cells_calb1_prph_pos$cell,]
pos_gene2=pos_gene[cells_calb1_pos_only$cell,]
pos_gene3=pos_gene[cells_prph_pos_only$cell,]


p1+geom_point(data = as.data.frame(pos_gene1),
              aes(x=UMAP_1,y=UMAP_2), 
              shape = 16,size = 2, color="#4de6dd")

p1+geom_point(data = as.data.frame(pos_gene2),
              aes(x=UMAP_1,y=UMAP_2), 
              shape = 16,size = 2, color="#ed158a")
p1+geom_point(data = as.data.frame(pos_gene3),
              aes(x=UMAP_1,y=UMAP_2), 
              shape = 16,size = 2, color="#353fcb")


p1+geom_point(data = as.data.frame(pos_gene1),
              aes(x=UMAP_1,y=UMAP_2), 
              shape = 16,size = 2.5, color="#4de6dd")+
             geom_point(data = as.data.frame(pos_gene2),
            aes(x=UMAP_1,y=UMAP_2), 
            shape = 16,size = 2.5, color="#ed158a")+
            geom_point(data = as.data.frame(pos_gene3),
            aes(x=UMAP_1,y=UMAP_2), 
            shape = 16,size = 2.5, color="#353fcb")


p1=multi_dittoDimPlot(sgn_human.combined.sct,c("CALB1","PRPH","TUBB3"),size = 4,
                      min.color = "#E0FFFF",max.color = "#C71585",order = "increasing",
                      nrow = 1,ncol = 3
)
p1
################################################################################################
################################################################################################
#########load P28 MOUSE DATA
sgn_mouse_combined <- readRDS("F:/PROJECTS/PROJECT_HUMAN FETAL COCHLEAE/WORKPLACE/R/sgn_mouse_combined.rds")
sgn_mouse.combined.sct=sgn_mouse_combined
# Visualization
p1 <- DimPlot(sgn_mouse.combined.sct, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(sgn_mouse.combined.sct, reduction = "umap", label = TRUE)
p3 <- DimPlot(sgn_mouse.combined.sct, reduction = "umap", label = TRUE,group.by = "cluster_label")
p1+p2+p3

dittoDimPlot(sgn_mouse.combined.sct,"PBX3",size = 2,
             min.color = "#E0FFFF",max.color = "#C71585",order = "increasing",assay = "SCT")

Idents(sgn_mouse.combined.sct) <- "orig.ident"
table(sgn_mouse.combined.sct$orig.ident)
sgn_mouse <- subset(sgn_mouse.combined.sct, idents=c("Mouse_P28_1","Mouse_P28_2"))
sgn_mouse


# re-cluster sgn cells by joint CCA-embedding
DefaultAssay(sgn_mouse) <- "RNA"
sgn_mouse$species <- "mouse"
combined.list <- SplitObject(sgn_mouse, split.by = "orig.ident")
for (i in 1:length(combined.list)) {
  combined.list[[i]] <- SCTransform(combined.list[[i]], vars.to.regress = "nCount_RNA",
                                    verbose = TRUE, method="glmGamPoi")
}
library(scrattch.hicat)
table(sgn_mouse$cluster_label)
table(sgn_mouse$seurat_clusters)

Var.genes.mouse_p28_1 <- select_markers(combined.list$Mouse_P28_1@assays$SCT@counts, 
                                        combined.list$Mouse_P28_1$seurat_clusters, 
                                        n.markers = 100)
Var.genes.mouse_p28_1.markers <- Var.genes.mouse_p28_1$markers

Var.genes.mouse_p28_2 <- select_markers(combined.list$Mouse_P28_2@assays$SCT@counts, 
                                        combined.list$Mouse_P28_2$seurat_clusters, 
                                        n.markers = 100)
Var.genes.mouse_p28_2.markers <- Var.genes.mouse_p28_2$markers


total.Var.genes <- unique(c(
                            Var.genes.mouse_p28_1$markers,Var.genes.mouse_p28_2$markers))


total.Var.genes <- total.Var.genes[which(total.Var.genes %in% rownames(combined.list$Mouse_P28_1@assays$SCT@counts))]
total.Var.genes <- total.Var.genes[which(total.Var.genes %in% rownames(combined.list$Mouse_P28_2@assays$SCT@counts))]
write.csv(total.Var.genes,file = "features used for cca anchors sgn mouse.csv")

features=total.Var.genes

combined.list <- PrepSCTIntegration(object.list =combined.list,
                                    anchor.features = features)
sgn_mouse.anchors <- FindIntegrationAnchors(object.list = combined.list,
                                            normalization.method = "SCT", 
                                            anchor.features = features,dims = 1:30,
                                            reduction="cca")
sgn_mouse.combined.sct <- IntegrateData(anchorset = sgn_mouse.anchors, 
                                        normalization.method = "SCT",dims = 1:30)
sgn_mouse.combined.sct <- RunPCA(sgn_mouse.combined.sct, features = features, npcs=30)
ElbowPlot(sgn_mouse.combined.sct , ndims = 30)
sgn_mouse.combined.sct <- FindNeighbors(sgn_mouse.combined.sct, 
                                        reduction = "pca", dims = 1:30,nn.eps = 0)
sgn_mouse.combined.sct <- FindClusters(sgn_mouse.combined.sct)


sgn_mouse.combined.sct <- FindNeighbors(sgn_mouse.combined.sct, reduction = "pca",
                                        dims = 1:15)
sgn_mouse.combined.sct <- FindClusters(sgn_mouse.combined.sct,resolution = 0.8)
sgn_mouse.combined.sct <- RunUMAP(sgn_mouse.combined.sct, reduction = "pca", dims = 1:10)

sgn_mouse.combined.sct$seurat_clusters.new <- as.integer(sgn_mouse.combined.sct$seurat_clusters)

# Visualization
p1 <- DimPlot(sgn_mouse.combined.sct, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(sgn_mouse.combined.sct, reduction = "umap", label = TRUE)
p3 <- DimPlot(sgn_mouse.combined.sct, reduction = "umap", label = TRUE,group.by = "cluster_label")
p1+p2+p3

sgn_mouse.combined.sct <- SCTransform(sgn_mouse.combined.sct, vars.to.regress = "nCount_RNA",
                                  verbose = TRUE, method="glmGamPoi")

DefaultAssay(sgn_mouse.combined.sct) <- "SCT"
markers.to.plot<-c("NEFL","SNAP25","TUBB3","NEFH","NEFM", "MAP2", #pan-SGN 
                   
                   
                   #TYPE I SGN
                   
                   "SYT2","RYR2","EPHA4","PROX1","KCNIP1","SCN4B","ESRRG", "KCNA1", # PAN-TYPE1
                   "MDGA1","B3GAT1", "CACNA1B","PALMD", "CD24","CALB2",    #IA "CALB2",
                   "SEMA3E", "CALB1","NTNG1","MGAT4C","PCP4L1",        #IB
                   "RUNX1", "POU4F1","LYPD1","GRM8" ,       #IC
                   
                   "GATA3","ATP2B4" ,"PDE1C","SORCS3",
                   "MAFB","NGFR","SNCA","FXYD6","EFNA5",
                   "RPH3A","TLE4","PRPH","SMAD6","ANO2"
                   
                   
                   #TYPE II "TH"
)

DotPlot(sgn_mouse.combined.sct, features = markers.to.plot, dot.scale = 8,
        
        cols  =c("white", "#2141B5")) + RotatedAxis()


##name subclass


Idents(sgn_mouse.combined.sct) <- "seurat_clusters"
table(Idents(sgn_mouse.combined.sct))
new.cluster.ids <- c("Type_I", "Type_I","Type_I", "Type_II")
names(new.cluster.ids) <- levels(sgn_mouse.combined.sct)
sgn_mouse.combined.sct<- RenameIdents(sgn_mouse.combined.sct, new.cluster.ids)
DimPlot(sgn_mouse.combined.sct, reduction = "umap",label = TRUE,repel = TRUE)
sgn_mouse.combined.sct$celltype <- Idents(sgn_mouse.combined.sct)
Idents(sgn_mouse.combined.sct) <- "celltype"
table(Idents(sgn_mouse.combined.sct))

##NAME SUBTYPE

Idents(sgn_mouse.combined.sct) <- "seurat_clusters"
table(Idents(sgn_mouse.combined.sct))
new.cluster.ids <- c("Type_IA", "Type_IC","Type_IB", "Type_II")
names(new.cluster.ids) <- levels(sgn_mouse.combined.sct)
sgn_mouse.combined.sct<- RenameIdents(sgn_mouse.combined.sct, new.cluster.ids)
DimPlot(sgn_mouse.combined.sct, reduction = "umap",label = TRUE,repel = TRUE)
sgn_mouse.combined.sct$subtype <- Idents(sgn_mouse.combined.sct)
Idents(sgn_mouse.combined.sct) <- "subtype"
table(Idents(sgn_mouse.combined.sct))

sgn_mouse.combined.sct <- FindNeighbors(sgn_mouse.combined.sct, reduction = "pca",
                                        dims = 1:15)
sgn_mouse.combined.sct <- FindClusters(sgn_mouse.combined.sct,resolution = 0.8)
sgn_mouse.combined.sct <- RunUMAP(sgn_mouse.combined.sct, reduction = "pca", dims = 1:30)
DimPlot(sgn_mouse.combined.sct, reduction = "umap",label = TRUE,repel = TRUE)
saveRDS(sgn_mouse.combined.sct,file = "mouse_p28_sgn.rds")
library(dittoSeq)
dittoDimPlot(sgn_mouse.combined.sct,"PRPH",size = 4,reduction.use = "umap",
             min.color = "#E0FFFF",max.color = "#C71585",order = "increasing",
             assay = "SCT")

#########################################################################################
Idents(sgn_mouse.combined.sct) <- factor(Idents(sgn_mouse.combined.sct), 
                                         levels = c("Type_IA","Type_IB","Type_IC",
                                                    "Type_II"
                                         ))

table(Idents(sgn_mouse.combined.sct))
sgn_mouse.combined.sct$subtype_2=Idents(sgn_mouse.combined.sct)
DefaultAssay(sgn_mouse.combined.sct) <- "SCT"
markers.to.plot<-c("NEFL","SNAP25","TUBB3","NEFH", "MAP2", #pan-SGN 
                   
                   
                   #TYPE I SGN
                   
                   "SYT2","RYR2","EPHA4","PROX1","KCNIP1","SCN4B","ESRRG", "KCNA1", # PAN-TYPE1
                   "MDGA1","B3GAT1", "CACNA1B","PALMD", "CD24","CALB2",    #IA "CALB2",
                   "SEMA3E", "CALB1","NTNG1","MGAT4C","PCP4L1",        #IB
                   "RUNX1", "POU4F1","LYPD1","GRM8" ,       #IC
                   
                   "GATA3","ATP2B4" ,"PDE1C","SORCS3",
                   "MAFB","NGFR","SNCA","FXYD6","EFNA5",
                   "RPH3A","TLE4","PRPH","SMAD6","ANO2"
                   
                   
                   #TYPE II "TH"
)

DotPlot(sgn_mouse.combined.sct, features = markers.to.plot, dot.scale = 8,
        
        cols  =c("white", "#2141B5")) + RotatedAxis()

# Visualization
p1 <- DimPlot(sgn_mouse.combined.sct, reduction = "umap", 
              group.by = "orig.ident",pt.size = 2)
p2 <- DimPlot(sgn_mouse.combined.sct, reduction = "umap", label = TRUE,
              group.by = "subtype",pt.size = 2)
p2+p1

FeaturePlot(sgn_mouse.combined.sct,features = "TLE4")



############################################################################



markers.to.plot<-c(
  "NEUROD1","SYT13","PTMS","EPHA5",#e14
  "GATA3","CHD7","SEMA3A","IGFBP2","VIM",#e15
  "DCX","CDH4","P2RX3", "CELF5","MARCKSL1", #e16
  "PTPRF","NTRK3","CDKN1C","PLXNA1",      #e17
  "CADM1","SERPINH1", "SERPINE2", "NR2F2",  #e18
  "NPTN","CPLX1","LDHB","SCN1A","KCNK1",   #p3
  "IGFBP6","FAM161A","BMP5","CHD8"  #p25
  
  
)

DotPlot(sgn_human.combined.sct, features = markers.to.plot,
        dot.scale = 8,cols  =c("white", "#ad9300")) + RotatedAxis()


markers.to.plot<-c(
  "NEUROD1","MFAP4","PPP1R14A","PPP1R17","CASP3",#Early SGN
  "UTS2B","IGFBP2","VIM","ALCAM","CALB1",#Intermediate SGN
  "LYPD1","PCDH9","JUN","PCP4",#IB/C precursor
  "FTH1.1","PVALB","ALDH1A3","RGS10", "RUNX1",          #Immature IB/C
  "ETV4","F2R","GATA3","CNR1","ID2",#IA/2 precursor
  "SPOCK3","CLU","KCTD12","PRPH",#IA precursor
  "CALB2","FXYD7","CNTN5","NEFH","MYL1",#Immature IA
  "EFNA5","RPH3A","TLE4","TH" # Immature Type2
)
DotPlot(sgn_human.combined.sct, features = markers.to.plot,
        dot.scale = 8,cols  =c("white", "#ad9300")) + RotatedAxis()
VlnPlot(sgn_human.combined.sct, features = markers.to.plot)
VlnPlot(sgn_human.combined.sct, features = markers.to.plot,stack = TRUE,flip = TRUE)

dev.off()

MOUSE_SGN_HARMONY_new <- readRDS("F:/PROJECTS/PROJECT_HUMAN FETAL COCHLEAE/WORKPLACE/R/MOUSE_SGN_HARMONY_new.rds")
mouse_sgn= MOUSE_SGN_HARMONY_new
DefaultAssay(mouse_sgn)<- "RNA"
DefaultAssay(sgn_human.combined.sct) <- "RNA"
Idents(mouse_sgn) <- "orig.ident"
all <- merge(x=sgn_human.combined.sct,y=mouse_sgn)
all <- SCTransform(all, verbose = T, conserve.memory = T,
                   vars.to.regress = "nCount_RNA")

table(Idents(all))
Idents(all) <- factor(Idents(all), levels = c("Mouse_E14","Mouse_E15","Mouse_E16","Mouse_E17",
                                              "Mouse_E18","Mouse_P3","Mouse_P25","2","1","0"))

table(Idents(all))
all$devtime=Idents(all)


cochlea_genes <- read.csv("degs_mouse_sgn_devtime.csv",row.names = 1)
top30 <- cochlea_genes %>% group_by(cluster) %>% top_n(n = 30, wt = avg_diff) 
data_to_plot <- subset(all, downsample = 50)
data_to_plot <- ScaleData(data_to_plot,assay = "SCT")

data_to_plot2 <- subset(mouse_sgn, downsample = 50)
data_to_plot2 <- ScaleData(data_to_plot2,assay = "SCT")



library(dittoSeq)
dittoHeatmap(data_to_plot, top30$gene,cluster_rows=FALSE,
             annot.by = c("devtime"),complex = TRUE)

dittoHeatmap(data_to_plot, top30$gene,
             annot.by = c("devtime"),complex = TRUE)


dittoHeatmap(data_to_plot, top30$gene,scaled.to.max = TRUE,cluster_rows=FALSE,
             annot.by = c("devtime"),complex = TRUE)
dittoHeatmap(data_to_plot, top30$gene,scaled.to.max = TRUE,
             annot.by = c("devtime"),complex = TRUE)
###
DefaultAssay(sgn_human.combined.sct) <- "SCT"
sgn_human.combined.sct <- PrepSCTFindMarkers(sgn_human.combined.sct)
all.markers <- FindAllMarkers(sgn_human.combined.sct, 
                              assay = "SCT", slot = "data", test.use = "roc")

cochlea_genes <- all.markers[which(all.markers$avg_diff > 0), ]
write.csv(cochlea_genes,file = "degs_HUMAN_sgn_devtime.csv")#remove MT RIBOSOME GENE
cochlea_genes <- read.csv("degs_HUMAN_sgn_devtime.csv",row.names = 1)
top100 <- cochlea_genes %>% group_by(cluster) %>% top_n(n = 100, wt = avg_diff) 

data_to_plot <- subset(sgn_human.combined.sct, downsample = 50)
data_to_plot <- ScaleData(data_to_plot,assay = "SCT")
DoHeatmap(data_to_plot,features = top100$gene,assay = "SCT", slot = "scale.data")

data_to_plot2 <- subset(all, downsample = 50)
data_to_plot2 <- ScaleData(data_to_plot2,assay = "SCT")

library(dittoSeq)
dittoHeatmap(data_to_plot2, top30$gene,cluster_rows=FALSE,
             annot.by = c("devtime"),complex = TRUE)

dittoHeatmap(data_to_plot2, top30$gene,
             annot.by = c("devtime"),complex = TRUE)


dittoHeatmap(data_to_plot2, top30$gene,scaled.to.max = TRUE,cluster_rows=TRUE,
             annot.by = c("devtime"),complex = TRUE)
dev.off()
###############trajectory analysis######
sce=sgn_human.combined.sct
library(CytoTRACE)
library(paletteer)
pal <- paletteer_d("ggsci::nrc_npg")[c(5,3,1,4)] 
Idents(sce) <- "subtype"
DimPlot(sce, reduction = "umap",label = TRUE,repel = TRUE,
        pt.size = 2,cols = pal)

exp1 <- as.matrix(sce@assays$RNA@counts)
exp1 <- exp1[apply(exp1 > 0,1,sum) >= 5,]
results <- CytoTRACE(exp1,ncores = 1)
phenot <- sce$subtype
phenot <- as.character(phenot)
names(phenot) <- rownames(sce@meta.data)
emb <- sce@reductions[["umap"]]@cell.embeddings
plotCytoTRACE(results, phenotype = phenot, emb = emb,outputDir = './')
plotCytoGenes(results, numOfGenes = 50, outputDir = './')
###cytotrace2
devtools::install_github("digitalcytometry/cytotrace2", subdir = "cytotrace2_r") #installing
library(CytoTRACE2) #loading
exp1 <- as.matrix(sce@assays$RNA@counts)
exp1 <- exp1[apply(exp1 > 0,1,sum) >= 5,]
cytotrace2_result <- cytotrace2(sce,species = "human",is_seurat = TRUE,
                                slot_type = "counts",ncores = 1)
cytotrace2_result2 <- cytotrace2(exp1,species = "human",
                                slot_type = "counts",ncores = 1)
annotation=as.data.frame(sce@meta.data)
annotation <- annotation[,c(27,28)]
plots <- plotData(cytotrace2_result = cytotrace2_result,is_seurat = TRUE,
                  annotation = annotation,
                  expression_data = as.matrix(sce@assays$SCT@counts)
)
plots
plots2 <- plotData(cytotrace2_result = cytotrace2_result2,
                  annotation = annotation,
                  expression_data = as.matrix(sce@assays$SCT@counts)
)
plots2
plots$CytoTRACE2_UMAP
plots$CytoTRACE2_Boxplot_byPheno
plots$CytoTRACE2_Potency_UMAP
plots$CytoTRACE2_Relative_UMAP
plots$Phenotype_UMAP
###################MONOCLE 3 ANALYSIS ########################################
library(monocle3)
sgn_monocle3=sce

Idents(sgn_monocle3) <- "orig.ident"
table(Idents(sgn_monocle3))
data <- GetAssayData(sgn_monocle3,assay = 'SCT',slot='counts')
cell_metadata <- sgn_monocle3@meta.data
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
                 color_cells_by = 'subtype')+
  ggtitle('cds.umap')
p1  
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(sgn_monocle3, reduction = 'umap')
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP<- int.embed
p2 <- plot_cells(cds,reduction_method = 'UMAP',
                 color_cells_by = 'subtype')+
  ggtitle('sce.umap')
p2
p= p1|p2
p

cds <- cluster_cells(cds)
plot_cells(cds,color_cells_by = 'partition',reduction_method = 'UMAP')
cds <- learn_graph(cds)
p= plot_cells(cds,color_cells_by = 'subtype',label_groups_by_cluster = FALSE,
              label_leaves = FALSE,label_branch_points = FALSE)
p
plot_cells(cds,color_cells_by = 'subtype',label_groups_by_cluster = FALSE,
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
              graph_label_size=1.5,cell_size = 1.5)
p3
saveRDS(cds,file = "human_SGN_MONOCLE3_new.rds")
cds=human_SGN_MONOCLE3_new
p2+p3


plot_genes_in_pseudotime(cds[c("CALB1","PRPH","LYPD1","RUNX1","CACNA1B","ESRRG"),],color_cells_by = "pseudotime",
                         ncol=2,cell_size = 2) 

plot_genes_in_pseudotime(cds[c("PROX1","PBX3","ESRRG","PCDH9"),],
                         color_cells_by = "pseudotime",
                         ncol=2,cell_size = 2)

plot_genes_in_pseudotime(cds[c("PROX1","PBX3","ESRRG","PCDH9"),],
                         color_cells_by = "subtype",
                         ncol=2,cell_size = 2)

########################################################################################
#######slingshot analysis
##shut down monocle3!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

library(monocle)
library(slingshot)
library(viridis)
library(RColorBrewer)
library(pheatmap)
library(tibble)
library(stats)
library(colorRamps)
library(scales)
library(colorspace)
library(rcartocolor)
library(Seurat)
library(dplyr)
library(reshape2)
library(tidyverse)
library(gtools)
library(furrr)
library(readr)
library(tibble)
library(tidyr)
library(stringr)
sce=sgn_human.combined.sct
table(sce$subtype)
#############
Idents(sce) <- "subtype"
Idents(sce)
sce$celltype <- Idents(sce)
table(sce$celltype)
##############
table(sce$celltype,sce$subtype)
include_celltype <- c("Type_IA/B",
                      "Type_II_intermediate" ,
                      "Type_II" ,
                      "Type_IC"
                      ) 
table(sce$celltype)
sce_ALL <- sce


## ---------------------------------- #
real_colors <- c("Type_IA/B" = "#33A02C",
                 "Type_II_intermediate" = "#B2DF8A",
                 "Type_II" = "#55A1B1",
                 "Type_IC" = "#8DD3C7"
                
) 


umap_colors <- c("Type_IA/B" = "#DC050C",
                 "Type_II_intermediate" = "#FB8072",
                 "Type_II" = "#1965B0",
                 "Type_IC" = "#7BAFDE"
                
)
################################################################################################
## ---------------------------------- #
sce <- subset(sce, celltype %in% include_celltype )
root_cell <- "Type_II"
Idents(sce) <- sce$celltype %>% as.character()
as.character()

## colors
## ---------------------------------- #
palette = plasma(100)
palette_celltype = brewer.pal(n = 10, name = "Set3")

C <- palette_celltype[as.factor(sce$celltype)]
names(C) <- as.factor(sce$celltype)

color <- unique(C)
names(color) <- unique(names(C))


## run Slingshots 
## ---------------------------------- #
start.clus <- root_cell
reduction = 'umap'
sds     = slingshot(Embeddings(sce, reduction), clusterLabels = Idents(sce),
                    start.clus = start.clus )
sce@tools[['slingshot']] = SlingshotDataSet(sds)
pseudotime = slingPseudotime(sds)

sds_all = slingshot(Embeddings(sce_ALL, reduction), clusterLabels = Idents(sce_ALL),
                    start.clus = start.clus )

## Plot slingshot curves#
curves = colnames(pseudotime)


## Plot slingshot curves
## ---------------------------------- #
sds_all$reducedDim %>% rownames() -> cell_id
umap_colors[sce_ALL$celltype] -> cell_color
names(cell_color) <- colnames(sce_ALL)
cell_color[cell_id] -> U_C

## Plot slingshot curve  by cell type
## ---------------------------------- #
pseudotime_orig <- pseudotime
sds_orig <- sds

R_TIME <- pseudotime_orig %>% as.data.frame() %>% dplyr::mutate("cell_id" = rownames(.))
R_META <- sce_ALL[["celltype"]] %>% dplyr::mutate("cell_id" = rownames(.))
R_UMAP <- Embeddings(object = sce_ALL, reduction = "umap") %>% data.frame() %>% 
  dplyr::mutate("cell_id" = rownames(.)) %>% 
  left_join(., R_META, by="cell_id") %>% 
  left_join(., R_TIME, by="cell_id") 


make_umap <- function(df, color_column, color_list) {
  
  print(ggplot(df, aes(x=UMAP_1, y=UMAP_2, group = eval(parse(text=color_column)),
                       colour = eval(parse(text=color_column)))) +
          geom_point(size=1.5, alpha=0.7) + 
          theme_light() + theme_classic() +
          scale_colour_manual(values=color_list) + 
          theme(legend.position = "none", panel.grid.major = element_blank(),
                panel.grid.minor = element_blank()) +
          scale_alpha(guide = 'none'))# to remove extra legend 
  
}
make_umap(R_UMAP, "celltype", umap_colors)



ggplot(R_UMAP, aes(x=UMAP_1, y=UMAP_2, group = eval(parse(text="celltype")), 
                   colour = eval(parse(text="celltype")))) +
  geom_point(size=2, alpha=0.7) + 
  theme_light() + theme_classic() +
  scale_colour_manual(values=umap_colors) + 
  theme(legend.position = "none", panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_alpha(guide = 'none')

P1=ggplot(R_UMAP, aes(x=UMAP_1, y=UMAP_2, colour = eval(parse(text="Lineage1")) )) +
  geom_point(size=2, alpha=0.7) +
  theme_light() + theme_classic() +
  scale_colour_viridis_c(na.value="#D3D3D3", option = "C") +
  theme(legend.title = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_alpha(guide = 'none')
P2=DimPlot(sce,group.by = "celltype")
P1+P2

make_umap_c <- function(df, color_column) {
  
  print(ggplot(R_UMAP, aes(x=UMAP_1, y=UMAP_2, colour = eval(parse(text=color_column)) )) +
          geom_point(size=2, alpha=0.7) +
          theme_light() + theme_classic() +
          scale_colour_viridis_c(na.value="#D3D3D3", option = "C") +
          theme(legend.title = element_blank(),panel.grid.major = element_blank(),
                panel.grid.minor = element_blank()) +
          scale_alpha(guide = 'none'))# to remove extra legend 
  
}  

#ggsave(filename = "R_UMAP.svg", width= 5, height = 5)
make_umap_c(R_UMAP, "Lineage1")
#make_umap_c(R_UMAP, "Lineage2")

## Plot slingshot curve by cell type
## ---------------------------------- #
pseudotime_orig <- pseudotime
sds_orig <- sds
#pdf(file = 'slingshot_curves.separate.colored.by.time.all.umap.NEW.pdf', width = 14, height=7)

par(mfrow = c(1, 1))
for ( c_num in seq(1, length(curves))) {
  
  sds <- sds[,curves[c_num]]
  
  sds$reducedDim %>% rownames() -> cell_id
  sds_all$reducedDim %>% rownames() %>% as.data.frame() -> all_cell_id
  colnames(all_cell_id) <- "cell_id"
  
  pseudotime = slingPseudotime(sds)
  colors = palette[cut(pseudotime[,1], breaks = 100)]
  names(colors) <- cell_id
  colors %>% as.data.frame() %>% rownames_to_column() -> colors_df  
  colnames(colors_df) <- c("cell_id", "color")
  
  all_cell_id %>% 
    left_join(., colors_df, by="cell_id") %>% 
    mutate(color = ifelse(is.na(.[["color"]])==TRUE, "#D3D3D3", color)) -> new_colors 
  
  print(plot(sds_all$reducedDim, col = new_colors$color, pch = 16, cex = 1,
             main = curves[c_num] ) +
          lines(SlingshotDataSet(sds), linInd = c_num, lwd = 2, col = 'black'))
  
  sds <- sds_orig
  pseudotime <- pseudotime_orig
  
}

## Add pseudotimes to meta data of R object 
## ---------------------------------- #
## > R[[c("Lineage1","Lineage2")]] %>% as.data.frame() -> b
## > pseudotime %>% as.data.frame() -> a
## > all.equal(a,b)
## [1] TRUE
## ---------------------------------- #
for ( curve in curves ) {
  pseudotime_sub <- pseudotime[colnames(sce),curve]
  sce <- AddMetaData(object = sce,
                     metadata = pseudotime_sub,
                     col.name = curve
  )
}
## Condition density along pseudotime
## ---------------------------------- #
df <- data.frame(sce[["celltype"]], sce[["Lineage1"]]) 
colnames(df) <- c("celltype", "Lineage")
na.omit(df) -> df
ggplot(df, aes(x=Lineage, fill=celltype)) +
  geom_density(alpha=0.4) + theme_classic()+
  scale_fill_manual(values=C)

dev.off()



############################################################lineage1####################
L <- sce[["Lineage1"]] %>% deframe()
names(L) <- rownames(sce[["Lineage1"]])
L[!is.na(L)] %>% names() -> L2_cell
#R$Lineage2[!is.na(R$Lineage2)] %>% names() -> L2_cell
sce[, L2_cell] -> new_sce 
Idents(new_sce) <- new_sce$subclasses %>% as.character()
table(new_sce$celltype)
table(sce$celltype)

# Transfer into cds
# ---------------------------------- #
library(monocle)
sample_ann <-  sce@meta.data  
gene_ann <- data.frame(
  gene_short_name = rownames(sce@assays$RNA) , 
  row.names =  rownames(sce@assays$RNA) 
)
pd <- new("AnnotatedDataFrame",
          data=sample_ann)
fd <- new("AnnotatedDataFrame",
          data=gene_ann)
ct<-as(as.matrix(sce@assays$RNA@counts),'sparseMatrix')
sc_cds <- newCellDataSet(
  ct, 
  phenoData = pd, 
  featureData =fd,
  expressionFamily = negbinomial.size(),
  lowerDetectionLimit=1)
sc_cds
sc_cds <- detectGenes(sc_cds, min_expr = 1) 
cds=sc_cds
# Estimate size factor
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)#ignore the errors


# call Monocle2
# ---------------------------------- #
# install https://www.bioconductor.org/packages/3.14/bioc/src/contrib/Archive/monocle/
# refer to : http://cole-trapnell-lab.github.io/monocle-release/docs/#constructing-single-cell-trajectories
# refer to : https://bookdown.org/ytliu13207/SingleCellMultiOmicsDataAnalysis/monocle2.html#monocle2-process
# ---------------------------------- #
# select superset of feature genes as genes expressed in at least 5% of all the cells.
# ---------------------------------- #

cds <- detectGenes(cds, min_expr = 0.1)
fData(cds)$use_for_ordering <- fData(cds)$num_cells_expressed > 0.05 * ncol(cds)
cds_genes <-  row.names(subset(fData(cds),num_cells_expressed >= 10))


# get genes used for ordering cells 
# ---------------------------------- #
# while removing batch by using fullModelFormulaStr 
# https://www.biostars.org/p/316204/
# ---------------------------------- #
clustering_DEG_genes <- differentialGeneTest(cds[cds_genes,],
                                             fullModelFormulaStr = '~celltype',
                                             cores = 1)
clustering_DEG_genes1 <- filter(clustering_DEG_genes,status == "OK")

cds_ordering_df <- clustering_DEG_genes1 %>% filter(qval < 0.01 & use_for_ordering == TRUE) %>% arrange(qval)
cds_ordering_df[1:1000, ] %>% pull(gene_short_name)  %>%  as.character() -> cds_ordering_genes


# Clustering Genes by Pseudotemporal Expression Pattern by Monocle2
# ---------------------------------- #
# df = 1 : a linear fit 
# df = 2 : affords a little nonlinearity
# df = 3 : VGAM
# http://www2.uaem.mx/r-mirror/web/packages/VGAM/vignettes/categoricalVGAM.pdf
# ---------------------------------- #  

pData(cds)[["Lineage1"]] -> pData(cds)$Pseudotime
diff_test_res <- differentialGeneTest(cds[cds_ordering_genes,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime, df=3)", 
                                      cores = 1)
diff_test_res_filter=diff_test_res %>% filter(qval <0.01) %>% arrange(qval)
write.csv(diff_test_res_filter,file = "human_sgn_slingshot.DEgenes.lineage1.csv")
sig_gene_names <- row.names(diff_test_res_filter)

head(sig_gene_names)

a <- as.matrix(rownames(cds))
b <- as.matrix(sig_gene_names)
c <- merge(x=a,y=b,all=FALSE)
Time_genes <- top_n(c, n = 100) %>% pull(V1) %>% as.character()
plot_pseudotime_heatmap(cds[Time_genes,], num_clusters=4, 
                        show_rownames=T, return_heatmap=T)

dev.off()

############################################################################
Time_diff <- diff_test_res_filter[,c(5,3,4)] 


Time_genes <- Time_diff %>% pull(gene_short_name) %>% as.character()

head(Time_genes[1:100])


p=plot_pseudotime_heatmap(cds[Time_genes[1:100],], num_clusters=6, 
                          show_rownames=T, return_heatmap=T)

p

dev.off()

########################################################################


sgn_genes <- row.names(subset(fData(cds),gene_short_name %in% c("PRPH","GATA3","CALB1","CACNA1B",
                                                                 "RUNX1","LYPD1","EFNA5","PROX1",
                                                                 "SYT2","PCDH9","PBX3","ESRRG")))
sgn_genes_subset <- cds[c("PBX3","ESRRG")]
sgn_genes_subset <- cds[lung_genes,]
plot_genes_in_pseudotime(sgn_genes_subset,color_by = "Pseudotime",ncol = 2)#celltype

plot_genes_in_pseudotime(sgn_genes_subset,color_by = "celltype",ncol = 2)+
  scale_color_manual( 
    values=viridis(4))



####################################################################################

source("get_pseudotime_matrix.R")
# make plots
# ---------------------------------- #
hm <- get_pseudotime_matrix(cds[Time_genes[1:100],],  
                            cluster_rows = TRUE,
                            hclust_method = "ward.D",
                            num_clusters = 6,
                            hmcols = NULL,
                            add_annotation_row = NULL,
                            add_annotation_col = NULL,
                            show_rownames = FALSE,
                            use_gene_short_name = TRUE,
                            norm_method = "log",
                            scale_max=3,
                            scale_min=-3,
                            trend_formula = "~sm.ns(Pseudotime, df=3)", 
                            return_heatmap=TRUE,
                            cores=1)


bks = c(seq(min(hm), 0, length.out=ceiling(200/2) + 1),
        seq(max(hm)/200, max(hm),length.out=floor(200/2)))

my_color4 = plasma(length(bks))
my_color5 = colorRampPalette(rev(rcartocolor::carto_pal(7, "Sunset")))(length(bks))
my_color6 = colorRampPalette(rcartocolor::carto_pal(7, "ag_Sunset"))(length(bks))
my_color7 = colorRampPalette(rev(rcartocolor::carto_pal(7, "SunsetDark")))(length(bks))

my_color_set <- list(my_color4, my_color5, my_color6, my_color7)
my_color_name <- c("plasma", "Sunset", "ag_Sunset", "SunsetDark")



# cluster and re-order rows
# ---------------------------------- #
# ALL_HCS <- c( "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")
# ---------------------------------- #

ALL_HCS <- c("ward.D")

for ( sub_color in seq(1,length(my_color_set)))  {
  for ( ALL_HC in c(ALL_HCS) ) {	
    
    
    print(monocle::plot_pseudotime_heatmap(cds[Time_genes[1:100],],

                                           cluster_rows = TRUE,
                                           trend_formula = "~sm.ns(Pseudotime, df=3)",
                                           hclust_method = ALL_HC, 
                                           num_clusters = 5,
                                           hmcols = my_color_set[sub_color][[1]],
                                           scale_max = 3, 
                                           scale_min = -3,
                                           cores = 1,
                                           show_rownames = T,
                                           return_heatmap = FALSE))
    
  }
  
  
}
dev.off()  

colors = palette[cut(pData(cds)$Pseudotime, breaks = 100)]
phenoData(cds)[["color"]] <- colors
GENE_OF_INTEREST <- c("PRPH","GATA3","CALB1","CACNA1B",
                      "RUNX1","LYPD1","EFNA5","PROX1",
                      "SYT2","PCDH9","PBX3","ESRRG")

print(plot_genes_in_pseudotime(cds[GENE_OF_INTEREST,], color_by = "Lineage1" ) +         
        scale_color_viridis(option = "C") 
      #scale_color_viridis()
)



print(plot_genes_in_pseudotime(cds[GENE_OF_INTEREST,], color_by = "Lineage1" ) +
        scale_color_viridis(option = "B")
      #scale_color_viridis()
)


ALL_HCS <- c("ward.D")

for ( sub_color in seq(1,length(my_color_set)))  {
  for ( ALL_HC in c(ALL_HCS) ) {	
    
    
    print(monocle::plot_pseudotime_heatmap(cds[GENE_OF_INTEREST,],
                                           #add_annotation_col = "HCtype",
                                           cluster_rows = TRUE,
                                           trend_formula = "~sm.ns(Pseudotime, df=3)",
                                           hclust_method = ALL_HC, 
                                           num_clusters = 6,
                                           hmcols = my_color_set[sub_color][[1]],
                                           scale_max = 3, 
                                           scale_min = -3,
                                           cores = 1,
                                           show_rownames = T,
                                           return_heatmap = FALSE))
    
  }
  
  
}
dev.off() 
 ################################################################################################
#####################################################
################################################################################################
#transfer seurat object into annadata-readable format h5ad

library(sceasy)
library(reticulate)
use_condaenv('EnvironmentName')
sce=human_sgn_annotation
DefaultAssay(sce) <- "RNA"
sceasy::convertFormat(sce, from="seurat", to="anndata",
                      outFile='Hu_sgn_python.h5ad')







