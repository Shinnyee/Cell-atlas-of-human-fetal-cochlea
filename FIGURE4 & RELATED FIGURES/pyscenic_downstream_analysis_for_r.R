rm(list=ls())
library(Seurat)
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(scRNAseq)
library(patchwork)
library(ggplot2) 
library(stringr)
library(circlize)

# 加载SeuratData
human_sgn_for_slinghsot_analysis_v2 <- readRDS("F:/PROJECTS/PROJECT_HUMAN_FETAL_COCHLEAE/WORKPLACE/R/sc_pipeline/SGN analysis/human_sgn_for_slinghsot_analysis_v2.rds")
#可视化
seurat.data=human_sgn_for_slinghsot_analysis_v2
regulonAUC <- read.csv("human_fetal_sgn_auc_matrix.csv",row.names = 1)
head(colnames(seurat.data))
head(rownames(regulonAUC))
regulonAUC=t(regulonAUC)
sub_regulonAUC <- regulonAUC[,match(colnames(seurat.data),colnames(regulonAUC))]
dim(sub_regulonAUC)
seurat.data
#确认是否一致
identical(colnames(sub_regulonAUC), colnames(seurat.data))

cellClusters <- data.frame(row.names = colnames(seurat.data), 
                           seurat_clusters = as.character(seurat.data$seurat_clusters))
cellTypes <- data.frame(row.names = colnames(seurat.data), 
                        celltype = seurat.data$integrated_cluster)
head(cellTypes)
head(cellClusters)
sub_regulonAUC[1:4,1:4] 
#保存
save(sub_regulonAUC,cellTypes,cellClusters,seurat.data,
     file = 'for_rss_and_visual.Rdata')

#look for reported tfs:IKZF2,TBX2
regulonsToPlot = c('GATA3(+)')
regulonsToPlot %in% row.names(sub_regulonAUC)
seurat.data@meta.data = cbind(seurat.data@meta.data ,t(assay(sub_regulonAUC[regulonsToPlot,])))
# VISUALIZATION
Idents(seurat.data) <- "integrated_subclass"
p1 = DotPlot(seurat.data, features = unique(regulonsToPlot)) + RotatedAxis()
p2 = RidgePlot(seurat.data, features = regulonsToPlot , ncol = 2) 
p3 = VlnPlot(seurat.data, features = regulonsToPlot,pt.size = 0.1)
p4 = FeaturePlot(seurat.data,features = regulonsToPlot,label = TRUE)
wrap_plots(p3,p4)
p1
p2
p3
p4
library(dittoSeq)
dittoPlot(seurat.data, "GATA3(+)", group.by = "integrated_subclass",
          plots = c("vlnplot", "jitter"),jitter.size = 0.2,vlnplot.scaling = "width")
DefaultAssay(seurat.data) <-"SCT"
p5= FeaturePlot(seurat.data,features = c("GATA3"),
                label = TRUE,pt.size = 1)
p5

p6 = VlnPlot(seurat.data, features = c("GATA3"),pt.size = 0)
p6
p7= DotPlot(seurat.data,features = c("GATA3"),cols = c("lightgrey", "red"))
p7
p8=FeaturePlot(seurat.data,features = c("GATA3"),
                label = TRUE,cols = c("lightgrey", "red"))
wrap_plots(p3,p4,p7,p8)
sub_regulonAUC[1:4,1:2]
dim(sub_regulonAUC)
selectedResolution <- "celltype" # select resolution
cellsPerGroup <- split(rownames(cellTypes), 
                       cellTypes[,selectedResolution])
sub_regulonAUC <- sub_regulonAUC[onlyNonDuplicatedExtended(rownames(sub_regulonAUC)),] # 去除extened regulons
dim(sub_regulonAUC)
# Calculate average expression:
regulonActivity_byGroup <- sapply(cellsPerGroup,
                                  function(cells) 
                                    rowMeans(getAUC(sub_regulonAUC)[,cells]))
write.csv(regulonActivity_byGroup,file = "regulonActivity_byGroup_macaque.csv")

regulonActivity_byGroup <- read.csv("regulonActivity_byGroup_macaque.csv",row.names = 1)

pheatmap(regulonActivity_byGroup,cluster_cols = FALSE,cluster_rows = TRUE,scale="column")


rss <- calcRSS(AUC=getAUC(sub_regulonAUC), 
               cellAnnotation=cellTypes[colnames(sub_regulonAUC), 
                                               selectedResolution])
rss=na.omit(rss)
rssPlot <- plotRSS(rss)
plotly::ggplotly(rssPlot$plot)
plotRSS_oneSet(rss, setName = "IPhIBC_1")
plotRSS_oneSet(rss, setName = "IHC")

library(dplyr) 
rss=regulonActivity_byGroup
head(rss)
library(dplyr) 
df = do.call(rbind,
             lapply(1:ncol(rss), function(i){
               dat= data.frame(
                 path  = rownames(rss),
                 cluster =   colnames(rss)[i],
                 sd.1 = rss[,i],
                 sd.2 = apply(rss[,-i], 1, median)  
               )
             }))
df$fc = df$sd.1 - df$sd.2
top20 <- df %>% group_by(cluster) %>% top_n(20, fc)
rowcn = data.frame(path = top20$cluster) 
n = rss[top20$path,] 
#rownames(rowcn) = rownames(n)
library(viridis)
pheatmap(n,
         annotation_row = rowcn,
         show_rownames = T,cluster_cols = FALSE,fontsize = 10)


b <- read.csv("adj_human.csv",header = TRUE)

gata3_target <- b[which(b$TF == "GATA3"),]
a <- read.csv("reg_human.csv",header = TRUE)
write.csv(gata3_target,file = "gata3_targets_macaque.csv")

library(dittoSeq)
library(viridis)
genes= kdm5b_target[which(kdm5b_target$importance > 1),]
genes2= tbx2_target[which(tbx2_target$importance > 1),]
genes3= ikzf2_target[which(ikzf2_target$importance > 1),]
genes4= pbx3_target[which(pbx3_target$importance > 3),]
genes6= gata3_target[which(gata3_target$importance > 1),]
write.csv(genes6,file = "gata3_targets_human.csv")
genes5 <- c("APC","CTNND1","GPC5","GNAQ","MLLT3","CDK14","LGR5","LGR6",
            "KREMEN1","PPP3CA","TNIK","SOX6")
genes5 <- as.data.frame(genes5)
colnames(genes5) <- "target"
genes5 <- merge(x=genes4, y=genes5, by="target", all=FALSE)
head(genes4)
DefaultAssay(macaque_hc_recluster) <- "RNA"
dittoHeatmap(macaque_hc_recluster,genes4$target,annot.by = "integrated_subclass",
             scaled.to.max = TRUE,main = "Predicted_PBX3_targets"
)
Idents(macaque_hc_recluster) <- "integrated_subclass"
table(macaque_hc_recluster$integrated_subclass)
hc <- subset(macaque_hc_recluster, idents= "HC")
iphibc <- subset(macaque_hc_recluster, idents= "IPhIBC")
dcopc <- subset(macaque_hc_recluster,idents="DCOPC")
tbc <- subset(macaque_hc_recluster,idents="TBC")
pc <- subset(macaque_hc_recluster,idents="PC")
hec <- subset(macaque_hc_recluster,idents="HeC")
DefaultAssay(hc) <- "RNA"
dittoHeatmap(iphibc,genes4$target,annot.by = "integrated_subclass",
             scaled.to.max = TRUE,main = "Predicted_PBX3_targets"
             
)


###########################################################################
# Required packages:
library(SCopeLoomR)
library(AUCell)
library(SCENIC)

# For some of the plots:
#library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
packageVersion("SCENIC")

inputDir='C:/Users/Dell/Desktop/tf_analysis'
scenicLoomPath=file.path(inputDir,'sce_SCENIC2.loom')
motifEnrichmentFile <- file.path(inputDir,'reg.csv')
file.exists(scenicLoomPath)
file.exists(motifEnrichmentFile)
list.files()
library(SCopeLoomR)
loom <- open_loom(scenicLoomPath)
# Read information from loom file:
exprMat <- get_dgem(loom)
exprMat_log <- log2(exprMat+1) # Better if it is logged/normalized
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(loom)
close_loom(loom)
cellClusters = human_hc_recluster$integrated_cluster %>% as.data.frame()
colnames(cellClusters) = "integrated_cluster"
embeddings <- Embeddings(human_hc_recluster, reduction = "umap")

length(regulons)
head(names(regulons))
regulonAUC
length(regulonAucThresholds)
plot(embeddings)
motifEnrichment <- data.table::fread(motifEnrichmentFile, header=T, skip=1)[-3,]
colnames(motifEnrichment)[1:2] <- c("TF", "MotifID")
regulonAUC

head(motifEnrichment)
length(regulonAUC)
head(cellClusters)
# Split the cells by cluster:
table(human_hc_recluster$integrated_cluster)
Idents(human_hc_recluster) <- "integrated_cluster"
human_hc_recluster <- subset(human_hc_recluster,idents=c("OHC","IHC"))
cellClusters = human_hc_recluster$integrated_cluster %>% as.data.frame()
cellClusters$.
cellsPerCluster <- split(rownames(cellClusters), cellClusters[,"."])
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]

# Calculate average expression:
regulonActivity_byCellType <- sapply(cellsPerCluster,
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))

topRegulators <- reshape2::melt(regulonActivity_byCellType)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators$CellType <- factor(as.character(topRegulators$CellType))
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
dim(topRegulators)
viewTable(topRegulators, options = list(pageLength = 10))
selectedResolution <- "celltype"
rss <- calcRSS(AUC=getAUC(regulonAUC), 
               cellAnnotation=cellTypes[colnames(regulonAUC), 
                                               selectedResolution])


## Showing regulons and cell types with any RSS > 0.01 
rssPlot <- plotRSS(rss)
plotly::ggplotly(rssPlot$plot)

plotRSS_oneSet(rss, setName = "IHC") # cluster ID
embeddings=Embeddings(sce, reduction = "umap")
selectedEmbedding <- embeddings

tfsToPlot <- c("TBX2", "KDM5B","IKZF2","PBX3")
regulonsToPlot = c("TBX2(+)", "KDM5B(+)","IKZF2(+)","PBX3(+)")

par(mfrow=c(2,4))
# Plot expression:
library(viridis)
DefaultAssay(sce)<- "SCT"
AUCell::AUCell_plotTSNE(selectedEmbedding, exprMat_log[tfsToPlot,], plots=c("Expression"), cex = 1.8,
                        exprCols = c("grey","blue"))
# Plot regulon activity:
AUCell::AUCell_plotTSNE(selectedEmbedding, exprMat_log, regulonAUC[regulonsToPlot,], plots=c("AUC"), cex = 1.8)


head(as.data.frame(regulonAucThresholds))
regulonsToPlot <- "TBX2(+)"

par(mfrow=c(1,3))
AUCell::AUCell_plotTSNE(selectedEmbedding, exprMat_log, 
                        regulonAUC[regulonsToPlot,], thresholds = regulonAucThresholds[regulonsToPlot],
                        plots=c("AUC", "histogram", "binary"), cex = .5)


aucellApp <- AUCell_createViewerApp(auc=regulonAUC,
                                    thresholds=regulonAucThresholds,
                                    tSNE=selectedEmbedding, 
                                    exprMat=exprMat_log)
savedSelections <- shiny::runApp(aucellApp)
binarizeAUC <- function(auc, thresholds)
{
  thresholds <- thresholds[intersect(names(thresholds), rownames(auc))]
  regulonsCells <- setNames(lapply(names(thresholds), 
                                   function(x) {
                                     trh <- thresholds[x]
                                     names(which(getAUC(auc)[x,]>trh))
                                   }),names(thresholds))
  
  regulonActivity <- reshape2::melt(regulonsCells)
  binaryRegulonActivity <- t(table(regulonActivity[,1], regulonActivity[,2]))
  class(binaryRegulonActivity) <- "matrix"  
  
  return(binaryRegulonActivity)
}
binaryRegulonActivity <- binarizeAUC(regulonAUC, regulonAucThresholds)
dim(binaryRegulonActivity)
########################################################################################
#########################################################################################










