rm(list = ls())
library(Seurat)
library(rtracklayer)
library(tibble)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(Hmisc)
human_sgn_for_slinghsot_analysis <- readRDS("F:/PROJECTS/PROJECT_HUMAN_FETAL_COCHLEAE/WORKPLACE/R/sc_pipeline/SGN analysis/human_sgn_for_slinghsot_analysis.rds")

sce=human_sgn_for_slinghsot_analysis

#Find DE genes

Idents(sce) <- "cell_type_final"
table(Idents(sce))
my_levels <- c('NEUROD1+_early_sgn','IGFBP2+_intermediate_sgn',
               'LYPD1+_immature_sgn','PRPH+_immature_sgn','ALDH1A3+_immature_sgn',
               'RUNX1+_mature_sgn',"CLU+_mature_sgn")
Idents(sce) <- factor(Idents(sce),levels = my_levels)
#Idents(hc_subtype) <- factor(Idents(hc_subtype), 
#                             levels = c("E14HC","E16HC","P1HC","P7HC","P14HC","P28HC"))#变换idents顺序
#Idents(hc_subtype)
sce$cell_type_final2 <- Idents(sce)
DefaultAssay(sce) <- "RNA"
sce <- SCTransform(sce, method = "glmGamPoi")
# 线粒体基因通常以 "MT-" 开头
mito_genes <- grep("^MT-", rownames(sce), value = TRUE)
# 核糖体基因通常以 "RPS" 或 "RPL" 开头
ribo_genes <- grep("^RPS|^RPL", rownames(sce), value = TRUE)
# 血细胞相关基因可以根据你的数据集和注释进行调整
# 这里假设你有特定的基因列表
blood_genes <- c("HBA1", "HBA2", "HBB", "CD3D",
                 "CD3E", "CD4", "CD8A", "CD14", "CD19", "CD20", "CD34")

genes_to_remove <- unique(c(mito_genes, ribo_genes, blood_genes))

gene_list <- read.csv("mouse_to_human_genes.csv")
gene_list <- gene_list$HGNC.symbol[gene_list$HGNC.symbol %in% rownames(sce)]
sce_filtered <- subset(sce, features = gene_list)
sce_filtered <- subset(sce_filtered, features = setdiff(rownames(sce_filtered), genes_to_remove))


sce_cells_markers <- FindAllMarkers(sce_filtered, assay = "RNA", 
                                    slot = "counts", test.use = "MAST")

my_cols <- c("#9b5792", "#575add", "#74c5e9",
             "#62edd0","#1db342", "#b7de27", 
             "#e8ba69","#a056d3","#ba1c1c"
             
)
my_cols2 <- my_cols[order(as.integer(names(my_cols)))]
scales::show_col(my_cols)
DimPlot(sce_filtered,
        cols = my_cols, label=FALSE , repel=TRUE,reduction = "umap",pt.size = 2
        ,group.by = "cell_type_final")#
sce_cells_markers_sub <- sce_cells_markers[which(sce_cells_markers$avg_log2FC > 0), ]

write.csv(sce_cells_markers_sub,
          file="human_sgn_DEGs_devtime_mast_test.csv")
sce_cells_markers_sub <- read.csv(file="human_sgn_DEGs_devtime_mast_test.csv",row.names = 1)
###################################################################################
BiocManager::install("ComplexHeatmap",force = T)

library(ClusterGVis)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(clusterProfiler)
#devtools::install_github('cole-trapnell-lab/monocle3')
#devtools::install_github("junjunlab/ClusterGVis")
group <- data.frame(gene=sce_cells_markers_sub$gene,
                    group= sce_cells_markers_sub$cluster)
Gene_ID <- bitr(sce_cells_markers_sub$gene,fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = "org.Hs.eg.db")
write.csv(Gene_ID,file = 'gene_id.csv')
data  <- merge(Gene_ID,group,by.x='SYMBOL',by.y='gene')
enrich2 <- compareCluster(
  ENTREZID~group, 
  data=data, 
  fun="enrichGO", 
  OrgDb="org.Hs.eg.db",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.5
)#qvalueCutoff = 0.05
enrich2_discription <- enrich2@compareClusterResult[!duplicated(enrich2@compareClusterResult[,3]),] %>%
  dplyr::group_by(Cluster) %>%
  dplyr::top_n(n=10,wt=Count) %>% dplyr::top_n(n=10,wt=ID)
enrich2_discription <- enrich2_discription[!duplicated(enrich2_discription[,3]),]
enrich2_discription <- enrich2_discription[,c(3,2,4,7,5)]


enrich2_discription1 <- enrich2@compareClusterResult[!duplicated(enrich2@compareClusterResult[,3]),] %>%
  dplyr::group_by(Cluster)

write.csv(enrich2_discription1,file = "goenrich_human_sgn_cell_type.csv")
# add gene name
markGenes = unique(sce_cells_markers_sub$gene)[sample(1:length(unique(sce_cells_markers_sub$gene)),50,
                                                  replace = F)]
# no average cells
cochlea.markers1 <- sce_cells_markers_sub %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 30, wt = avg_log2FC)
# retain duplicate diff gene in multiple clusters
sce2=subset(sce_filtered2, downsample=30)
DefaultAssay(sce2) <- "RNA"
sce2 <- SCTransform(sce2)
sce2 <- ScaleData(sce2,assay = "SCT")
Idents(sce2)
table(sce2$cell_type_final2)

st.data <- prepareDataFromscRNA(object = sce2,
                                diffData = cochlea.markers1,
                                showAverage = FALSE)
str(st.data)
st.data1 <- prepareDataFromscRNA(object = sce_filtered2,
                                 diffData = cochlea.markers1,
                                 showAverage = FALSE)
str(st.data1)
enrich <- enrichCluster(object = st.data,
                        OrgDb = org.Hs.eg.db,
                        type = "BP",
                        organism = "Hsa",
                        pvalueCutoff = 0.5,
                        topn = 10,
                        seed = 5201314,
                        readable = TRUE)
# add GO annotation
pdf('HEATMAP_human_sgn_DEGS_GO.pdf',height = 30,width = 30,onefile = F)
library(viridis)
visCluster(object = st.data,
           plot.type = "both",
           column_title_rot = 45,
           markGenes = unique(cochlea.markers1$gene),
           markGenes.side = "left",
           annoTerm.data = enrich,
           ht.col.list = list(col_range = seq(-1,1,length.out = 100),
                              col_color = inferno(100)),
           go.size = 15,
           border = FALSE,
           line.side = "left",
           cluster.order = c(1:7),go.col = "black",
           add.bar = T)


dev.off()
pdf('HEATMAP_human_sgn_DEGS_GO2.pdf',height = 30,width = 35,onefile = F)
visCluster(object = st.data,
           plot.type = "both",
           column_title_rot = 45,
           markGenes = unique(cochlea.markers1$gene),
           markGenes.side = "left",
           annoTerm.data = enrich,
           go.size = 15,
           border = FALSE,
           line.side = "left",
           cluster.order = c(1:7),go.col = "black",
           add.bar = T)
dev.off()
