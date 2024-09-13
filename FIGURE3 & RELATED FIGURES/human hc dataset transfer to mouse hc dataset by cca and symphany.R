library(Seurat)
library(rtracklayer)
library(tibble)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(Hmisc)
options(stringsAsFactors = FALSE)
################################################################################################
################################################################################################
# In this script, we are going to anchor human hc snRNA-seq data to mouse hc snRNA-seq data and transfer the mouse cell type labels to each of the human cells
# Adapted from Seurat tutorial: https://satijalab.org/seurat/articles/integration_mapping.html

# ============ Load data ============ 
human_hc_GW_ALL <- readRDS("C:/Users/Dell/Desktop/WORKPLACE/R/HC TRAJOTORY/human_hc_GW_ALL.rds")# query
human_hc_GW_ALL <- readRDS("F:/Chai Lab/PROJECT_INTEGRATION ACROSS SPECIES/Cochlea_Evo_new/HC TRAJOTORY/human_hc_GW_ALL.rds") # query

mouse_hc_days_ALL <- readRDS("F:/Chai Lab/PROJECT_INTEGRATION ACROSS SPECIES/Cochlea_Evo_new/HC TRAJOTORY/mouse data/mouse_hc_days_ALL.rds") # reference
mouse_hc_days_ALL <- readRDS("C:/Users/Dell/Desktop/WORKPLACE/R/HC TRAJOTORY/mouse data/mouse_hc_days_ALL.rds")# reference
human_hc=human_hc_GW_ALL
mouse_hc=mouse_hc_days_ALL
human_hc$orig.ident2 <- "Human"
mouse_hc$orig.ident2 <- "Mouse"
# Names of HUMAN genes 
matrix_human=human_hc@assays$RNA@counts
meta_human=human_hc@meta.data

mouse_counts <- mouse_hc@assays$RNA@counts
mouse_counts <- data.frame(gene=rownames(mouse_counts), mouse_counts, check.names = F)
dim(mouse_counts)
# mouse genes transfer into one-to-one orthological human genes
genesV2 <- read.table("mouse_to_human_genes.txt", sep="\t", header=T)
mouse_counts$Gene <- genesV2[match(mouse_counts$gene, genesV2[,1]),2]
mouse_counts <- subset(mouse_counts, Gene!='NA')
mouse_counts <- dplyr::select(mouse_counts, Gene, everything())
mouse_counts <- mouse_counts[, !(colnames(mouse_counts) %in% 'gene')]
dim(mouse_counts)
write.csv(mouse_counts,file="rawcounts_mouse_to_human_transformed.csv")
meta_mouse=mouse_hc@meta.data

rownames(mouse_counts)<-mouse_counts[,1]
mouse_counts<-mouse_counts[,-1]
head(rownames(mouse_counts))


seurat_mouse_hc <- CreateSeuratObject(counts = mouse_counts)
seurat_mouse_hc
#preprocessing of human data
seurat_mouse_hc <- AddMetaData(seurat_mouse_hc, metadata = meta_mouse)
seurat_mouse_hc <- SCTransform(object = seurat_mouse_hc,vars.to.regress = c("nCount_RNA"),
                              conserve.memory = T)

seurat_mouse_hc[["harmony"]] <- mouse_hc[["harmony"]]
seurat_mouse_hc <- RunUMAP(seurat_mouse_hc,dims = 1:100,
                   reduction = "harmony",reduction.name = "umap_harmony"
)
DimPlot(seurat_mouse_hc,reduction = "umap_harmony",label = TRUE,repel = TRUE,
        group.by = "orig.ident",pt.size = 1.5)



# create a seurat object from human dataset that is compatible with mouse dataset
seurat_human_hc <- CreateSeuratObject(counts = matrix_human)
seurat_human_hc
#preprocessing of human data
seurat_human_hc <- AddMetaData(seurat_human_hc, metadata = meta_human)
seurat_human_hc <- SCTransform(object = seurat_human_hc,
                                   vars.to.regress = c("nCount_RNA"),conserve.memory = T)
seurat_human_hc <- RunPCA(seurat_human_hc)
seurat_mouse_hc <- RunPCA(seurat_mouse_hc)
seurat_human_hc[["harmony"]] <- human_hc[["harmony"]]
seurat_human_hc <- RunUMAP(seurat_human_hc,dims = 1:3,
                           reduction = "harmony",reduction.name = "umap_harmony"
)
DimPlot(seurat_human_hc,reduction = "umap_harmony",label = TRUE,repel = TRUE,
        group.by = "orig.ident",pt.size = 2)



# ============ projection of reference data onto query object ============ 
# We use all default parameters here for identifying anchors
transfer.anchors <- FindTransferAnchors(reference = seurat_mouse_hc, 
                                        query = seurat_human_hc, normalization.method = "SCT",
                                        k.filter = 100,
                                        reference.assay = "SCT",reference.reduction = "pca", 
                                        query.assay = "SCT",dims = 1:30
                                        )# reduction = "cca",

# TransferData returns a matrix with predicted IDs and prediction scores, which are written into meta.data
celltype.predictions <- TransferData(anchorset = transfer.anchors, 
                                     refdata = seurat_mouse_hc$orig.ident, 
                                     weight.reduction = "pcaproject",#k.weight = 20,
                                     dims = 1:30)#weight.reduction = seurat_human_hc[["pca"]]
seurat_human_hc <- AddMetaData(seurat_human_hc, metadata = celltype.predictions)


# Cells with low anchoring scores (<=0.25) are either low-quality cells or cells with ambigous identity.
# They are marked as low-confidence cells and removed later.
cutoff=0.25
seurat_human_hc$hiConfi=ifelse(seurat_human_hc$prediction.score.max > cutoff,'TRUE','FALSE')

# summary report of anchoring results
table(seurat_human_hc$predicted.id == seurat_human_hc$orig.ident)
seurat_human_hc@meta.data%>%mutate(agree.orig.ident=(predicted.id==orig.ident))%>%group_by(orig.ident)%>%summarise(n=sum(agree.orig.ident)/n())

saveRDS(seurat_human_hc,file="seurat_human_label_transfer_mouse_hc.rds")

cells <- rownames(seurat_human_hc@meta.data[seurat_human_hc@meta.data$prediction.score.max>0.25,])
DimPlot(seurat_human_hc, cells= cells, group.by="predicted.id", 
        label=T,reduction = "umap_harmony",pt.size = 2)

# ============ Intergration of two datasets for visualization ============ 

# note that we restrict the imputation to variable genes from scRNA-seq, but could impute the
# full transcriptome if we wanted to

genes.use <- VariableFeatures(seurat_mouse_hc)
refdata <- GetAssayData(seurat_mouse_hc, assay = "SCT", slot = "data")[genes.use, ]
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, 
                           weight.reduction = seurat_human_hc[["pca"]],dims = 1:30)
# this line adds the imputed data matrix to the seurat_macaque object
seurat_human_hc[["SCT"]] <- imputation
coembed <- merge(x = seurat_mouse_hc, y = seurat_human_hc)

# Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
# datasets
coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE,dims = 1:30)

coembed_harmony <- RunHarmony(coembed,group.by.vars = 'orig.ident2',reduction = "pca",
                              dims.use = 1:5,assay.use = "RNA")
coembed[["harmony"]] <- coembed_harmony[["harmony"]]
coembed <- RunUMAP(coembed,dims = 1:5,
                   reduction = "harmony",reduction.name = "umap_harmony")

DimPlot(object = coembed, reduction = "pca", pt.size = .1, group.by = "orig.ident")
DimPlot(object = coembed, reduction = "umap_harmony", pt.size = .1, group.by = "orig.ident")
# umap cluster colors
Idents(coembed) <- "orig.ident"
table(Idents(coembed))
levels(Idents(coembed))

my_cols <- c('iIHC'='#2B8DBF','IHC'='#F09150',
             'iOHC'='#C7EAE5',
             'OHC'='#A6D96A'
)
my_cols2 <- my_cols[order(as.integer(names(my_cols)))]
scales::show_col(my_cols2)
library(harmony)

DimPlot(coembed,
         label=TRUE , repel=TRUE,reduction = "umap_harmony",pt.size = 1)



my_cols3 <- c( "royalblue1","maroon4")
p1 <- DimPlot(coembed,reduction = "umap_harmony",label = TRUE,repel = TRUE, group.by = "orig.ident2",cols = my_cols3,pt.size = 1.5)
p2 <- DimPlot(coembed,reduction = "umap_harmony",label = TRUE,repel = TRUE,
              group.by = "celltype",cols = my_cols2,pt.size = 1.5)
p3 <- DimPlot(coembed,reduction = "umap_harmony",label = TRUE,repel = TRUE,
              group.by = "orig.ident",pt.size = 1.5)
p1+p2+p3



# Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
# datasets


coembed = AddMetaData(coembed,Embeddings(coembed[["umap_harmony"]]),
                      colnames(Embeddings(coembed[["umap_harmony"]])))

coembed$species=ifelse(!is.na(coembed$predicted.id), 'Human', 'Mouse')

# summary plots
d1 <- DimPlot(coembed, group.by = "species",pt.size = 1)
d2 <- DimPlot(coembed, cells= coembed@meta.data%>%rownames_to_column('V1')%>%filter(species=='Mouse')%>%pull(V1),group.by = "orig.ident", label = TRUE, repel = TRUE,pt.size = 1)+labs(title='Mouse subclasses')

d3 <- DimPlot(coembed, cells= coembed@meta.data%>%rownames_to_column('V1')%>%filter(species=='Human')%>%pull(V1),group.by = "orig.ident", label = TRUE, repel = TRUE,pt.size = 1)+labs(title='Human subclasses',pt.size = 1.5)
d4 <- DimPlot(coembed, cells= coembed@meta.data%>%rownames_to_column('V1')%>%filter(species=='Human')%>%pull(V1),group.by = "predicted.id", label = TRUE, repel = TRUE,pt.size = 1)+labs(title='Human predicted.id',pt.size = 1.5)

d5 = FeaturePlot(coembed, 'prediction.score.max',reduction = "umap_harmony",cols = c("grey","blue"),pt.size = 1)
d6 <- DimPlot(coembed, group.by = "hiConfi", label = TRUE, repel = TRUE, cols=c('TRUE'='green','FALSE'='red','NA'='transparent'),pt.size = 1)


grid.arrange(d1,d2,d3,d4,d5,d6,ncol=2)
DimPlot(coembed, cells= coembed@meta.data%>%rownames_to_column('V1')%>%pull(V1),group.by = "orig.ident", label = TRUE, repel = TRUE)+labs(title='subclasses')

saveRDS(coembed,file='Seurat_coembed_coch_human_mouse_label_transfer.rds')

table(seurat_human_hc$celltype)
table(seurat_human_hc$predicted.id)
table(human_hc$orig.ident,human_hc$celltype)
dev.off()

library(Seurat)
library(rtracklayer)
library(tibble)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(Hmisc)
grid.arrange(d1,d2,d3,d4,d5,d6,ncol=2)


table(coembed$orig.ident)
table(coembed$orig.ident,coembed$predicted.id)

dev.off()
################################################################################################
################################################################################################
results <- table(coembed@meta.data$orig.ident, coembed@meta.data$predicted.id)

library(Seurat)
library(dplyr)
library(cowplot)
library(zoo)
library(ggplot2)
library(loomR)
library(scater)
library(stringr)
library(corrplot)
library(matrixStats)

results.norm <- 100*(results/rowSums(results))
mus.ids <- colnames(results.norm)
results.norm <- results.norm[,mus.ids]
library(corrplot)
hist(coembed@meta.data$prediction.score.max, breaks=20)
corrplot(results.norm, order="original",tl.pos="lt", method="color", tl.col="black", is.corr=F)
dev.off()
################################################################################################
################################################################################################
############Symphony for harmonizing spatial-seq into P1 scRNA-seq data (from NC)
################################################################################################
################################################################################################
Idents(seurat_mouse_hc) <- "orig.ident"
table(Idents(seurat_mouse_hc))
P1=DimPlot(seurat_mouse_hc, reduction = "umap_harmony",label = TRUE,repel = TRUE)
Idents(seurat_human_hc) <- "orig.ident"
table(Idents(seurat_human_hc))
P2=DimPlot(seurat_human_hc, reduction = "umap_harmony",label = TRUE,repel = TRUE)
P1+P2
library(symphony)
library(Matrix)
library(ggplot2)
library(Seurat)
library(BGmix)
#prepare the reference dataset of their expression matrix and metadata info
ref_exp = seurat_mouse_hc@assays[["RNA"]]@counts
ref_metadata = seurat_mouse_hc@meta.data
ref_UMAP = seurat_mouse_hc@reductions[["umap_harmony"]]@cell.embeddings
colnames(ref_UMAP) = c("UMAP1","UMAP2")
ref_metadata = cbind(ref_UMAP,ref_metadata)
colnames(ref_metadata)
ref_metadata = ref_metadata[,c("UMAP1","UMAP2","orig.ident","celltype")]
ref_metadata$Sample = NA
ref_metadata$Sample = "reference"
################################################################################################
# Build reference
set.seed(0)
reference = symphony::buildReference(
  ref_exp,
  ref_metadata,
  #vars = c('batch'),         # variables to integrate over
  K = 100,                    # number of Harmony clusters
  verbose = TRUE,             # verbose output
  do_umap = TRUE,             # can set to FALSE if want to run umap separately later
  do_normalize = TRUE,        # set to TRUE if input counts are not normalized yet
  vargenes_method = 'vst',    # method for variable gene selection ('vst' or 'mvp')
  #vargenes_groups = 'batch', # metadata column specifying groups for variable gene selection 
  topn = 2000,                # number of variable genes to choose per group
  d = 20,                     # number of PCs
  save_uwot_path = './testing_uwot_model_2'  #
)
#replace embeddings from the 'testing_uwot_model_2'

testing_uwot_model_2 = uwot::load_uwot("testing_uwot_model_2", verbose = FALSE)
head(testing_uwot_model_2[["embedding"]])
head(seurat_mouse_hc@reductions[["umap_harmony"]]@cell.embeddings)
colnames(seurat_mouse_hc@reductions[["umap_harmony"]]@cell.embeddings) = c("UMAP1","UMAP2")
testing_uwot_model_2[["embedding"]] = seurat_mouse_hc@reductions[["umap_harmony"]]@cell.embeddings
#then delete the original 'testing_uwot_model_2'
uwot::save_uwot(model = testing_uwot_model_2,file = "testing_uwot_model_2",verbose = FALSE)

head(testing_uwot_model_2[["embedding"]])
head(seurat_mouse_hc@reductions[["umap_harmony"]]@cell.embeddings)

###########
reference$normalization_method = 'log(CP10k+1)' # optionally save normalization method in custom slot
# Save reference (modify with your desired output path)
#saveRDS(reference, './testing_reference2.rds')

head(reference[["umap"]][["embedding"]])
head(ref_UMAP)
reference[["umap"]][["embedding"]] = ref_UMAP

umap_labels = ref_metadata
#plotBasic(reference,ybar=umap_labels,ss =umap_labels)
ggplot(data=umap_labels,aes(UMAP1,UMAP2,colour=orig.ident))+
  geom_point(size=1,alpha=0.7)
#ggsave(file = "reference_umap.jpeg", width= 6, height = 4)
#ggsave(file = "reference_umap.pdf", width= 6, height = 2)



query_exp = seurat_human_hc@assays[["RNA"]]@counts
query_meta = seurat_human_hc@meta.data

table(query_meta$orig.ident)
query_exp[1:5,1:5]
query_meta = query_meta[,c("orig.ident","celltype")]


query_meta$Sample = NA
query_meta$Sample = "query"
colnames(ref_metadata)
colnames(query_meta) 
colnames(query_meta) = c("orig.ident","celltype","Sample")
# Map query
query = mapQuery(query_exp,             # query gene expression (genes x cells)
                 query_meta,            # query metadata (cells x attributes)
                 reference,             # Symphony reference object
                 do_normalize = TRUE,   # perform log(CP10k+1) normalization on query
                 do_umap = TRUE)        # project query cells into reference UMAP

query = knnPredict(query, reference, reference$meta_data$orig.ident, k = 5)

#################
table(query$meta_data$cell_type_pred_knn)
reference$meta_data$cell_type_pred_knn = NA
reference$meta_data$cell_type_pred_knn_prob = NA
umap=query$meta_data
umap=cbind(query$umap,umap)
query$meta_data = umap
colnames(reference[["meta_data"]])[1:2] = c("UMAP1","UMAP2")
meta_data_combined = rbind(query$meta_data,reference$meta_data)

#################

p1=ggplot(data=query$meta_data,aes(UMAP1,UMAP2,colour=orig.ident))+
  geom_point(size=2,alpha=0.7)
#ggsave(file = "query_umap.jpeg", width= 6, height = 4)
#ggsave(file = "query_umap.pdf", width= 6, height = 4)

p2=ggplot(data=query$meta_data,aes(UMAP1,UMAP2,colour=query$meta_data$cell_type_pred_knn))+
  geom_point(size=2,alpha=0.7)
#ggsave(file = "query_umap_pred.jpeg", width= 6, height = 4)
#ggsave(file = "query_umap_pred.pdf", width= 6, height = 4)
p1+p2

########################################
#correlation plot
table(query$meta_data$cell_type_pred_knn)
table(query$meta_data$orig.ident)
plot_heatmap = query$meta_data
write.csv(plot_heatmap,file = "plot_heatmap_integrated_symphony.csv")


############################################# Confusion heatmap  #############################################################################
##############################################################################################################################################
library(Seurat)
library(dplyr)
library(Matrix)
library(matrixStats)
library(feather)
sample.combined = query
sample.combined$label_for_heatmap <- paste(
  sample.combined$meta_data$orig.ident, sep = "_")

reorder_matrix <- function(matrix1, by.rows = TRUE) {
  if (by.rows == TRUE) {
    conf.order <- order(apply(matrix1, 1, which.max))
    matrix1.reordered <- matrix1[conf.order, ]
  } else {
    conf.order <- order(apply(matrix1, 2, which.max))
    matrix1.reordered <- matrix1[, conf.order]
  }
}

compare_cl <- function(cl, ref.cl,
                       plot.title = NA, plot.silent = TRUE,
                       heat.colors = colorRampPalette(c("white", "grey70", "black"))(100),
                       row.cl.num = min(length(unique(cl)),
                                        length(unique(ref.cl)))) {
  library(grid)
  library(pheatmap)
  
  conf1 <- table(cl, ref.cl)
  conf1 <- sweep(conf1, 1, rowSums(conf1), "/")
  conf2 <- reorder_matrix(conf1)
  
  # Cluster co-occurence
  cl.prop.cocl <- apply(conf1, 2, function(x) {
    grid1 <- expand.grid(x, x)
    min.prop <- apply(grid1, 1, min)
  })
  
  cl.prop.cocl.total <- apply(cl.prop.cocl, 1, sum)
  cl.prop.cocl.m <- matrix(cl.prop.cocl.total, nrow(conf1), nrow(conf1),
                           dimnames = list(rownames(conf1), rownames(conf1)))
  
  ph1 <- pheatmap(conf2, cutree_rows = row.cl.num, clustering_method = "ward.D2",
                  # annotation_row = ref.cl.anno[, -grep("cluster_label", colnames(ref.cl.anno))],
                  color = heat.colors, fontsize = 6,
                  main = plot.title, silent = plot.silent)
  return(list(conf = conf2, cocl = cl.prop.cocl.m, ph = ph1))
}

# Heatmap palette
#heat.colors <- colorRampPalette(c("grey99", "orange", "red"))(100)
#heat.colors <- colorRampPalette(c("white", "#95f0f0", "#00fcfc"))(100)
heat.colors <- colorRampPalette(c("white", "grey70", "black"))(100)
#"white", "grey70", "black"
# Compare clustering
ref.cl <- sample.combined$label_for_heatmap
cca.cl <- sample.combined$meta_data$cell_type_pred_knn

cl.conf <- compare_cl(ref.cl, cca.cl)
cocl <- cl.conf$conf
a <- as.data.frame.array(cl.conf[["conf"]])
write.csv(a, file = "cell_proportion_integrated_clusters_symphony.csv")
a <- read.csv(file = "cell_proportion_integrated_clusters_symphony.csv",row.names = 1)
library(fpc)
clus.method <- "single"

clus.num <- pamk(cocl, 1:(min(nrow(cocl), ncol(cocl)) - 1))$nc

ph1 <- pheatmap(a, clustering_method = clus.method,
                cutree_cols = clus.num, cutree_rows = clus.num,
                color = heat.colors,
                fontsize = 6)

pheatmap(a, cluster_rows = FALSE, cluster_cols = FALSE,
         color = heat.colors,
         fontsize = 8) 

results <- table(query$meta_data$orig.ident,query$meta_data$cell_type_pred_knn)

library(Seurat)
library(dplyr)
library(cowplot)
library(zoo)
library(ggplot2)
library(loomR)
library(scater)
library(stringr)
library(corrplot)
library(matrixStats)

results.norm <- 100*(results/rowSums(results))
mus.ids <- colnames(results.norm)
results.norm <- results.norm[,mus.ids]
library(corrplot)
hist(coembed@meta.data$prediction.score.max, breaks=20)
corrplot(results.norm, order="original",tl.pos="lt", method="color", tl.col="black", is.corr=F)
dev.off()

