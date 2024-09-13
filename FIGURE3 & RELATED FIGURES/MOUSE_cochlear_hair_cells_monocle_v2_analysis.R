
library(Seurat)
library(ggplot2)
library(sctransform)
library(tidyverse)
options(stringsAsFactors = FALSE)
library(dittoSeq)
mouse_e14_hc <- readRDS("F:/PROJECTS/PROJECT_HUMAN FETAL COCHLEAE/WORKPLACE/R/HC TRAJOTORY/mouse data/mouse_e14_hc.rds")
mouse_e16_hc <- readRDS("F:/PROJECTS/PROJECT_HUMAN FETAL COCHLEAE/WORKPLACE/R/HC TRAJOTORY/mouse data/mouse_e16_hc.rds")
hc_P1 <- readRDS("F:/PROJECTS/PROJECT_HUMAN FETAL COCHLEAE/WORKPLACE/R/HC TRAJOTORY/mouse data/hc_P1.rds")
mouse_p7_hc <- readRDS("F:/PROJECTS/PROJECT_HUMAN FETAL COCHLEAE/WORKPLACE/R/HC TRAJOTORY/mouse data/mouse_p7_hc.rds")
hc_P14 <- readRDS("F:/PROJECTS/PROJECT_HUMAN FETAL COCHLEAE/WORKPLACE/R/HC TRAJOTORY/mouse data/hc_P14.rds")
HC_P28 <- readRDS("F:/PROJECTS/PROJECT_HUMAN FETAL COCHLEAE/WORKPLACE/R/HC TRAJOTORY/mouse data/HC_P28.rds")
################################################################################################
################################################################################################
################################################################################################

e14 = mouse_e14_hc
e16 = mouse_e16_hc
p1 = hc_P1
p7 = mouse_p7_hc
p14 = hc_P14
p28 = HC_P28
table(e14$celltype)
table(e16$celltype)
table(p1$celltypes)
table(p7$celltype)
table(p14$celltypes)
table(p28$celltypes)
#rename cluster_label
Idents(p1) <- "celltypes"
table(Idents(p1))
new.cluster.ids <- c("iOHC", "iOHC","iIHC", "iOHC"
)
names(new.cluster.ids) <- levels(p1)
p1<- RenameIdents(p1, new.cluster.ids)
p1$celltype <- Idents(p1)
Idents(p1) <- "celltype"
Idents(p7) <- "celltype"
table(Idents(p7))

Idents(p14) <- "celltypes"
table(Idents(p14))
new.cluster.ids <- c("OHC", "OHC","IHC", "OHC"
)
names(new.cluster.ids) <- levels(p14)
p14<- RenameIdents(p14, new.cluster.ids)
p14$celltype <- Idents(p14)
Idents(p14) <- "celltype"
Idents(p28) <- "celltypes"
table(Idents(p28))
new.cluster.ids <- c("OHC","IHC", "OHC"
)
names(new.cluster.ids) <- levels(p28)
p28<- RenameIdents(p28, new.cluster.ids)
p28$celltype <- Idents(p28)
Idents(p28) <- "celltype"
table(e14$orig.ident)
table(e16$orig.ident)
table(p1$orig.ident)
table(p7$orig.ident)
table(p14$orig.ident)
table(p28$orig.ident)

p28$orig.ident <- "mouse_P28"
#####we used harmony to batch 6 timepoints mouse hc dataset
library(Seurat)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(sctransform)
library(tidyverse)
library(harmony)
options(stringsAsFactors = FALSE)
library(tidyverse)

matrix_mouse_e14=e14@assays$RNA@counts
meta_mouse_e14=e14@meta.data
matrix_mouse_e16=e16@assays$RNA@counts
meta_mouse_e16=e16@meta.data
matrix_mouse_p1=p1@assays$RNA@counts
meta_mouse_p1=p1@meta.data
matrix_mouse_p7=p7@assays$RNA@counts
meta_mouse_p7=p7@meta.data
matrix_mouse_p14=p14@assays$RNA@counts
meta_mouse_p14=p14@meta.data
matrix_mouse_p28=p28@assays$RNA@counts
meta_mouse_p28=p28@meta.data

mouse_hc_e14 <- CreateSeuratObject(counts = matrix_mouse_e14)
mouse_hc_e14
mouse_hc_e14 <- AddMetaData(mouse_hc_e14, metadata = meta_mouse_e14)

mouse_hc_e16 <- CreateSeuratObject(counts = matrix_mouse_e16)
mouse_hc_e16
mouse_hc_e16<- AddMetaData(mouse_hc_e16, metadata = meta_mouse_e16)

mouse_hc_p1 <- CreateSeuratObject(counts = matrix_mouse_p1)
mouse_hc_p1
mouse_hc_p1 <- AddMetaData(mouse_hc_p1, metadata = meta_mouse_p1)
mouse_hc_p7 <- CreateSeuratObject(counts = matrix_mouse_p7)
mouse_hc_p7
mouse_hc_p7 <- AddMetaData(mouse_hc_p7, metadata = meta_mouse_p7)
mouse_hc_p14 <- CreateSeuratObject(counts = matrix_mouse_p14)
mouse_hc_p14
mouse_hc_p14 <- AddMetaData(mouse_hc_p14, metadata = meta_mouse_p14)
mouse_hc_p28 <- CreateSeuratObject(counts = matrix_mouse_p28)
mouse_hc_p28
mouse_hc_p28 <- AddMetaData(mouse_hc_p28, metadata = meta_mouse_p28)




sce_all <- merge(x = mouse_hc_e14, 
                 y = c(mouse_e16_hc,mouse_hc_p1,mouse_hc_p7,mouse_hc_p14,mouse_hc_p28 ))

sce_all
sce_all <- SCTransform(sce_all,vars.to.regress = "nCount_RNA")

DefaultAssay(sce_all) <- "SCT"
sce_all <- sce_all %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData()

sce_all <- RunPCA(sce_all, features = VariableFeatures(sce_all), npcs = 100)

DimPlot(object = sce_all, reduction = "pca", pt.size = .1, group.by = "celltype")

sce_all_harmony <- RunHarmony(sce_all,group.by.vars = 'orig.ident',reduction = "pca",
                              dims.use = 1:20,assay.use = "SCT"
                             )

sce_all[["harmony"]] <- sce_all_harmony[["harmony"]]
sce_all <- RunUMAP(sce_all,dims = 1:5,
                   reduction = "harmony",reduction.name = "umap_harmony"
                   )

p3 <- DimPlot(sce_all, reduction = "umap_harmony", group.by = "orig.ident",pt.size = 1) + 
  ggtitle("UMAP Harmony")
p3
dev.off()
p4 <- DimPlot(sce_all, reduction = "umap_harmony",pt.size = 1,
              group.by = "celltype",repel = TRUE) + 
  ggtitle("UMAP Harmony")
p4
p3+p4
p5 <- DimPlot(sce_all, reduction = "umap_harmony",pt.size = 1,
              split.by = "celltype",repel = TRUE,group.by = "orig.ident") + 
  ggtitle("UMAP Harmony")
p5
p3+p4
saveRDS(sce_all,file = "mouse_hc_days_ALL.rds")
mouse_hc_days_ALL <- readRDS("F:/PROJECTS/PROJECT_HUMAN FETAL COCHLEAE/WORKPLACE/R/HC TRAJOTORY/mouse data/mouse_hc_days_ALL.rds")
sce_all=mouse_hc_days_ALL
sce_all <- SCTransform(sce_all, verbose = T, conserve.memory = T,vars.to.regress = 
                         c("nCount_RNA"))
DefaultAssay(sce_all) <- "SCT"
markers.to.plot<-c(
  "Sox2","Atoh1","Insm1",#Early HC
  "Fgfr3","Slc26a5","Ocm","Ikzf2",   
  "Otof","Slc17a8","Slc7a14","Atp2a3" ,
  "Myo7a","Pvalb","Myo6","Pcp4","Tmc1","Atp2b2",
  "Kdm5b","Zfp462","Tox","Npas3","Tbx2"
  
)
markers.to.plot<-c(
  "Calb2","Otof","Ikzf2",
  "Atoh1","Tbx2","Slc26a5","Fgfr3",   
  "Slc17a8","Tmc1"
  
)

Idents(sce_all) <- "celltype"
VlnPlot(sce_all,features = markers.to.plot,group.by = "celltype")

DotPlot(sce_all, features = markers.to.plot,
        dot.scale = 8,cols  =c("white", "#ad9300")
        ,group.by = "celltype") + RotatedAxis()

FeaturePlot(sce_all,features = markers.to.plot,label = TRUE)
Idents(sce_all) <- "orig.ident"
FeaturePlot(sce_all,features =c("Sorbs2","Ikzf2","Tbx2","Kdm5b"),label = TRUE)
saveRDS(sce_all, file = "mouse_hc_days_ALL.rds")
####umap time point##############################################
#UMAP species
sce=mouse_hc_days_ALL
sce=sce_all
data_to_plot <- data.frame(sce@reductions$umap@cell.embeddings)
sce$new_id <- colnames(sce)
data_to_plot$species <- sce$orig.ident[match(rownames(data_to_plot), sce$new_id)]
data_to_plot$plot_order <- 0
for(i in 1:nrow(data_to_plot)){
  data_to_plot$plot_order[i] <- sample(1:100000, 1)
  print(i)
}
data_to_plot <- arrange(data_to_plot, plot_order)

p1=ggplot(data_to_plot,
          aes(x = UMAP_1, y = UMAP_2, color = species)) +
  geom_point(size = 1.5)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_color_manual(values = c("#aa9ac6", "#3583ec", "#23d8b7","#d6ec7d","#f2c8a7","#b35555"))
p1
DimPlot(sce, reduction = "umap_harmony",label = TRUE,pt.size = 1.5,group.by = "celltype")
Idents(sce) <- "celltype"
DimPlot(sce, reduction = "umap_harmony",label = TRUE,pt.size = 1.5,split.by = "orig.ident",
        group.by = "celltype")
# umap cluster colors
table(sce$celltype)
levels(Idents(sce))

my_cols <- c("iIHC"="#e39a8c","IHC"="#a056d3",
             
             
             "iOHC"="#8cb6e3","OHC"="#cc4512"
             
)



my_cols2 <- my_cols[order(as.integer(names(my_cols)))]
scales::show_col(my_cols2)
p2=DimPlot(sce,
           cols = my_cols2, label=TRUE , repel=TRUE,reduction = "umap_harmony",pt.size = 1.5
           ,group.by = "celltype")
p2
p1+p2
#"royalblue1","maroon4","sienna2", "#ae5a8c"
################################################################################

###################MONOCLE 3 ANALYSIS ########################################
library(monocle3)

#sce=sce_all
Idents(sce) <- "celltype"
table(Idents(sce))

data <- GetAssayData(sce,assay = 'RNA',slot='counts')
cell_metadata <- sce@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 100)
plot_pc_variance_explained(cds)
cds
cds <- reduce_dimension(cds,preprocess_method = 'PCA',reduction_method = 'UMAP')
plot_cells(cds)
colnames(colData(cds))
p1 <- plot_cells(cds,reduction_method = 'UMAP',
                 color_cells_by = 'celltype')+
  ggtitle('cds.tsne')
p1  
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(sce, reduction = 'umap_harmony')
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP<- int.embed
p2 <- plot_cells(cds,reduction_method = 'UMAP',
                 color_cells_by = 'celltype')+
  ggtitle('sce.tsne')
p2
p= p1|p2
p
genes_of_interest <- c("Slc26a5","Slc17a8","Tmc1",
                       "Sox2","Insm1","Atoh1")
plot_cells(cds, genes = genes_of_interest)
cds <- cluster_cells(cds)
plot_cells(cds,color_cells_by = 'partition',reduction_method = 'UMAP')
cds <- learn_graph(cds)
p= plot_cells(cds,color_cells_by = 'celltype',label_groups_by_cluster = FALSE,
              label_leaves = FALSE,label_branch_points = FALSE)
p
plot_cells(cds,color_cells_by = 'celltype',label_groups_by_cluster = FALSE,
           label_leaves = TRUE,label_branch_points = TRUE,graph_label_size = 2,
           label_principal_points = TRUE)
cds <- order_cells(cds)
#cds_sub <- choose_graph_segments(cds)
p=plot_cells(cds,color_cells_by = 'pseudotime',label_cell_groups = FALSE,
             label_leaves = FALSE,label_branch_points = FALSE,cell_size = 1.2)
p



plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1)
#########DE analysis in monocle3
genes <- c("Slc26A5", "Ocm", "Ikzf2",
           "Slc17a8","Slc7a14","Tbx2","Kdm5b",
           "Tmc1",
           "Sox2","Insm1","Atoh1")
cds_subset <- cds[rowData(cds)$gene_short_name %in% genes,]

gene_fits <- fit_models(cds_subset, model_formula_str = "~pseudotime")
fit_coefs <- coefficient_table(gene_fits)
time_terms <- fit_coefs %>% filter(term == "pseudotime")
time_terms %>% filter (q_value < 0.05) %>%
  select(gene_short_name, term, q_value, estimate)
plot_genes_violin(cds_subset, group_cells_by="celltype", ncol=2) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
####################################MONOCLE 2 ANALYSIS########################
mouse_hc_days_ALL <- readRDS("F:/PROJECTS/PROJECT_HUMAN FETAL COCHLEAE/WORKPLACE/R/HC TRAJOTORY/mouse data/mouse_hc_days_ALL.rds")
sce=mouse_hc_days_ALL
sce=sce_all
library(monocle)
table(Idents(sce))
Idents(sce) <- 'orig.ident'
sample_ann <-  sce@meta.data  
head(sample_ann)
gene_ann <- data.frame(
  gene_short_name = rownames(sce@assays$RNA) , 
  row.names =  rownames(sce@assays$RNA) 
)
head(gene_ann)
pd <- new("AnnotatedDataFrame",
          data=sample_ann)
fd <- new("AnnotatedDataFrame",
          data=gene_ann)
ct=as.data.frame(sce@assays$RNA@counts)
ct[1:4,1:4]
sc_cds <- newCellDataSet(
  as.matrix(ct), 
  phenoData = pd,
  featureData =fd,
  expressionFamily = negbinomial.size(),
  lowerDetectionLimit=1)
sc_cds
sc_cds <- detectGenes(sc_cds, min_expr = 0.1) 
sc_cds <- sc_cds[fData(sc_cds)$num_cells_expressed >= 5, ]
sc_cds
cds <- sc_cds
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds) 
# genes selective for clustering
disp_table <- dispersionTable(cds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.5)
cds <- setOrderingFilter(cds, unsup_clustering_genes$gene_id)
plot_ordering_genes(cds) 
plot_pc_variance_explained(cds, return_all = F) # norm_method='log'

cds <- reduceDimension(cds, max_components = 2, num_dim = 10,
                       reduction_method = 'tSNE', verbose = T)
cds <- clusterCells(cds, num_clusters =10) 
plot_cell_clusters(cds, 1, 2 )
p1<-plot_cell_clusters(cds, 1, 2 , color_by = "celltype")
p2<-plot_cell_clusters(cds, 1, 2 , color_by = "orig.ident")
p1+p2
table(pData(cds)$Cluster) 
colnames(pData(cds)) 

table(pData(cds)$Cluster,pData(cds)$new)
plot_cell_clusters(cds, 1, 2 )

library(monocle)
colnames(pData(cds))
table(pData(cds)$Cluster)
table(pData(cds)$Cluster,pData(cds)$celltype)
plot_cell_clusters(cds, 1, 2 )
plot_cell_clusters(cds, 1, 2 ,color_by = "celltype")

pData(cds)$Cluster=pData(cds)$celltype
table(pData(cds)$Cluster)

Sys.time()
diff_test_res_mouse <- differentialGeneTest(cds,
                                            fullModelFormulaStr = "~Cluster")
Sys.time()
# Select genes that are significant at an FDR < 10%
sig_genes <- subset(diff_test_res_mouse, qval < 0.01)
sig_genes=sig_genes[order(sig_genes$pval),]
head(sig_genes[,c("gene_short_name", "pval", "qval")] ) 
cg=as.character(head(sig_genes$gene_short_name,n=10)) 

plot_genes_jitter(cds[cg,],
                  grouping = "Cluster",
                  color_by = "Cluster",
                  nrow= 3,
                  ncol = NULL )
cg2=as.character(tail(sig_genes$gene_short_name,n=10)) 
plot_genes_jitter(cds[cg2,],
                  grouping = "Cluster",
                  color_by = "Cluster",
                  nrow= 3,
                  ncol = NULL )

# step1
ordering_genes <- row.names (subset(diff_test_res_mouse, qval < 1e-150))#20;150;160;220
ordering_genes 
cds <- setOrderingFilter(cds, ordering_genes)
plot_ordering_genes(cds)
# step2
# we used defaul setting “DDRTree”
cds <- reduceDimension(cds, max_components = 2,
                       method = 'DDRTree')
# step3
cds <- orderCells(cds)
# 最后两个可视化函数，对结果进行可视化
plot_cell_trajectory(cds, color_by = "Cluster")  

p2=plot_cell_trajectory(cds, color_by = "Pseudotime",cell_size = 1)  
p2
length(cg)
plot_genes_in_pseudotime(cds[cg,],
                         color_by = "Pseudotime") 

plot_genes_in_pseudotime(cds[cg,],
                         color_by = "Cluster")
# https://davetang.org/muse/2017/10/01/getting-started-monocle/

my_cds_subset=cds
# pseudotime is now a column in the phenotypic data as well as the cell state
head(pData(my_cds_subset))

my_pseudotime_de <- differentialGeneTest(my_cds_subset,
                                         fullModelFormulaStr = "~sm.ns(Pseudotime)",
                                         cores = 1)
head(my_pseudotime_de)
p1=plot_cell_clusters(cds, 1, 2, color = c("Pseudotime"))
p2=plot_cell_clusters(cds, 1, 2, color = "celltype")
p3=plot_cell_clusters(cds, 1, 2, color = "orig.ident")
p1+p2+p3
save( my_cds_subset,my_pseudotime_de,
      file = 'mouse_hc_monocle2.Rdata')
saveRDS(my_cds_subset,file = 'mouse_hc_monocle2.rds')
################################################################################################
################################################################################################
################################################################################################
library(dplyr)
cds=mouse_hc_monocle2
cds=my_cds_subset
colnames(pData(cds))
table(pData(cds)$State,pData(cds)$Cluster)
library(ggsci)
p1=plot_cell_trajectory(cds, color_by = "Cluster",cell_size = 1)  + scale_color_nejm()
p1


p2=plot_cell_trajectory(cds, color_by = "Pseudotime",cell_size = 1)  
p2

p3=plot_cell_trajectory(cds, color_by = "State",cell_size = 1)  + scale_color_npg()
p3
p4=plot_cell_trajectory(cds, color_by = "celltype",cell_size = 1)
p4
p5=plot_cell_trajectory(cds, color_by = "orig.ident",cell_size = 1)
p5
library(patchwork)
p2+p4
p1+p2+p3+p4+p5
plot_cell_trajectory(cds, markers_linear=F, show_branch_points=F, color_by = "Pseudotime")  + 
  scale_color_gradient(low="grey", high="sienna2")

plot_cell_trajectory(cds, markers = c("Sox2", "Atoh1","Otof","Tbx2",
                                      "Slc17a8","Slc26a5","Ikzf2","Fgfr3",
                                      "Tmc1","Sorbs2","Kdm5b"), 
                     markers_linear=F, use_color_gradient=T, 
                     show_branch_points=F,cell_size = 1)#"Magi3","Ankfn1"
markers.to.plot<-c(
  "Sox2", "Atoh1","Otof","Tbx2",
  "Slc17a8","Slc26a5","Ikzf2","Fgfr3",
  "Tmc1","Chrna10"
  
)
sce=sce_all
Idents(sce) <- "celltype"
markers.to.plot<-c(
  "Ankfn1","Magi3"
  
)
VlnPlot(sce,features = markers.to.plot,assay = "SCT")
Idents(sce) <- "orig.ident"
VlnPlot(sce,features = markers.to.plot,assay = "SCT")



p5=plot_cell_trajectory(cds, color_by = "celltype",cell_size = 2) +
  scale_color_manual(breaks = c("iOHC", "iIHC","OHC", "IHC"), 
                     values=c("#c8a485", "#9ac7ac", "#f7b741","#55efc8"))
p6=plot_cell_trajectory(cds, color_by = "Pseudotime",cell_size = 2)  
p5+p6
#Performing differential gene expression analysis using the pseudo-time variable 
#created in the construction of the manifold.
diff_test_res2 <- differentialGeneTest(cds, fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names2 <- rownames(diff_test_res2[head(order(diff_test_res2$qval),100),])

plot_genes_in_pseudotime(cds[sig_gene_names2[1:6],], 
                         ncol=3, color_by = "celltype",cell_size = 2) +
  scale_color_manual(breaks = c("iOHC", "iIHC","OHC", "IHC"), 
                     values=c("#c8a485", "#9ac7ac", "#f7b741","#55efc8"))
plot_genes_in_pseudotime(cds[c("Slc26a5", "Ocm", "Ikzf2",
                               "Slc17a8","Slc7a14","Tbx2","Kdm5b",
                               "Tmc1",
                               "Sox2","Insm1","Atoh1"),],
                         ncol=3, color_by = "celltype",cell_size = 2) +
  scale_color_manual(breaks = c("iOHC", "iIHC","OHC", "IHC"), 
                     values=c("#c8a485", "#9ac7ac", "#f7b741","#55efc8"))


Time_diff <- diff_test_res2[,c(5,2,3,4,1,6,7)] #把gene放前面，也可以不改
write.csv(Time_diff, "Time_diff_all_mouse_HC.csv", row.names = F)

Time_genes <- Time_diff %>% pull(gene_short_name) %>% as.character()

p=plot_pseudotime_heatmap(cds[Time_genes,], num_clusters=4, show_rownames=T,
                          return_heatmap=T)

#extract top100 genes
Time_genes <- top_n(Time_diff, n = 100, desc(qval)) %>% pull(gene_short_name) %>% as.character()
p=plot_pseudotime_heatmap(cds[Time_genes,], num_clusters=3, 
                          show_rownames=T, return_heatmap=T)

BEAM_res <- BEAM(cds, branch_point = 1, cores = 1)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
BEAM_res2 = subset(BEAM_res,qval < 1e-50)
plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res,
                                                  qval < 1e-50)),],
                            branch_point = 1,
                            num_clusters = 3,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T)
hc_genes <- row.names(subset(fData(cds),
                               gene_short_name %in% c("Slc26a5", "Ocm", "Ikzf2",
                                                      "Slc17a8","Slc7a14","Tbx2","Kdm5b",
                                                      "Tmc1",
                                                      "Sox2","Insm1","Atoh1")))
plot_genes_violin(cds[c("Slc26a5", "Ocm", "Ikzf2",
                        "Slc17a8","Slc7a14","Tbx2","Kdm5b",
                        "Tmc1",
                        "Sox2","Insm1","Atoh1"),],
                  grouping = "celltype",
                  color_by = "celltype",
                  cell_size = 3,
                  ncol = 2) + 
  theme_classic()
Idents(sce) <- "celltype"
FeaturePlot(sce,features = c("Slc26a5", "Ocm", "Ikzf2",
                             "Slc17a8","Slc7a14","Tbx2","Kdm5b",
                             "Tmc1",
                             "Sox2","Insm1","Atoh1"),
            reduction = "umap_harmony",label = TRUE,pt.size = 1)

plot_cell_clusters(cds, 1, 2, color = c("Pseudotime"))
plot_cell_clusters(cds, 1, 2, color = "celltype")
plot_cell_clusters(cds, 1, 2,markers = c("Slc26a5", "Ocm", "Ikzf2",
                                         "Slc17a8","Slc7a14","Tbx2","Kdm5b",
                                         "Tmc1",
                                         "Sox2","Insm1","Atoh1"))
saveRDS(cds,file = "mouse_hc_monocle2.rds")









