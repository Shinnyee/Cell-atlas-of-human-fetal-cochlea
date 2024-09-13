library(Seurat)
library(ggplot2)
library(sctransform)
library(tidyverse)
options(stringsAsFactors = FALSE)
library(dittoSeq)
CoE_integrated <- readRDS("F:/Chai Lab/PROJECT_INTEGRATION ACROSS SPECIES/Cochlea_Evo_new/HOMOLOGY OF COCHLEA EPIEHELIUM ACROSS SPECIES/CoE_integrated.rds")
sce=CoE_integrated
# extract haircell subgroup to reintegrated, for identifying subtype across species.
table(Idents(sce))
Idents(sce) <- "orig.ident"
sce <- subset(sce,idents = "human")
Idents(sce) <- "integrated_cluster"
hc = sce[,sce$integrated_cluster  %in% c('OHC','IHC')]
table(Idents(hc))

hc$orig.ident <- "human_GW23"# for simplicity, this data includes BOTH GW23/26 timepoints
saveRDS(hc,file = "human_GW23_hc.rds")
#####we used harmony to batch 3 timepoints human hc dataset
library(Seurat)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(sctransform)
library(tidyverse)
library(harmony)
options(stringsAsFactors = FALSE)
library(tidyverse)
human_GW23_hc <- readRDS("F:/Chai Lab/PROJECT_INTEGRATION ACROSS SPECIES/Cochlea_Evo_new/HC TRAJOTORY/human_GW23_hc.rds")
hc=human_GW23_hc
matrix_human_23=hc@assays$RNA@counts
meta_human_23=hc@meta.data
human_hc_23 <- CreateSeuratObject(counts = matrix_human_23)
human_hc_23
human_hc_23 <- AddMetaData(human_hc_23, metadata = meta_human_23)
human_GW15_hc <- readRDS("F:/Chai Lab/PROJECT_INTEGRATION ACROSS SPECIES/human fetal_rawdata_reference/rawdata_GW15/human_GW15_hc.rds")
matrix_human_15=human_GW15_hc@assays$RNA@counts
meta_human_15=human_GW15_hc@meta.data
human_GW17_hc <- readRDS("F:/Chai Lab/PROJECT_INTEGRATION ACROSS SPECIES/human fetal_rawdata_reference/rawdata_GW17/human_GW17_hc.rds")
matrix_human_17=human_GW17_hc@assays$RNA@counts
meta_human_17=human_GW17_hc@meta.data
human_hc_15 <- CreateSeuratObject(counts = matrix_human_15)
human_hc_15
human_hc_15 <- AddMetaData(human_hc_15, metadata = meta_human_15)
human_hc_17 <- CreateSeuratObject(counts = matrix_human_17)
human_hc_17
human_hc_17 <- AddMetaData(human_hc_17, metadata = meta_human_17)
human_hc_23$celltype <- human_hc_23$integrated_cluster
human_hc_23 <- PercentageFeatureSet(human_hc_23, pattern = "^MT-", col.name = "pMT")
sce_all <- merge(x = human_hc_15, 
                 y = c( human_hc_17,human_hc_23))

sce_all
DefaultAssay(sce_all) <- "RNA"
sce_all <- sce_all %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData()

sce_all <- RunPCA(sce_all, features = VariableFeatures(sce_all), npcs = 5)

sce_all_harmony <- RunHarmony(sce_all,group.by.vars = 'orig.ident',reduction = "pca",
                              dims.use = 1:4,assay.use = "RNA")
sce_all[["harmony"]] <- sce_all_harmony[["harmony"]]
sce_all <- RunUMAP(sce_all,dims = 1:3,
                   reduction = "harmony",reduction.name = "umap_harmony")

p1 <- DimPlot(human_hc_GW_ALL, reduction = "umap_harmony", group.by = "orig.ident",pt.size = 2) + 
  ggtitle("UMAP Harmony")
p1
p2 <- DimPlot(human_hc_GW_ALL, reduction = "umap_harmony",pt.size = 2,
              group.by = "celltype",repel = TRUE) + 
  ggtitle("UMAP Harmony")
p2
p3 <- DimPlot(sce_all, reduction = "umap_harmony",pt.size = 1,
              split.by = "celltype",repel = TRUE) + 
  ggtitle("UMAP Harmony")
p3
p1+p2
####umap time point##############################################
#UMAP species
sce=human_hc_GW_ALL
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
  scale_color_manual(values = c("royalblue1", "maroon4", "sienna2"))
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
sce_all=human_hc_GW_ALL
sce_all <- SCTransform(sce_all, verbose = T, conserve.memory = T,vars.to.regress = 
                         c("nCount_RNA","pMT"))
DefaultAssay(sce_all) <- "SCT"
markers.to.plot<-c(
  "SOX2","ATOH1","INSM1",#Early HC
  "FGFR3","SLC26A5","OCM","IKZF2",   
  "OTOF","SLC17A8","SLC7A14","ATP2A3","CALB2",  
  "MYO7A","MYO6","PCP4","TMC1",
  "KDM5B","ZNF462","TOX","NPAS3","TBX2"
  
)
Idents(sce_all) <- "celltype"
VlnPlot(sce_all,features = markers.to.plot)

DotPlot(sce_all, features = markers.to.plot,
        dot.scale = 8,cols  =c("white", "#ad9300")
        ,group.by = "celltype") + RotatedAxis()

FeaturePlot(sce_all,features = markers.to.plot,label = TRUE,pt.size = 2)
FeaturePlot(sce_all,features = c("SOX2","ATOH1","PCP4","INSM1"),label = TRUE,pt.size = 2)

FeaturePlot(sce_all,features = c("TMC1","SLC26A5","SLC17A8","ATP2B2"),label = TRUE,pt.size = 2)


saveRDS(sce_all, file = "human_hc_GW_ALL.rds")
sce_all=human_hc_GW_ALL
DimPlot(sce_all)
################################################################################################
################################################################################################
###################MONOCLE 3 ANALYSIS ##########################################################
################################################################################################

library(monocle3)
sce=sce_all
sce=human_hc_GW_ALL
Idents(sce) <- "celltype"
table(Idents(sce))

data <- GetAssayData(sce,assay = 'RNA',slot='counts')
cell_metadata <- sce@meta.data
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
genes_of_interest <- c("SLC26A5","SLC17A8","TMC1",
                       "SOX2","INSM1","ATOH1")
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
             label_leaves = FALSE,label_branch_points = FALSE,cell_size = 0.8)
p

# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds, time_bin="130-170"){
  cell_ids <- which(colData(cds)[, "celltype"] == time_bin)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))


plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 0.8)
################################################################################################
#########DE analysis in monocle3
genes <- c("SLC26A5", "OCM", "IKZF2","FGFR3",
           "SLC17A8","SLC7A14","TBX2","KDM5B",
           "TMC1","SOX2","INSM1","ATOH1"
           )
cds_subset <- cds[rowData(cds)$gene_short_name %in% genes,]

gene_fits <- fit_models(cds_subset, model_formula_str = "~pseudotime")
fit_coefs <- coefficient_table(gene_fits)
time_terms <- fit_coefs %>% filter(term == "pseudotime")
time_terms %>% filter (q_value < 0.05) %>%
  select(gene_short_name, term, q_value, estimate)
plot_genes_violin(cds_subset, group_cells_by="celltype", ncol=2) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
################################################################################################
################################################################################################
####################################MONOCLE 2 ANALYSIS##########################################
################################################################################################
sce=human_hc_GW_ALL
sce=sce_all
library(monocle)# to lower the version below 4.1
table(Idents(sce))
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
sc_cds <- detectGenes(sc_cds, min_expr = 1) 
sc_cds <- sc_cds[fData(sc_cds)$num_cells_expressed > 5, ]
sc_cds
cds <- sc_cds
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds) 

disp_table <- dispersionTable(cds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 1)
cds <- setOrderingFilter(cds, unsup_clustering_genes$gene_id)
plot_ordering_genes(cds) 
plot_pc_variance_explained(cds, return_all = F) # norm_method='log'

cds <- reduceDimension(cds, max_components = 2, num_dim = 6,
                       reduction_method = 'tSNE', verbose = T,perplexity=24)
cds <- clusterCells(cds, num_clusters =4) 
plot_cell_clusters(cds, 1, 2 )
plot_cell_clusters(cds, 1, 2 , color_by = "celltype")
table(pData(cds)$Cluster) 
colnames(pData(cds)) 

table(pData(cds)$Cluster,pData(cds)$new)
plot_cell_clusters(cds, 1, 2 )
save(cds,file = 'hc_human_input_cds.Rdata')

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
sig_genes <- subset(diff_test_res_mouse, qval < 0.1)
sig_genes=sig_genes[order(sig_genes$pval),]
head(sig_genes[,c("gene_short_name", "pval", "qval")] ) 
cg=as.character(head(sig_genes$gene_short_name,n=6)) 

plot_genes_jitter(cds[cg,],
                  grouping = "Cluster",
                  color_by = "Cluster",
                  nrow= 3,
                  ncol = NULL )
cg2=as.character(tail(sig_genes$gene_short_name,n=6)) 
plot_genes_jitter(cds[cg2,],
                  grouping = "Cluster",
                  color_by = "Cluster",
                  nrow= 3,
                  ncol = NULL )

# step1
ordering_genes <- row.names (subset(diff_test_res_mouse, pval < 0.00000000000000000000003))#0.000000000000000000001
ordering_genes #0.00000000000000000000003
cds <- setOrderingFilter(cds, ordering_genes)
plot_ordering_genes(cds)
#step2

cds <- reduceDimension(cds, max_components = 2,
                       method = 'DDRTree')
# step3
cds <- orderCells(cds)

plot_cell_trajectory(cds, color_by = "Cluster")  


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
plot_cell_clusters(cds, 1, 2, color = c("Pseudotime"))
plot_cell_clusters(cds, 1, 2, color = "celltype")
save( my_cds_subset,my_pseudotime_de,
      file = 'human_hc_monocle2.Rdata')

################################################################################################ 

library(dplyr)
cds=my_cds_subset
colnames(pData(cds))
table(pData(cds)$State,pData(cds)$Cluster)
library(ggsci)
p1=plot_cell_trajectory(cds, color_by = "Cluster",cell_size = 2)  + scale_color_nejm()
p1


p2=plot_cell_trajectory(cds, color_by = "Pseudotime",cell_size = 2)  
p2

p3=plot_cell_trajectory(cds, color_by = "State",cell_size = 2)  + scale_color_npg()
p3
p4=plot_cell_trajectory(cds, color_by = "celltype",cell_size = 2)
p4
p5=plot_cell_trajectory(cds, color_by = "orig.ident",cell_size = 2)
p5
library(patchwork)
p2+p4
p1+p2+p3+p4+p5
plot_cell_trajectory(cds, markers_linear=F, show_branch_points=F, color_by = "Pseudotime")  + 
  scale_color_gradient(low="grey", high="sienna2")

p5=plot_cell_trajectory(cds, markers = c("SLC26A5", "OCM", "IKZF2","FGFR3",
                                      "SLC17A8","SLC7A14","TBX2","KDM5B",
                                      "TMC1","SOX2","INSM1","ATOH1"), 
                     markers_linear=F, use_color_gradient=T, 
                     show_branch_points=F,cell_size = 2)
p5
plot_cell_trajectory(cds, markers = c("SOX2", "ATOH1","OTOF","TBX2",
                                         "SLC17A8","SLC26A5","IKZF2","FGFR3",
                                         "TMC1","MAGI3","ANKFN1"), 
                        markers_linear=F, use_color_gradient=T, 
                        show_branch_points=F,cell_size = 3)

markers.to.plot<-c(
  "SOX2", "ATOH1","OTOF","TBX2",
  "SLC17A8","SLC26A5","IKZF2","FGFR3",
  "TMC1","MAGI3","ANKFN1"
  
)
Idents(sce) <- "celltype"
VlnPlot(sce,features = markers.to.plot)


p6=plot_cell_trajectory(cds, color_by = "celltype",cell_size = 2) +
  scale_color_manual(breaks = c("iOHC", "iIHC","OHC", "IHC"), 
                     values=c("#c8a485", "#9ac7ac", "#f7b741","#55efc8"))
p7=plot_cell_trajectory(cds, color_by = "Pseudotime",cell_size = 2)  
p6+p7

#Performing differential gene expression analysis using the pseudo-time variable 
#created in the construction of the manifold.
diff_test_res2 <- differentialGeneTest(cds, fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names2 <- rownames(diff_test_res2[head(order(diff_test_res2$qval),100),])

plot_genes_in_pseudotime(cds[sig_gene_names2[1:6],], 
                         ncol=3, color_by = "celltype",cell_size = 2) +
  scale_color_manual(breaks = c("iOHC", "iIHC","OHC", "IHC"), 
                     values=c("#c8a485", "#9ac7ac", "#f7b741","#55efc8"))
plot_genes_in_pseudotime(cds[c("SLC26A5", "OCM", "IKZF2",
                               "SLC17A8","SLC7A14","TBX2","KDM5B",
                               "TMC1",
                               "SOX2","INSM1","ATOH1"),],
                         ncol=3, color_by = "celltype",cell_size = 2) +
  scale_color_manual(breaks = c("iOHC", "iIHC","OHC", "IHC"), 
                     values=c("#c8a485", "#9ac7ac", "#f7b741","#55efc8"))


Time_diff <- diff_test_res2[,c(5,2,3,4,1,6,7)] #把gene放前面，也可以不改
write.csv(Time_diff, "Time_diff_all_HUMAN_HC.csv", row.names = F)

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
plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res,
                                                  pval < 0.05)),],
                            branch_point = 1,
                            num_clusters = 4,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T)
hc_genes <- row.names(subset(fData(cds),
                               gene_short_name %in% c("SLC26A5", "OCM", "IKZF2",
                                                      "SLC17A8","SLC7A14","TBX2","KDM5B",
                                                      "TMC1","OTOF",
                                                      "SOX2","INSM1","ATOH1")))
plot_genes_violin(cds[c("SLC26A5", "OCM", "IKZF2",
                        "SLC17A8","SLC7A14","TBX2","KDM5B",
                        "TMC1",
                        "SOX2","INSM1","ATOH1"),],
                  grouping = "celltype",
                  color_by = "celltype",
                  cell_size = 3,
                  ncol = 2) + 
  theme_classic()
Idents(sce) <- "celltype"
FeaturePlot(sce,features = c("SLC26A5", "OCM", "IKZF2",
                            "SLC17A8","SLC7A14","TBX2","KDM5B",
                            "TMC1",
                            "SOX2","INSM1","ATOH1"),
            reduction = "umap_harmony",label = TRUE,pt.size = 1)

plot_cell_clusters(cds, 1, 2, color = c("Pseudotime"))
plot_cell_clusters(cds, 1, 2, color = "celltype")
plot_cell_clusters(cds, 1, 2,markers = c("SLC26A5", "OCM", "IKZF2",
                                         "SLC17A8","SLC7A14","TBX2","KDM5B",
                                         "TMC1","OTOF",
                                         "SOX2","INSM1","ATOH1"))
Idents(sce) <- "celltype"
DefaultAssay(sce) <- "RNA"
markers.to.plot<-c(
  "ANKFN1","MAGI3","CHRNA10","SLC26A5", "OCM", "IKZF2",
  "SLC17A8","SLC7A14","TBX2","KDM5B",
  "TMC1","OTOF",
  "SOX2","INSM1","ATOH1","ZBTB20","CALB1"
  
)
VlnPlot(sce,features = markers.to.plot,assay = "SCT")
Idents(sce) <- "orig.ident"
VlnPlot(sce,features = markers.to.plot,assay = "SCT")









