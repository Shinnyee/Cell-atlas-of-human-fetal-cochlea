rm(list=ls())
#skip to 587 line
#remotes::install_github("zh542370159/SCP")
#remotes::install_version("ggplot2", version = "3.5.0")
#devtools::install_github("zhanghao-njmu/SCP")
#library(SCP)
#dev.off()
library(Seurat)
library(rtracklayer)
library(tibble)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(Hmisc)
library(cowplot)
library(zoo)
#remotes::install_version("ggplot2", version = "3.4.2")
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

#####################################################
#transfer H5AD into SEURAT OBJECT

library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(patchwork)
library(zellkonverter)
#BiocManager::install("zellkonverter")
#BiocManager::install("loomR",force = TRUE)
library(loomR)
library(Seurat)
library(SeuratDisk)
library(data.table)
adata_loom <- connect(filename = "human_sgn_scANVI_annotation.loom",
                      mode = "r+",skip.validate = TRUE)
matrix=adata_loom[["matrix"]][,]
matrix=t(matrix)
dim(matrix)
gene = adata_loom$row.attrs$var_names[]
barcode = adata_loom$col.attrs$obs_names[]

meta_data = read.csv('human_sgn_scANVI_annotation_obs.csv',row.names = 1) # as form as dataframe format
meta_feature = read.csv('human_sgn_scANVI_annotation_var.csv',row.names = 1)

colnames(matrix)= barcode
row.names(matrix)= gene
x_scvi = adata_loom$col.attrs$X_scVI[,]
x_umap = adata_loom$col.attrs$X_umap[,]
x_scanvi = adata_loom$col.attrs$X_scANVI[,]
x_pca=adata_loom$col.attrs$X_pca[,]
x_harmony=adata_loom$col.attrs$X_harmony[,]
#x_scvi = adata_loom$col.attrs$X_scVI[,]


seurat_object= CreateSeuratObject(counts = matrix,meta.data = meta_data,
                                  project = 'human_sgn_loom',
                                  min.cells = 0, 
                                  min.features = 0)
seurat_object@assays[["RNA"]]@meta.features <- meta_feature
x_scvi = t(x_scvi)
x_umap = t(x_umap)
x_scanvi = t(x_scanvi)
x_pca = t(x_pca)
rownames(x_scvi) = barcode
rownames(x_umap) = barcode
rownames(x_scanvi) = barcode
rownames(x_pca) = barcode
colnames(x_scvi) = c("scVI_1","scVI_2","scVI_3","scVI_4","scVI_5","scVI_6","scVI_7","scVI_8","scVI_9","scVI_10","scVI_11","scVI_12","scVI_13","scVI_14","scVI_15","scVI_16","scVI_17","scVI_18","scVI_19","scVI_20","scVI_21","scVI_22","scVI_23","scVI_24","scVI_25","scVI_26","scVI_27","scVI_28","scVI_29","scVI_30")
colnames(x_umap) = c('UMAP_1','UMAP_2')
colnames(x_scanvi) = c("scANVI_1","scANVI_2","scANVI_3","scANVI_4","scANVI_5","scANVI_6","scANVI_7","scANVI_8","scANVI_9","scANVI_10","scANVI_11","scANVI_12","scANVI_13","scANVI_14","scANVI_15","scANVI_16","scANVI_17","scANVI_18","scANVI_19","scANVI_20","scANVI_21","scANVI_22","scANVI_23","scANVI_24","scANVI_25","scANVI_26","scANVI_27","scANVI_28","scANVI_29","scANVI_30")
colnames(x_pca) = c('PCA_1','PCA_2','PCA_3','PCA_4','PCA_5','PCA_6','PCA_7','PCA_8',
                    'PCA_9','PCA_10','PCA_11','PCA_12','PCA_13','PCA_14','PCA_15','PCA_16',
                    'PCA_17','PCA_18','PCA_19','PCA_20','PCA_21','PCA_22','PCA_23','PCA_24',
                    'PCA_25','PCA_26','PCA_27','PCA_28','PCA_29','PCA_30','PCA_31','PCA_32',
                    'PCA_33','PCA_34','PCA_35','PCA_36','PCA_37','PCA_38','PCA_39','PCA_40',
                    'PCA_41','PCA_42','PCA_43','PCA_44','PCA_45','PCA_46','PCA_47','PCA_48',
                    'PCA_49','PCA_50')

sce1 = as.SingleCellExperiment(seurat_object)

sce1@int_colData@listData[["reducedDims"]]@listData[["scVI"]] = x_scvi
sce1@int_colData@listData[["reducedDims"]]@listData[["umap"]] = x_umap
sce1@int_colData@listData[["reducedDims"]]@listData[["scANVI"]] = x_scanvi
sce1@int_colData@listData[["reducedDims"]]@listData[["PCA"]] = x_pca
sgn_seurat <- as.Seurat(sce1)
DimPlot(sgn_seurat,reduction = "umap",group.by = "cell_type_final")
DimPlot(sgn_seurat,reduction = "umap",group.by = "batch")

my_cols <- c("#a056d3","#e1ea8b",
             "#118913",
             "#ce4a1b",
             "#6892bc",
             
             "#25b74e",
             "#62edd0", "#c4876e", "#ba1c1c"
             
)

values = c("#9b5792", "#575add", "#74c5e9","#62edd0","#1db342", "#b7de27", "#e8ba69")
my_cols <- c("#9b5792", "#575add", "#74c5e9",
             "#62edd0","#1db342", "#b7de27", 
             "#e8ba69","#a056d3","#ba1c1c"
             
)
my_cols2 <- my_cols[order(as.integer(names(my_cols)))]
scales::show_col(my_cols)
DimPlot(human_sgn_for_slinghsot_analysis,
        cols = my_cols, label=FALSE , repel=TRUE,reduction = "umap",pt.size = 2
        ,group.by = "cell_type_final_v2")#
saveRDS(sgn_seurat,file = "human_sgn_for_slinghsot_analysis.rds")
######################################################################################
##########################SLINGSHOT ANALYSIS########################################
#####################################################################################
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
## Read RDS
## ---------------------------------- #
human_sgn_for_slinghsot_analysis <- readRDS("F:/PROJECTS/PROJECT_HUMAN_FETAL_COCHLEAE/WORKPLACE/R/sc_pipeline/SGN analysis/human_sgn_for_slinghsot_analysis.rds")
sce=human_sgn_for_slinghsot_analysis
table(sce$cell_type_final)
#############
Idents(sce) <- "cell_type_final"
Idents(sce)
sce$cell_type_final <- Idents(sce)
table(sce$cell_type_final)


Idents(sce) <- "cell_type_final"
Idents(sce)
sce$timepoint =sce$batch
table(sce$timepoint)
##############
table(sce$timepoint,sce$cell_type_final)
include_celltype <- c("CLU+_mature_sgn",
                      "ALDH1A3+_immature_sgn" ,
                      "RUNX1+_mature_sgn" ,
                      "LYPD1+_immature_sgn",
                      "PRPH+_immature_sgn",
                      "NEUROD1+_early_sgn",
                  
                      "IGFBP2+_intermediate_sgn"
                      ) #
sce_ALL <- sce

## ---------------------------------- #
real_colors <- c("CLU+_mature_sgn" = "#33A02C",
                 #"NEFH+_mature_sgn"  = "#B2DF8A",
                 "RUNX1+_mature_sgn" = "#55A1B1",
                 "LYPD1+_immature_sgn" = "#8DD3C7",
                 "PRPH+_immature_sgn" = "#56e8d7",
                 #"POU4F1+_immature_sgn"="#9c7da0",
                 "IGFBP2+_intermediate_sgn"="#692872",
                 "NEUROD1+_early_sgn"="#e55ad5",
                 "ALDH1A3+_immature_sgn"="#8548e5"
                 
                 
) 


umap_colors <- c("CLU+_mature_sgn" = "#e2b129",
                 #"NEFH+_mature_sgn"  = "#af9161",
                 "RUNX1+_mature_sgn" = "#2ce5e5",
                 "LYPD1+_immature_sgn" = "#1965B0",
                 "PRPH+_immature_sgn" = "#7BAFDE",
                 #"POU4F1+_immature_sgn"="#dddd76",
                 "IGFBP2+_intermediate_sgn"="#dd9136",
                 "NEUROD1+_early_sgn"="#FB8072",
                 "ALDH1A3+_immature_sgn"="#DC050C"
  
)
## Subset 
## ---------------------------------- #
sce <- subset(sce, cell_type_final %in% include_celltype )
root_cell <- "NEUROD1+_early_sgn"
Idents(sce) <- sce$cell_type_final %>% as.character()
as.character()
## colors
## ---------------------------------- #
palette = plasma(100)
palette_celltype = brewer.pal(n = 10, name = "Set3")

C <- palette_celltype[as.factor(sce$cell_type_final)]
names(C) <- as.factor(sce$cell_type_final)

color <- unique(C)
names(color) <- unique(names(C))
## run Slingshots 
## ---------------------------------- #
start.clus <- root_cell
reduction = 'umap' #
sds= slingshot(Embeddings(sce, reduction), clusterLabels = Idents(sce),
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
umap_colors[sce_ALL$cell_type_final] -> cell_color
names(cell_color) <- colnames(sce_ALL)
cell_color[cell_id] -> U_C
## Plot slingshot curve II by cell type
## ---------------------------------- #
pseudotime_orig <- pseudotime
sds_orig <- sds
R_TIME <- pseudotime_orig %>% as.data.frame() %>% dplyr::mutate("cell_id" = rownames(.))
R_META <- sce_ALL[["cell_type_final"]] %>% dplyr::mutate("cell_id" = rownames(.))
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
make_umap(R_UMAP, "cell_type_final", umap_colors)
#
ggplot(R_UMAP, aes(x=UMAP_1, y=UMAP_2, group = eval(parse(text="cell_type_final")), 
                   colour = eval(parse(text="cell_type_final")))) +
  geom_point(size=2, alpha=0.7) + 
  theme_light() + theme_classic() +
  scale_colour_manual(values=umap_colors) + 
  theme(legend.position = "none", panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_alpha(guide = 'none')
#
P1=ggplot(R_UMAP, aes(x=UMAP_1, y=UMAP_2, colour = eval(parse(text="Lineage1")) )) +
  geom_point(size=2, alpha=0.7) +
  theme_light() + theme_classic() +
  scale_colour_viridis_c(na.value="#D3D3D3", option = "C") +
  theme(legend.title = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_alpha(guide = 'none')
P2=DimPlot(sce,group.by = "cell_type_final")
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
make_umap_c(R_UMAP, "Lineage2")
make_umap_c(R_UMAP, "Lineage3")
#
## Plot slingshot curve by cell type
## ---------------------------------- #
pseudotime_orig <- pseudotime
sds_orig <- sds
#pdf(file = 'slingshot_curves.separate.colored.by.time.all.umap.NEW.pdf', width = 14, height=7)

par(mfrow = c(1, 3))
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
  
  print(plot(sds_all$reducedDim, col = new_colors$color, pch = 16, cex =1,
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
df <- data.frame(sce[["cell_type_final"]], sce[["Lineage3"]]) 
colnames(df) <- c("cell_type_final", "Lineage")
na.omit(df) -> df
ggplot(df, aes(x=Lineage, fill=cell_type_final)) +
  geom_density(alpha=0.4) + theme_classic()+
  scale_fill_manual(values=C)

dev.off()

#
############################################################lineage1####################
L <- sce[["Lineage1"]] %>% deframe()
names(L) <- rownames(sce[["Lineage1"]])
L[!is.na(L)] %>% names() -> L2_cell
#R$Lineage2[!is.na(R$Lineage2)] %>% names() -> L2_cell
sce[, L2_cell] -> new_sce 
Idents(new_sce) <- new_sce$cell_type_final %>% as.character()
table(new_sce$cell_type_final)
#remove certain cell_types that not merely concentrate in lineage 1.
new_sce <- subset(new_sce, idents=c("RUNX1+_mature_sgn","PRPH+_immature_sgn"), inver=TRUE)
table(sce$cell_type_final)
# lineage 1 in new_sce, we remove sparse iIHC
#new_sce <- subset(new_sce, idents="iIHC", inver=TRUE)
# lineage 2 in new_sce, we remove sparse iOHC
#new_sce <- subset(new_sce, idents="iOHC", inver=TRUE)
Idents(new_sce)
table(Idents(new_sce))
# Transfer into cds
# ---------------------------------- #
cds <- as.CellDataSet(new_sce)
# Estimate size factor
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)#

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
                                             fullModelFormulaStr = '~cell_type_final',
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
write.csv(diff_test_res_filter,file = "slingshot.DEgenes.lineage1.csv")
sig_gene_names <- row.names(diff_test_res_filter)

head(sig_gene_names)
a <- as.matrix(rownames(cds))
b <- as.matrix(sig_gene_names)
c <- merge(x=a,y=b,all=FALSE)
Time_genes <- top_n(c, n = 500) %>% pull(V1) %>% as.character()

plot_pseudotime_heatmap(cds[Time_genes,], num_clusters=3, 
                        show_rownames=T, return_heatmap=T)

dev.off()
############################################################################
Time_diff <- diff_test_res_filter[,c(5,3,4)] #

Time_diff=Time_diff[!duplicated(Time_diff$gene_short_name),]
Time_genes <- Time_diff %>% pull(gene_short_name) %>% as.character()

head(Time_genes[1:200])


p=plot_pseudotime_heatmap(cds[Time_genes[1:200],], num_clusters=6, 
                          show_rownames=T, return_heatmap=T)

p

dev.off()
########################################################################

#
lung_genes <- row.names(subset(fData(cds),gene_short_name %in% c("NEUROD1","IGFBP2","TUBB3","PRPH",
                                                                 "CLU","LYPD1","PVALB",'GRIA2')))
lung_genes_subset <- cds[c("NEUROD1","CLU")]
lung_genes_subset <- cds[lung_genes,]
plot_genes_in_pseudotime(lung_genes_subset,color_by = "Pseudotime",ncol = 2)#celltype

plot_genes_in_pseudotime(lung_genes_subset,color_by = "timepoint",ncol = 2,
                         cell_size = 1.5)+
  scale_color_manual( 
    values=viridis(9))
plot_genes_in_pseudotime(lung_genes_subset,color_by = "cell_type_final",ncol = 2,
                         
                         cell_size = 1.5)+
  scale_color_manual(values=viridis(5))
  
dev.off()
####################################################################################

source("trajectory_Slingshot_to_Monocle.R")
# make plots
# ---------------------------------- #
hm <- get_pseudotime_matrix(cds[Time_genes[1:200],],  
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
                                           #add_annotation_col = "HCtype",
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
GENE_OF_INTEREST <- c("NEUROD1","IGFBP2","VIM","TUBB3","PRPH",
                      "EYA2","PVALB","CLU","LYPD1")

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
                                           num_clusters = 1,
                                           hmcols = my_color_set[sub_color][[1]],
                                           scale_max = 3, 
                                           scale_min = -3,
                                           cores = 1,
                                           show_rownames = T,
                                           return_heatmap = FALSE))
    
  }
  
  
}
dev.off()  


######################################################################################################
########################################################################################################
#######################################################################################################
##########===========================VERSION 2==========================================#############################
######################################################################################################
########################################################################################################
#transfer H5AD into SEURAT OBJECT

library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(patchwork)
library(zellkonverter)
#BiocManager::install("zellkonverter")
#BiocManager::install("loomR",force = TRUE)
library(loomR)
library(Seurat)
library(SeuratDisk)
library(data.table)
adata_loom <- connect(filename = "human_sgn_scANVI_annotation_V2.loom",
                      mode = "r+",skip.validate = TRUE)
matrix=adata_loom[["matrix"]][,]
matrix=t(matrix)
dim(matrix)
gene = adata_loom$row.attrs$var_names[]
barcode = adata_loom$col.attrs$obs_names[]

meta_data = read.csv('human_sgn_scANVI_annotation_obs_V2.csv',row.names = 1) # as form as dataframe format
meta_feature = read.csv('human_sgn_scANVI_annotation_var_V2.csv',row.names = 1)

colnames(matrix)= barcode
row.names(matrix)= gene
x_scvi = adata_loom$col.attrs$X_scVI[,]
x_umap = adata_loom$col.attrs$X_umap[,]
x_scanvi = adata_loom$col.attrs$X_scANVI[,]
x_pca=adata_loom$col.attrs$X_pca[,]
x_harmony=adata_loom$col.attrs$X_harmony[,]
#x_scvi = adata_loom$col.attrs$X_scVI[,]


seurat_object= CreateSeuratObject(counts = matrix,meta.data = meta_data,
                                  project = 'human_sgn_loom',
                                  min.cells = 0, 
                                  min.features = 0)
seurat_object@assays[["RNA"]]@meta.features <- meta_feature
x_scvi = t(x_scvi)
x_umap = t(x_umap)
x_scanvi = t(x_scanvi)
x_pca = t(x_pca)
rownames(x_scvi) = barcode
rownames(x_umap) = barcode
rownames(x_scanvi) = barcode
rownames(x_pca) = barcode
colnames(x_scvi) = c("scVI_1","scVI_2","scVI_3","scVI_4","scVI_5","scVI_6","scVI_7","scVI_8","scVI_9","scVI_10","scVI_11","scVI_12","scVI_13","scVI_14","scVI_15","scVI_16","scVI_17","scVI_18","scVI_19","scVI_20","scVI_21","scVI_22","scVI_23","scVI_24","scVI_25","scVI_26","scVI_27","scVI_28","scVI_29","scVI_30")
colnames(x_umap) = c('UMAP_1','UMAP_2')
colnames(x_scanvi) = c("scANVI_1","scANVI_2","scANVI_3","scANVI_4","scANVI_5","scANVI_6","scANVI_7","scANVI_8","scANVI_9","scANVI_10","scANVI_11","scANVI_12","scANVI_13","scANVI_14","scANVI_15","scANVI_16","scANVI_17","scANVI_18","scANVI_19","scANVI_20","scANVI_21","scANVI_22","scANVI_23","scANVI_24","scANVI_25","scANVI_26","scANVI_27","scANVI_28","scANVI_29","scANVI_30")
colnames(x_pca) = c('PCA_1','PCA_2','PCA_3','PCA_4','PCA_5','PCA_6','PCA_7','PCA_8',
                    'PCA_9','PCA_10','PCA_11','PCA_12','PCA_13','PCA_14','PCA_15','PCA_16',
                    'PCA_17','PCA_18','PCA_19','PCA_20','PCA_21','PCA_22','PCA_23','PCA_24',
                    'PCA_25','PCA_26','PCA_27','PCA_28','PCA_29','PCA_30','PCA_31','PCA_32',
                    'PCA_33','PCA_34','PCA_35','PCA_36','PCA_37','PCA_38','PCA_39','PCA_40',
                    'PCA_41','PCA_42','PCA_43','PCA_44','PCA_45','PCA_46','PCA_47','PCA_48',
                    'PCA_49','PCA_50')

sce1 = as.SingleCellExperiment(seurat_object)

sce1@int_colData@listData[["reducedDims"]]@listData[["scVI"]] = x_scvi
sce1@int_colData@listData[["reducedDims"]]@listData[["umap"]] = x_umap
sce1@int_colData@listData[["reducedDims"]]@listData[["scANVI"]] = x_scanvi
sce1@int_colData@listData[["reducedDims"]]@listData[["PCA"]] = x_pca
sgn_seurat <- as.Seurat(sce1)
DimPlot(sgn_seurat,reduction = "umap",group.by = "cell_type_final_v2")
FeaturePlot(sgn_seurat,features = c("NEUROD1",'PRPH','PBX3',
                                    "ESRRG",'PCDH9','OPCML',
                                    "PROX1",'PVALB','NEFH'))

my_cols <- c("#a056d3","#e1ea8b",
             "#118913",
             "#ce4a1b",
             "#6892bc",
             
             "#25b74e",
             "#62edd0", "#c4876e", "#ba1c1c"
             
)

values = c("#9b5792", "#575add", "#74c5e9","#62edd0","#1db342", "#b7de27", "#e8ba69")
my_cols <- c("#9b5792", "#575add", "#74c5e9",
             "#62edd0","#1db342", "#b7de27", 
             "#e8ba69","#a056d3","#ba1c1c"
             
)
my_cols2 <- my_cols[order(as.integer(names(my_cols)))]
scales::show_col(my_cols)
DimPlot(sgn_seurat,
        cols = my_cols, label=FALSE , repel=TRUE,reduction = "umap",pt.size = 2
        ,group.by = "cell_type_final_v2")#
saveRDS(sgn_seurat,file = "human_sgn_for_slinghsot_analysis_v2.rds")
######################################################################################
##########################SLINGSHOT ANALYSIS########################################
#####################################################################################
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
## Read RDS
## ---------------------------------- #
human_sgn_for_slinghsot_analysis_v2 <- readRDS("F:/PROJECTS/PROJECT_HUMAN_FETAL_COCHLEAE/WORKPLACE/R/sc_pipeline/SGN analysis/human_sgn_for_slinghsot_analysis_v2.rds")
sce=human_sgn_for_slinghsot_analysis_v2
# 
mito_genes <- grep("^MT-", rownames(sce), value = TRUE)
# 
ribo_genes <- grep("^RPS|^RPL", rownames(sce), value = TRUE)
# 
blood_genes <- c("HBA1", "HBA2", "HBB", "CD3D",
                 "CD3E", "CD4", "CD8A", "CD14", "CD19", "CD20", "CD34")

genes_to_remove <- unique(c(mito_genes, ribo_genes, blood_genes))

gene_list <- read.csv("mouse_to_human_genes.csv")
gene_list <- gene_list$HGNC.symbol[gene_list$HGNC.symbol %in% rownames(sce)]
sce_filtered <- subset(sce, features = gene_list)
sce_filtered <- subset(sce_filtered, features = setdiff(rownames(sce_filtered), genes_to_remove))
sce=sce_filtered
table(sce$cell_type_final_v2)
#############
Idents(sce) <- "cell_type_final_v2"
Idents(sce)
#sce$cell_type_final_v2 <- Idents(sce)
table(sce$cell_type_final_v2)


Idents(sce) <- "cell_type_final_v2"
Idents(sce)
sce$timepoint =sce$batch
table(sce$timepoint)
##############
table(sce$timepoint,sce$cell_type_final_v2)
include_celltype <- c("NEUROD1+/CDH2_early_sgn",
                      "GATA3+/IGFBP2+_intermediate_sgn" ,
                      "SLIT3+/LYPD1+_immature_sgn" ,
                      "ROR2+/GPC6+_immature_sgn",
                      "NTNG1+/GRM8+_immature_sgn",
                      "RUNX1+/CAMK2B_mature_sgn",
                      
                      "CLU+/SYT11+_mature_sgn"
) #
sce_ALL <- sce
## 
## ---------------------------------- #
real_colors <- c("NEUROD1+/CDH2_early_sgn" = "#33A02C",
                 #"NEFH+_mature_sgn"  = "#B2DF8A",
                 "GATA3+/IGFBP2+_intermediate_sgn" = "#55A1B1",
                 "SLIT3+/LYPD1+_immature_sgn" = "#8DD3C7",
                 "ROR2+/GPC6+_immature_sgn" = "#56e8d7",
                 #"POU4F1+_immature_sgn"="#9c7da0",
                 "NTNG1+/GRM8+_immature_sgn"="#692872",
                 "RUNX1+/CAMK2B_mature_sgn"="#e55ad5",
                 "CLU+/SYT11+_mature_sgn"="#8548e5"
                 
                 
) 


umap_colors <- c("NEUROD1+/CDH2_early_sgn" = "#e2b129",
                 #"NEFH+_mature_sgn"  = "#af9161",
                 "GATA3+/IGFBP2+_intermediate_sgn" = "#2ce5e5",
                 "SLIT3+/LYPD1+_immature_sgn" = "#1965B0",
                 "ROR2+/GPC6+_immature_sgn" = "#7BAFDE",
                 #"POU4F1+_immature_sgn"="#dddd76",
                 "NTNG1+/GRM8+_immature_sgn"="#dd9136",
                 "RUNX1+/CAMK2B_mature_sgn"="#FB8072",
                 "CLU+/SYT11+_mature_sgn"="#DC050C"
                 
)
## Subset 
## ---------------------------------- #
sce <- subset(sce, cell_type_final_v2 %in% include_celltype )
root_cell <- "NEUROD1+/CDH2_early_sgn"
Idents(sce) <- sce$cell_type_final_v2 %>% as.character()
as.character()
## colors
## ---------------------------------- #
palette = plasma(100)
palette_celltype = brewer.pal(n = 10, name = "Set3")

C <- palette_celltype[as.factor(sce$cell_type_final_v2)]
names(C) <- as.factor(sce$cell_type_final_v2)

color <- unique(C)
names(color) <- unique(names(C))
## run Slingshots 
## ---------------------------------- #
start.clus <- root_cell
reduction = 'umap' #
sds= slingshot(Embeddings(sce, reduction), clusterLabels = Idents(sce),
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
umap_colors[sce_ALL$cell_type_final_v2] -> cell_color
names(cell_color) <- colnames(sce_ALL)
cell_color[cell_id] -> U_C
## Plot slingshot curve II by cell type
## ---------------------------------- #
pseudotime_orig <- pseudotime
sds_orig <- sds
R_TIME <- pseudotime_orig %>% as.data.frame() %>% dplyr::mutate("cell_id" = rownames(.))
R_META <- sce_ALL[["cell_type_final_v2"]] %>% dplyr::mutate("cell_id" = rownames(.))
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
make_umap(R_UMAP, "cell_type_final_v2", umap_colors)
#
ggplot(R_UMAP, aes(x=UMAP_1, y=UMAP_2, group = eval(parse(text="cell_type_final_v2")), 
                   colour = eval(parse(text="cell_type_final_v2")))) +
  geom_point(size=2, alpha=0.7) + 
  theme_light() + theme_classic() +
  scale_colour_manual(values=umap_colors) + 
  theme(legend.position = "none", panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_alpha(guide = 'none')
#
P1=ggplot(R_UMAP, aes(x=UMAP_1, y=UMAP_2, colour = eval(parse(text="Lineage1")) )) +
  geom_point(size=2, alpha=0.7) +
  theme_light() + theme_classic() +
  scale_colour_viridis_c(na.value="#D3D3D3", option = "C") +
  theme(legend.title = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_alpha(guide = 'none')
P2=DimPlot(sce,group.by = "cell_type_final_v2")
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
make_umap_c(R_UMAP, "Lineage2")
make_umap_c(R_UMAP, "Lineage3")
#
## Plot slingshot curve by cell type
## ---------------------------------- #
pseudotime_orig <- pseudotime
sds_orig <- sds
#pdf(file = 'slingshot_curves.separate.colored.by.time.all.umap.NEW.pdf', width = 14, height=7)

par(mfrow = c(1, 3))
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
  
  print(plot(sds_all$reducedDim, col = new_colors$color, pch = 16, cex =1,
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
df <- data.frame(sce[["cell_type_final_v2"]], sce[["Lineage1"]]) 
colnames(df) <- c("cell_type_final_v2", "Lineage")
na.omit(df) -> df
ggplot(df, aes(x=Lineage, fill=cell_type_final_v2)) +
  geom_density(alpha=0.4) + theme_classic()+
  scale_fill_manual(values=C)

dev.off()

#
############################################################lineage1####################
L <- sce[["Lineage1"]] %>% deframe()
names(L) <- rownames(sce[["Lineage1"]])
L[!is.na(L)] %>% names() -> L2_cell
#R$Lineage2[!is.na(R$Lineage2)] %>% names() -> L2_cell
sce[, L2_cell] -> new_sce 
Idents(new_sce) <- new_sce$cell_type_final_v2 %>% as.character()
table(new_sce$cell_type_final_v2)
#remove certain cell_types that not merely concentrate in lineage 1.
new_sce <- subset(new_sce, idents=c("ROR2+/GPC6+_immature_sgn"), inver=TRUE)
table(sce$cell_type_final_v2)
# lineage 1 in new_sce, we remove sparse iIHC
#new_sce <- subset(new_sce, idents="iIHC", inver=TRUE)
# lineage 2 in new_sce, we remove sparse iOHC
#new_sce <- subset(new_sce, idents="iOHC", inver=TRUE)
Idents(new_sce)
table(Idents(new_sce))

df <- data.frame(new_sce[["cell_type_final_v2"]], new_sce[["Lineage1"]]) 
colnames(df) <- c("cell_type_final_v2", "Lineage")
na.omit(df) -> df
ggplot(df, aes(x=Lineage, fill=cell_type_final_v2)) +
  geom_density(alpha=0.4) + theme_classic()+
  scale_fill_manual(values=C)
#save pdf
dev.off()

# Transfer into cds
# ---------------------------------- #
cds <- as.CellDataSet(new_sce)
# Estimate size factor
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)#

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
                                             fullModelFormulaStr = '~cell_type_final_v2',
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
write.csv(diff_test_res_filter,file = "slingshot.DEgenes.lineage1.csv")
sig_gene_names <- row.names(diff_test_res_filter)

head(sig_gene_names)
a <- as.matrix(rownames(cds))
b <- as.matrix(sig_gene_names)
c <- merge(x=a,y=b,all=FALSE)
Time_genes <- top_n(c, n = 500) %>% pull(V1) %>% as.character()

plot_pseudotime_heatmap(cds[Time_genes,], num_clusters=3, 
                        show_rownames=T, return_heatmap=T)

dev.off()
############################################################################
Time_diff <- diff_test_res_filter[,c(5,3,4)] #

Time_diff=Time_diff[!duplicated(Time_diff$gene_short_name),]
Time_genes <- Time_diff %>% pull(gene_short_name) %>% as.character()

head(Time_genes[1:200])


p=plot_pseudotime_heatmap(cds[Time_genes[1:200],], num_clusters=6, 
                          show_rownames=T, return_heatmap=T)

p

dev.off()
########################################################################

#
lung_genes <- row.names(subset(fData(cds),gene_short_name %in% c("NEUROD1","IGFBP2","TUBB3","PRPH",
                                                                 "CLU","LYPD1","PVALB",'GRIA2',
                                                                 'RUNX1','ESRRG','PROX1','PBX3',
                                                                 'RORB','PCDH9')))
lung_genes_subset <- cds[c("NEUROD1","CLU")]
lung_genes_subset <- cds[lung_genes,]
plot_genes_in_pseudotime(lung_genes_subset,color_by = "Pseudotime",ncol = 2)#celltype

plot_genes_in_pseudotime(lung_genes_subset,color_by = "timepoint",ncol = 2,
                         cell_size = 1.5)+
  scale_color_manual( 
    values=viridis(9))
plot_genes_in_pseudotime(lung_genes_subset,color_by = "cell_type_final_v2",ncol = 2,
                         
                         cell_size = 1.5)+
  scale_color_manual(values=viridis(5))
lung_genes <- row.names(subset(fData(cds),gene_short_name %in% c("NEUROD1","PBX1","GATA3","EYA2",
                                                                 "RUNX1","PROX1")))
lung_genes_subset <- cds[lung_genes,]
plot_genes_in_pseudotime(lung_genes_subset,color_by = "cell_type_final_v2",ncol = 2,
                         
                         cell_size = 1.5)+
  scale_color_manual(values=viridis(5))
my_colors <- c("NEUROD1+/CDH2_early_sgn" = "#1db342", 
               "GATA3+/IGFBP2+_intermediate_sgn" = "#74c5e9", 
               "SLIT3+/LYPD1+_immature_sgn" = "#62edd0", 
               "NTNG1+/GRM8+_immature_sgn" = "#9b5792", 
               "RUNX1+/CAMK2B_mature_sgn" = "#e8ba69")

plot_genes_in_pseudotime(lung_genes_subset, 
                         color_by = "cell_type_final_v2",
                         ncol = 2,
                         cell_size = 1.5) +
  scale_color_manual(values = my_colors)
dev.off()
####################################################################################

source("trajectory_Slingshot_to_Monocle.R")
# make plots
# ---------------------------------- #
hm <- get_pseudotime_matrix(cds[Time_genes[1:200],],  
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
                                           #add_annotation_col = "HCtype",
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
GENE_OF_INTEREST <- c("NEUROD1","IGFBP2","TUBB3","PRPH",
                      "CLU","LYPD1","PVALB",'GRIA2',
                      'RUNX1','ESRRG','PROX1','PBX3',
                      'RORB','PCDH9')

print(plot_genes_in_pseudotime(cds[GENE_OF_INTEREST,], color_by = "Lineage1" ) +         
        scale_color_viridis(option = "C") 
      #scale_color_viridis()
)



print(plot_genes_in_pseudotime(cds[GENE_OF_INTEREST,], color_by = "Lineage1",ncol=2 ) +
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
                                           num_clusters = 1,
                                           hmcols = my_color_set[sub_color][[1]],
                                           scale_max = 3, 
                                           scale_min = -3,
                                           cores = 1,
                                           show_rownames = T,
                                           return_heatmap = FALSE))
    
  }
  
  
}
dev.off()  


############################################################lineage2####################
L <- sce[["Lineage2"]] %>% deframe()
names(L) <- rownames(sce[["Lineage2"]])
L[!is.na(L)] %>% names() -> L2_cell
#R$Lineage2[!is.na(R$Lineage2)] %>% names() -> L2_cell
sce[, L2_cell] -> new_sce 
Idents(new_sce) <- new_sce$cell_type_final_v2 %>% as.character()
table(new_sce$cell_type_final_v2)
#remove certain cell_types that not merely concentrate in lineage 1.
new_sce <- subset(new_sce, idents=c("ROR2+/GPC6+_immature_sgn","RUNX1+/CAMK2B_mature_sgn"), inver=TRUE)
table(sce$cell_type_final_v2)
# lineage 1 in new_sce, we remove sparse iIHC
#new_sce <- subset(new_sce, idents="iIHC", inver=TRUE)
# lineage 2 in new_sce, we remove sparse iOHC
#new_sce <- subset(new_sce, idents="iOHC", inver=TRUE)
Idents(new_sce)
table(Idents(new_sce))

df <- data.frame(new_sce[["cell_type_final_v2"]], new_sce[["Lineage2"]]) 
colnames(df) <- c("cell_type_final_v2", "Lineage")
na.omit(df) -> df
ggplot(df, aes(x=Lineage, fill=cell_type_final_v2)) +
  geom_density(alpha=0.4) + theme_classic()+
  scale_fill_manual(values=C)
#save pdf
dev.off()

# Transfer into cds
# ---------------------------------- #
cds <- as.CellDataSet(new_sce)
# Estimate size factor
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)#

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
                                             fullModelFormulaStr = '~cell_type_final_v2',
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

pData(cds)[["Lineage2"]] -> pData(cds)$Pseudotime
diff_test_res <- differentialGeneTest(cds[cds_ordering_genes,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime, df=3)", 
                                      cores = 1)

diff_test_res_filter=diff_test_res %>% filter(qval <0.01) %>% arrange(qval)
write.csv(diff_test_res_filter,file = "slingshot.DEgenes.lineage2.csv")
sig_gene_names <- row.names(diff_test_res_filter)

head(sig_gene_names)
a <- as.matrix(rownames(cds))
b <- as.matrix(sig_gene_names)
c <- merge(x=a,y=b,all=FALSE)
Time_genes <- top_n(c, n = 500) %>% pull(V1) %>% as.character()

plot_pseudotime_heatmap(cds[Time_genes,], num_clusters=3, 
                        show_rownames=T, return_heatmap=T)

dev.off()
############################################################################
Time_diff <- diff_test_res_filter[,c(5,3,4)] #

Time_diff=Time_diff[!duplicated(Time_diff$gene_short_name),]
Time_genes <- Time_diff %>% pull(gene_short_name) %>% as.character()

head(Time_genes[1:200])


p=plot_pseudotime_heatmap(cds[Time_genes[1:200],], num_clusters=6, 
                          show_rownames=T, return_heatmap=T)

p

dev.off()
########################################################################

#
lung_genes <- row.names(subset(fData(cds),gene_short_name %in% c("NEUROD1","IGFBP2","TUBB3","PRPH",
                                                                 "CLU","LYPD1","PVALB",'GRIA2',
                                                                 'RUNX1','ESRRG','PROX1','PBX3',
                                                                 'RORB','PCDH9','CLU',"SYT11")))
lung_genes_subset <- cds[c("NEUROD1","CLU")]
lung_genes_subset <- cds[lung_genes,]
plot_genes_in_pseudotime(lung_genes_subset,color_by = "Pseudotime",ncol = 2)#celltype

plot_genes_in_pseudotime(lung_genes_subset,color_by = "timepoint",ncol = 2,
                         cell_size = 1.5)+
  scale_color_manual( 
    values=viridis(9))
plot_genes_in_pseudotime(lung_genes_subset,color_by = "cell_type_final_v2",ncol = 2,
                         
                         cell_size = 1.5)+
  scale_color_manual(values=viridis(5))

dev.off()
lung_genes <- row.names(subset(fData(cds),gene_short_name %in% c("EYA2", "GATA3",
                                                                 "NEUROD1","PBX1",
                                                                 "RUNX1","ESRRG")))
lung_genes_subset <- cds[lung_genes,]
plot_genes_in_pseudotime(lung_genes_subset,color_by = "cell_type_final_v2",ncol = 2,
                         
                         cell_size = 1.5)+
  scale_color_manual(values=viridis(5))
plot_genes_in_pseudotime(lung_genes_subset,color_by = "Pseudotime",ncol = 2)#celltype
my_colors <- c("NEUROD1+/CDH2_early_sgn" = "#1db342", 
               "GATA3+/IGFBP2+_intermediate_sgn" = "#74c5e9", 
               "SLIT3+/LYPD1+_immature_sgn" = "#62edd0", 
               "NTNG1+/GRM8+_immature_sgn" = "#9b5792", 
               "CLU+/SYT11+_mature_sgn" = "#575add")

plot_genes_in_pseudotime(lung_genes_subset, 
                         color_by = "cell_type_final_v2",
                         ncol = 2,
                         cell_size = 1.5) +
  scale_color_manual(values = my_colors)


####################################################################################

source("trajectory_Slingshot_to_Monocle.R")
# make plots
# ---------------------------------- #
hm <- get_pseudotime_matrix(cds[Time_genes[1:200],],  
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
                                           #add_annotation_col = "HCtype",
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
GENE_OF_INTEREST <- c("EYA2", "GATA3",
                      "NEUROD1","PBX1",
                      "RUNX1","ESRRG")

print(plot_genes_in_pseudotime(cds[GENE_OF_INTEREST,], color_by = "Lineage2" ) +         
        scale_color_viridis(option = "C") 
      #scale_color_viridis()
)



print(plot_genes_in_pseudotime(cds[GENE_OF_INTEREST,], color_by = "Lineage2",ncol=2 ) +
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
                                           num_clusters = 1,
                                           hmcols = my_color_set[sub_color][[1]],
                                           scale_max = 3, 
                                           scale_min = -3,
                                           cores = 1,
                                           show_rownames = T,
                                           return_heatmap = FALSE))
    
  }
  
  
}
dev.off()  


############################################################lineage3####################
L <- sce[["Lineage3"]] %>% deframe()
names(L) <- rownames(sce[["Lineage3"]])
L[!is.na(L)] %>% names() -> L2_cell
#R$Lineage2[!is.na(R$Lineage2)] %>% names() -> L2_cell
sce[, L2_cell] -> new_sce 
Idents(new_sce) <- new_sce$cell_type_final_v2 %>% as.character()
table(new_sce$cell_type_final_v2)
#remove certain cell_types that not merely concentrate in lineage 1.
new_sce <- subset(new_sce, idents=c("NTNG1+/GRM8+_immature_sgn","RUNX1+/CAMK2B_mature_sgn"), inver=TRUE)
table(sce$cell_type_final_v2)
# lineage 1 in new_sce, we remove sparse iIHC
#new_sce <- subset(new_sce, idents="iIHC", inver=TRUE)
# lineage 2 in new_sce, we remove sparse iOHC
#new_sce <- subset(new_sce, idents="iOHC", inver=TRUE)
Idents(new_sce)
table(Idents(new_sce))

df <- data.frame(new_sce[["cell_type_final_v2"]], new_sce[["Lineage3"]]) 
colnames(df) <- c("cell_type_final_v2", "Lineage")
na.omit(df) -> df
ggplot(df, aes(x=Lineage, fill=cell_type_final_v2)) +
  geom_density(alpha=0.4) + theme_classic()+
  scale_fill_manual(values=C)
#save pdf
dev.off()

# Transfer into cds
# ---------------------------------- #
cds <- as.CellDataSet(new_sce)
# Estimate size factor
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)#

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
                                             fullModelFormulaStr = '~cell_type_final_v2',
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

pData(cds)[["Lineage3"]] -> pData(cds)$Pseudotime
diff_test_res <- differentialGeneTest(cds[cds_ordering_genes,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime, df=3)", 
                                      cores = 1)

diff_test_res_filter=diff_test_res %>% filter(qval <0.01) %>% arrange(qval)
write.csv(diff_test_res_filter,file = "slingshot.DEgenes.lineage3.csv")
sig_gene_names <- row.names(diff_test_res_filter)

head(sig_gene_names)
a <- as.matrix(rownames(cds))
b <- as.matrix(sig_gene_names)
c <- merge(x=a,y=b,all=FALSE)
Time_genes <- top_n(c, n = 500) %>% pull(V1) %>% as.character()

plot_pseudotime_heatmap(cds[Time_genes,], num_clusters=3, 
                        show_rownames=T, return_heatmap=T)

dev.off()
############################################################################
Time_diff <- diff_test_res_filter[,c(5,3,4)] #

Time_diff=Time_diff[!duplicated(Time_diff$gene_short_name),]
Time_genes <- Time_diff %>% pull(gene_short_name) %>% as.character()

head(Time_genes[1:200])


p=plot_pseudotime_heatmap(cds[Time_genes[1:200],], num_clusters=6, 
                          show_rownames=T, return_heatmap=T)

p

dev.off()
########################################################################

#
lung_genes <- row.names(subset(fData(cds),gene_short_name %in% c("NEUROD1","IGFBP2","TUBB3","PRPH",
                                                                 "CLU","LYPD1","PVALB",'GRIA2',
                                                                 'RUNX1','ESRRG','PROX1','PBX3',
                                                                 'RORB','PCDH9','CLU',"SYT11")))
lung_genes_subset <- cds[c("NEUROD1","CLU")]
lung_genes_subset <- cds[lung_genes,]
plot_genes_in_pseudotime(lung_genes_subset,color_by = "Pseudotime",ncol = 2)#celltype

plot_genes_in_pseudotime(lung_genes_subset,color_by = "timepoint",ncol = 2,
                         cell_size = 1.5)+
  scale_color_manual( 
    values=viridis(9))
plot_genes_in_pseudotime(lung_genes_subset,color_by = "cell_type_final_v2",ncol = 2,
                         
                         cell_size = 1.5)+
  scale_color_manual(values=viridis(5))
lung_genes <- row.names(subset(fData(cds),gene_short_name %in% c("EYA2", "GATA3",
                                                                 "NEUROD1","PBX1",
                                                                 "ESRRG",
                                                                 'PBX3')))
lung_genes_subset <- cds[lung_genes,]
plot_genes_in_pseudotime(lung_genes_subset,color_by = "cell_type_final_v2",ncol = 2,
                         
                         cell_size = 1.5)+
  scale_color_manual(values=viridis(5))
plot_genes_in_pseudotime(lung_genes_subset,color_by = "Pseudotime",ncol = 2)#celltype
my_colors <- c("NEUROD1+/CDH2_early_sgn" = "#1db342", 
               "GATA3+/IGFBP2+_intermediate_sgn" = "#74c5e9", 
               "SLIT3+/LYPD1+_immature_sgn" = "#62edd0", 
               "ROR2+/GPC6+_immature_sgn" = "#b7de27", 
               "CLU+/SYT11+_mature_sgn" = "#575add")

plot_genes_in_pseudotime(lung_genes_subset, 
                         color_by = "cell_type_final_v2",
                         ncol = 2,
                         cell_size = 1.5) +
  scale_color_manual(values = my_colors)

dev.off()
####################################################################################

source("trajectory_Slingshot_to_Monocle.R")
# make plots
# ---------------------------------- #
hm <- get_pseudotime_matrix(cds[Time_genes[1:200],],  
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
                                           #add_annotation_col = "HCtype",
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
GENE_OF_INTEREST <- c("NEUROD1","IGFBP2","TUBB3","PRPH",
                      "CLU","LYPD1","PVALB",'GRIA2',
                      'RUNX1','ESRRG','PROX1','PBX3',
                      'RORB','PCDH9',"SYT11")

print(plot_genes_in_pseudotime(cds[GENE_OF_INTEREST,], color_by = "Lineage3" ) +         
        scale_color_viridis(option = "C") 
      #scale_color_viridis()
)



print(plot_genes_in_pseudotime(cds[GENE_OF_INTEREST,], color_by = "Lineage3",ncol=2 ) +
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
                                           num_clusters = 1,
                                           hmcols = my_color_set[sub_color][[1]],
                                           scale_max = 3, 
                                           scale_min = -3,
                                           cores = 1,
                                           show_rownames = T,
                                           return_heatmap = FALSE))
    
  }
  
  
}
dev.off()  


