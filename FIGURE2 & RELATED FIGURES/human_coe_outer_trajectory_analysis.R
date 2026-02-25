rm(list=ls())
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
adata_loom <- connect(filename = "human_coe_outer.loom",
                      mode = "r+",skip.validate = TRUE)
matrix=adata_loom[["matrix"]][,]
matrix=t(matrix)
dim(matrix)
gene = adata_loom$row.attrs$var_names[]
barcode = adata_loom$col.attrs$obs_names[]

meta_data = read.csv('human_coe_outer_obs.csv',row.names = 1) # as form as dataframe format
meta_feature = read.csv('human_coe_outer_var.csv',row.names = 1)

colnames(matrix)= barcode
row.names(matrix)= gene
#x_scvi = adata_loom$col.attrs$X_scVI[,]
x_umap=read.csv("human_coe_outer_umap.csv",row.names = 1)
x_scANVI=read.csv("human_coe_outer_scANVI.csv",row.names = 1)

seurat_object= CreateSeuratObject(counts = matrix,meta.data = meta_data,
                                  project = 'human_cochlea_loom',
                                  min.cells = 0, 
                                  min.features = 0)
seurat_object@assays[["RNA"]]@meta.features <- meta_feature

rownames(x_umap) = barcode
colnames(x_umap) = c('UMAP_1','UMAP_2')
rownames(x_scANVI) = barcode
colnames(x_scANVI) = c('scANVI_1','scANVI_2','scANVI_3','scANVI_4','scANVI_5','scANVI_6',
                       'scANVI_7','scANVI_8','scANVI_9','scANVI_10','scANVI_11','scANVI_12',
                       'scANVI_13','scANVI_14','scANVI_15','scANVI_16','scANVI_17','scANVI_18',
                       'scANVI_19','scANVI_20','scANVI_21','scANVI_22','scANVI_23','scANVI_24',
                       'scANVI_25','scANVI_26','scANVI_27','scANVI_28','scANVI_29','scANVI_30')

seurat_object$cellid <- rownames(seurat_object@meta.data)
head(seurat_object$cellid)
rownames(x_umap)=seurat_object$cellid
rownames(x_scANVI)=seurat_object$cellid

seurat_object <-FindVariableFeatures(seurat_object)
seurat_object <- ScaleData(seurat_object)
seurat_object <- RunPCA(seurat_object)
seurat_object<- RunUMAP(seurat_object,dims = 1:30)
seurat_object <- FindNeighbors(seurat_object,dims = 1:30)
seurat_object <- FindClusters(seurat_object,resolution = 0.8)
DimPlot(seurat_object,reduction = "umap",group.by = "cell_type_final_v4")

seurat_object@reductions[["umap"]]@cell.embeddings=as.matrix(x_umap)
DimPlot(seurat_object,reduction = "umap",group.by = "cell_type_final_v4")

sce1 = as.SingleCellExperiment(seurat_object)
colnames(x_umap) = c('umap_1','umap_2')
sce1@int_colData@listData[["reducedDims"]]@listData[["umap"]] = x_umap
sce1@int_colData@listData[["reducedDims"]]@listData[["scANVI"]] = x_scANVI
seurat_object <- as.Seurat(sce1)
seurat_object <- RunUMAP(seurat_object,reduction = "scANVI",dims = 1:30)
DimPlot(seurat_object,reduction = "umap",group.by = "cell_type_final_v4")
DimPlot(seurat_object,reduction = "umap",group.by = "gw")
#devtools::install_version("ggplot2", version = "3.4.2")
#devtools::install_version("Matrix", version = "1.5.4")


#human_cochlea_seurat <- as.Seurat(sce1)
#DimPlot(human_cochlea_seurat,reduction = "umap",group.by = "age_bins")

#sce1 = as.SingleCellExperiment(human_cochlea_seurat)
#sce1@int_colData@listData[["reducedDims"]]@listData[["scVI"]] = x_scvi
#sce1@int_colData@listData[["reducedDims"]]@listData[["umap"]] = x_umap
##############################################################################
############################################################################
#MILOR TO TEST DIFERENTIAL ABUNDANCE BETWEEN CONTROL AND NEOMYCIN CONDITIONS

########################################################################################
#######################################################################################
sce=seurat_object
#sce <- SCTransform(sce)
Idents(sce) <- "cell_type_final_v4"
#rename cluster_label
Idents(sce) <- "cell_type_final_v4"
table(Idents(sce))
new.cluster.ids <- c("DC",'Lateral_PSD','PC',
                     'DC','Outer_HC ','HEC','PC','HEC'
)
names(new.cluster.ids) <- levels(sce)
sce<- RenameIdents(sce, new.cluster.ids)
sce$celltype <- Idents(sce)
Idents(sce) <- "celltype"
table(Idents(sce))

all.markers <- FindAllMarkers(sce, assay = "RNA", slot = "data", test.use = "roc")
write.csv(all.markers,file = "human_coe_outer_trajetory_roctest.csv")
all.markers <- all.markers[which(all.markers$avg_diff > 0), ]
top50 <- all.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_diff)
DoHeatmap(sce, features = top50$gene) + NoLegend()

sce_monocle2=sce

Idents(sce_monocle2) <- "celltype"


library(monocle)
table(Idents(sce_monocle2))
sample_ann <-  sce_monocle2@meta.data  
head(sample_ann)
gene_ann <- data.frame(
   gene_short_name = rownames(sce_monocle2@assays$RNA) , 
   row.names =  rownames(sce_monocle2@assays$RNA) 
)
head(gene_ann)
pd <- new("AnnotatedDataFrame",
          data=sample_ann)
fd <- new("AnnotatedDataFrame",
          data=gene_ann)
ct=as.data.frame(sce_monocle2@assays$RNA@counts)
ct[1:4,1:4]
sc_cds <- newCellDataSet(
   as.matrix(ct), 
   phenoData = pd,
   featureData =fd,
   expressionFamily = negbinomial.size(),
   lowerDetectionLimit=1)
sc_cds
sc_cds <- detectGenes(sc_cds, min_expr = 0.1) 
sc_cds <- sc_cds[fData(sc_cds)$num_cells_expressed > 5, ]
sc_cds
cds <- sc_cds
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds) 
# 并不是所有的基因都有作用，所以先进行挑选，合适的基因用来进行聚类。
disp_table <- dispersionTable(cds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
cds <- setOrderingFilter(cds, unsup_clustering_genes$gene_id)
plot_ordering_genes(cds) 
plot_pc_variance_explained(cds, return_all = F) # norm_method='log'
# 其中 num_dim 参数选择基于上面的PCA图
cds <- reduceDimension(cds, max_components = 2, num_dim = 10,
                       reduction_method = 'tSNE', verbose = T)
cds <- clusterCells(cds, num_clusters =8) 
plot_cell_clusters(cds, 1, 2 )
plot_cell_clusters(cds, 1, 2 , color_by = "celltype")
table(pData(cds)$Cluster) 
colnames(pData(cds)) 

table(pData(cds)$Cluster,pData(cds)$celltype)
plot_cell_clusters(cds, 1, 2 )
# 接下来很重要，到底是看哪个性状的轨迹

colnames(pData(cds))
table(pData(cds)$Cluster)
table(pData(cds)$Cluster,pData(cds)$celltype)
plot_cell_clusters(cds, 1, 2 )
plot_cell_clusters(cds, 1, 2 ,color_by = "celltype")
## 我们这里并不能使用 monocle的分群
# 还是依据前面的 seurat分群, 其实取决于自己真实的生物学意图
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
cg=as.character(head(sig_genes$gene_short_name,n=6)) 
#  挑选差异最显著的基因可视化
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
# 第一步: 挑选合适的基因. 有多个方法，例如提供已知的基因集，这里选取统计学显著的差异基因列表
dev.off()
Idents(sce_monocle2) <- "celltype"
#sce_monocle2 <- SCTransform(sce_monocle2, method = "glmGamPoi")
sce_monocle2 <- PrepSCTFindMarkers(sce_monocle2)
sce_monocle2 <- SCTransform(sce_monocle2)
cellmarker <- FindAllMarkers(sce_monocle2, assay = "SCT", slot = "data", test.use = "roc" )
#write.csv(cellmarker,file = "COE_integrated_clusters_DEGs.csv")
cellmarker <- cellmarker[which(cellmarker$avg_diff > 0), ]
ordering_genes <- cellmarker$gene[1:565]#565
ordering_genes <- unique(ordering_genes)
#ordering_genes <- row.names (subset(diff_test_res_mouse, pval < 1e-50))#20,50,

cds <- setOrderingFilter(cds, ordering_genes)
plot_ordering_genes(cds)
# 第二步: 降维。降维的目的是为了更好的展示数据。函数里提供了很多种方法,
# 不同方法的最后展示的图都不太一样, 其中“DDRTree”是Monocle2使用的默认方法
cds <- reduceDimension(cds, max_components = 2,
                       method = 'DDRTree',norm_method = "log")#,norm_method = "log"
# 第三步: 对细胞进行排序
cds <- orderCells(cds)
# 最后两个可视化函数，对结果进行可视化

library(paletteer)
d_palettes<- palettes_d_names
d_palettes
colcors_merge <- data.frame(d_palettes)
mycols_3 <- paletteer_d("awtools::mpalette", n = 9)
mycols_4 <- paletteer_d("palettetown::magmar", n = 12)
table(sce_monocle2$celltype)
p1=plot_cell_trajectory(cds, color_by = "celltype",cell_size = 1)+  
   scale_color_manual(breaks = c( "Lateral_PSD","DC",
                                  "PC","HEC","Outer_HC"
                                  
   ), 
   values=mycols_3)



p1

p2=plot_cell_trajectory(cds, markers_linear=F, show_branch_points=F, color_by = "Pseudotime",
                        cell_size = 1)  + 
   scale_color_gradient(low="grey", high="sienna2")
p2
p3=plot_cell_trajectory(cds, color_by = "gw",cell_size = 1) 


p3

table(sce_monocle2$cell_type_final_v4)
p4=plot_cell_trajectory(cds, color_by = "cell_type_final_v4",cell_size = 1)+  
   scale_color_manual(breaks = c( "EGFR+ELMO1+HEC","FGFR3+PROX1+Lateral_PSD",
                                  "FST+AFF3+HEC","LGR6+TMPRSS3+_PC",
                                  "Outer_HC","PPP1R2+HIST1H2AC+_PC",
                                  "PTGDS+CEMIP+DC","SERPINE2+JAG1+DC"
                                  
   ), 
   values=c("#9a9ae8","#cc9997","#5757dd","#843aa3",
            "#992b28","#a48baf","#38a3b2","#3d707c"
            
            
            
   ))
p4
p5=plot_cell_trajectory(cds, color_by = "State",cell_size = 1)
p5
table(sce_monocle2$sample)
mycols_5 <- paletteer_d("awtools::a_palette")


library(patchwork)
p1+p2+p3+p4




#WNT PATHWAY
plot_cell_trajectory(cds, markers = c("LGR5", "FGFR3","PROX1",
                                      "JAG1","HES5","OCM","SLC26A5"), 
                     markers_linear=F, use_color_gradient=T, 
                     show_branch_points=F,cell_size = 0.8)+
   scale_color_gradient2(low="#88f6f2", mid="#cbe0e0",high="#f77146")


saveRDS(cds,file="human_coe_outer.rds")
save.image("human_coe_outer.RData")


HC_genes <- row.names(subset(fData(cds),
                             gene_short_name %in% c("LGR5", "ATOH1", "SOX2")))
plot_genes_branched_pseudotime(cds[HC_genes,],
                               branch_point = 1,
                               color_by = "celltype",
                               ncol = 1)
colour=c("#9a9ae8","#cc9997","#5757dd","#843aa3",
         "#992b28","#a48baf","#38a3b2","#3d707c"
)
p2 <- plot_complex_cell_trajectory(cds, x = 1, y = 2,size=1.2,
                                   color_by = "cell_type_final_v4")+
   scale_color_manual(values = colour) +
   theme(legend.title = element_blank())
p2
p3 <- plot_complex_cell_trajectory(cds, x = 1, y = 2,
                                   color_by = "Pseudotime")

p3
p3 <- plot_complex_cell_trajectory(cds, x = 1, y = 2,
                                   color_by = "gw")

p3
#提取数据=======================================================================
data_df <- t(reducedDimS(cds)) %>% as.data.frame() %>% #提取坐标  
   select_(Component_1 = 1, Component_2 = 2) %>% #重命名  
   rownames_to_column("cells") %>% #rownames命名  
   mutate(pData(cds)$State) %>% #添加State  
   mutate(pData(cds)$Pseudotime,         
          pData(cds)$cell_type_final_v4,          
          pData(cds)$celltype)#将这些需要作图的有用信息都添加上
colnames(data_df) <- c("cells","Component_1","Component_2","State",                     
                       "Pseudotime","cell_type_final_v4","celltype")


#轨迹数据提取---完全摘录于monocle包原函数
dp_mst <- minSpanningTree(cds)
reduced_dim_coords <- reducedDimK(cds)
ica_space_df <- Matrix::t(reduced_dim_coords) %>% as.data.frame() %>%   
   select_(prin_graph_dim_1 = 1, prin_graph_dim_2 = 2) %>%   
   mutate(sample_name = rownames(.), sample_state = rownames(.))


#构建一个做轨迹线图的数据
edge_df <- dp_mst %>% igraph::as_data_frame() %>%   
   select_(source = "from", target = "to") %>%   
   left_join(ica_space_df %>% select_(source = "sample_name",                                     
                                      source_prin_graph_dim_1 = "prin_graph_dim_1",                                      
                                      source_prin_graph_dim_2 = "prin_graph_dim_2"), by = "source") %>%   
   left_join(ica_space_df %>% select_(target = "sample_name",                                      
                                      target_prin_graph_dim_1 = "prin_graph_dim_1",                                     
                                      target_prin_graph_dim_2 = "prin_graph_dim_2"), by = "target")



#计算细胞比例
Cellratio <- prop.table(table(data_df$State, data_df$cell_type_final_v4), margin = 2)#计算各组样本不同细胞群比例
Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio) <- c('State',"cell_type_final_v4","Freq")


#ggplot作图
library(ggplot2)
library(tidydr)
library(ggforce)
library(ggrastr)

g <- ggplot() +  
   geom_point_rast(data = data_df, aes(x = Component_1,                                 
                                       y = Component_2,                                
                                       color =Pseudotime)) + #散点图  
   scale_color_viridis()+ #密度色  
   geom_segment(aes_string(x = "source_prin_graph_dim_1",                          
                           y = "source_prin_graph_dim_2",                          
                           xend = "target_prin_graph_dim_1",                          
                           yend = "target_prin_graph_dim_2"),               
                linewidth = 1,                
                linetype = "solid", na.rm = TRUE, data = edge_df)+#添加轨迹线  
   theme_dr(arrow = grid::arrow(length = unit(0, "inches")))+#坐标轴主题修改  
   theme(panel.grid.major = element_blank(),        
         panel.grid.minor = element_blank())+  
   geom_arc(arrow = arrow(length = unit(0.15, "inches"), #曲线箭头                         
                          type = "closed",angle=30),           
            aes(x0=0,y0=-3,r=5, start=-0.4*pi, end=0.4*pi),lwd=1)+  
   geom_arc_bar(data=subset(Cellratio,State=='1'),stat = "pie",#添加饼图             
                aes(x0=-10,y0=0,r0=0,r=1.5,amount=Freq,fill=cell_type_final_v4))+ 
   geom_arc_bar(data=subset(Cellratio,State=='2'),stat = "pie",               
                aes(x0=5,y0=2,r0=0,r=1.5,amount=Freq,fill=cell_type_final_v4))+  
   geom_arc_bar(data=subset(Cellratio,State=='3'),stat = "pie",               
                aes(x0=-3,y0=8,r0=0,r=1.5,amount=Freq,fill=cell_type_final_v4))+  
   scale_fill_manual(values = c("#9a9ae8","#cc9997","#5757dd","#843aa3",
                                "#992b28","#a48baf","#38a3b2","#3d707c"))
g
########################################################################################
#######################################################################################
library(CytoTRACE)
monocle_meta <- data.frame(t(cds@reducedDimS),                          
                           cds$Pseudotime,                         
                           cds$State,                          
                           cds$cell_type_final_v4)
colnames(monocle_meta) <- c("C1", "C2", "Pseudotime", "State", "cell_type_final_v4")
phenot1 <- monocle_meta$cell_type_final_v4
phenot1 <- as.character(phenot1)
names(phenot1) <- rownames(monocle_meta)
emb_monocle<-monocle_meta[,1:2]
exp1 <- as.matrix(sce_monocle2@assays$RNA@counts)
exp1<-exp1[apply(exp1>0,1,sum)>=5,]
results<-CytoTRACE(exp1,ncores=1)
plotCytoTRACE(results, phenotype = phenot1, emb = emb_monocle)
########################################################################################
##########################################################################################
expressed_genes=row.names(subset(fData(cds),num_cells_expressed>=10)) #在部分基因里面找
pseudotime_de <- differentialGeneTest(cds[expressed_genes,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
pseudotime_de <- pseudotime_de[order(pseudotime_de$qval), ]
states_de <- differentialGeneTest(cds[expressed_genes,],
                                  fullModelFormulaStr = "~State")
states_de <- states_de[order(states_de$qval), ]


write.table(pseudotime_de, file = "pseudotime_de_coe_outer.rds", quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
write.table(states_de, file = "states_de_coe_outer.rds", quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)

#cell fate 2===>outer_HC
tmp1=plot_genes_branched_heatmap(cds[row.names(subset(states_de,qval<1e-5)),],
                                 branch_point = 1,
                                 num_clusters = 10, 
                                 cores = 1,
                                 branch_labels = c("Cell fate 1", "Cell fate 2"),
                                 #hmcols = NULL, #默认值
                                 hmcols = colorRampPalette(rev(brewer.pal(9, "PRGn")))(62),
                                 branch_colors = c("#979797", "#F05662", "#7990C8"), #pre-branch, Cell fate 1, Cell fate 2
                                 use_gene_short_name = T,
                                 show_rownames = F,
                                 return_heatmap = T #是否返回一些重要信息
)

dev.off()

gene_group=tmp1$annotation_row
gene_group$gene=rownames(gene_group)
write.csv(gene_group,file = "human_coe_outer_by_state_qval_e5.csv")



BEAM_res=BEAM(cds,branch_point = 1,cores = 1)
BEAM_res=BEAM_res[,c("gene_short_name","pval","qval")]
saveRDS(BEAM_res, file = "BEAM_res_coe_inner.rds")

library(RColorBrewer)
a=cds[row.names(subset(BEAM_res,qval<1e-50)),]
#cell fate 2===>inner_HC
tmp1=plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res,qval<1e-4)),],
                                 branch_point = 1,
                                 num_clusters = 10, 
                                 cores = 1,
                                 branch_labels = c("Cell fate 1", "Cell fate 2"),
                                 #hmcols = NULL, #默认值
                                 hmcols = colorRampPalette(rev(brewer.pal(9, "PRGn")))(62),
                                 branch_colors = c("#979797", "#F05662", "#7990C8"), #pre-branch, Cell fate 1, Cell fate 2
                                 use_gene_short_name = T,
                                 show_rownames = F,
                                 return_heatmap = T #是否返回一些重要信息
)
dev.off()

saveRDS(cds,file="human_coe_outer.rds")
save.image("human_coe_outer.RData")

library(ClusterGVis)
#use state to plot
df<-plot_genes_branched_heatmap2(cds[row.names(subset(states_de,qval<1e-2)),],
                                 branch_point = 1,
                                 num_clusters = 10, 
                                 cores = 1,
                                
                                 use_gene_short_name = T,
                                 show_rownames = F,
                                 
                               )
str(df)


visCluster(object=df,plot.type="heatmap")

gene=c("MEIS2","ADAMTSL1","NR2F1","NR2F2","GATA3",
       "SIX1","SOX2","EYA1","BRG1","ATOH1","POU4F3","GFI1","NEUROD1","INSM1","IKZF2",
       "TBX2","NEUROG1","FOXG1","FGF10",
         "LGR5","RORB","RXRA",
       
       "JAG1",  "LFNG","NOTCH1","DLL1","JAG2","HES1","HES5","HEY2","SPRY2",
       "PROX1","FGFR3","BMP2","NGFR","NRCAM","TGFBR1","FZD9",
       "SERPINE2",
       "TMC1","OCM","SLC26A5"
       )
visCluster(object=df,plot.type="heatmap",
           markGenes=gene)
gene_group=df$long.res
gene_group2=gene_group[!duplicated(gene_group$gene),]
write.csv(gene_group2,file = "human_coe_outer_by_state_qval_e2_10clusters_v2.csv")
saveRDS(cds,file="human_coe_outer.rds")
save.image("human_coe_outer.RData")
library(clusterProfiler)
library(org.Hs.eg.db)
allcluster_go=data.frame()
for (i in unique(gene_group$Cluster)) {
  small_gene_group=filter(gene_group,gene_group$Cluster==i)
  df_name=bitr(small_gene_group$gene, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
  go <- enrichGO(gene         = unique(df_name$ENTREZID),
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENTREZID',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.2,
                 readable      = TRUE)
  go_res=go@result
  if (dim(go_res)[1] != 0) {
    go_res$cluster=i
    allcluster_go=rbind(allcluster_go,go_res)
  }
}
head(allcluster_go[,c("ID","Description","qvalue","cluster")])

#####################################################################################
######################################################################################

###################################################################################
####################################################################################
# 批次效应分析函数
batch_effect_analysis <- function(seurat_obj, heatmap_genes, batch_var = "sample") {
  
  # 提取表达矩阵和元数据
  expr_matrix <- GetAssayData(seurat_obj, assay = "SCT", slot = "data")
  metadata <- seurat_obj@meta.data
  
  # 确保批次变量存在
  if (!batch_var %in% colnames(metadata)) {
    stop(paste("Batch variable", batch_var, "not found in metadata"))
  }
  
  # 过滤热图基因中实际存在的基因
  valid_genes <- intersect(heatmap_genes, rownames(expr_matrix))
  cat("Analyzing", length(valid_genes), "genes out of", length(heatmap_genes), "provided genes\n")
  
  # 初始化结果数据框
  batch_effect_results <- data.frame(
    gene = valid_genes,
    p_value = NA,
    adj_p_value = NA,
    significant = FALSE,
    stringsAsFactors = FALSE
  )
  
  # 对每个基因进行批次效应检验
  for (i in seq_along(valid_genes)) {
    gene <- valid_genes[i]
    
    # 提取基因表达量
    gene_expr <- as.numeric(expr_matrix[gene, ])
    
    # 构建数据框
    test_df <- data.frame(
      expression = gene_expr,
      batch = metadata[[batch_var]]
    )
    
    # 使用Kruskal-Wallis检验（非参数ANOVA）
    if (length(unique(test_df$batch)) > 1) {
      kw_test <- kruskal.test(expression ~ batch, data = test_df)
      batch_effect_results$p_value[i] <- kw_test$p.value
    }
    
    # 进度显示
    if (i %% 100 == 0) {
      cat("Processed", i, "genes...\n")
    }
  }
  
  # 多重检验校正
  batch_effect_results$adj_p_value <- p.adjust(batch_effect_results$p_value, method = "BH")
  batch_effect_results$significant <- batch_effect_results$adj_p_value < 0.05
  
  return(batch_effect_results)
}

# 对于inner轨迹分析（对应图2e）
inner_heatmap_genes <- read.csv("human_coe_outer_by_state_qval_e2_10clusters_v2.csv")
inner_batch_results <- batch_effect_analysis(sce_monocle2, inner_heatmap_genes$gene, batch_var = "sample")


# 结果汇总
summarize_batch_effects <- function(results, analysis_name) {
  total_genes <- nrow(results)
  significant_genes <- sum(results$significant, na.rm = TRUE)
  percentage <- round(significant_genes / total_genes * 100, 2)
  
  cat(analysis_name, "分析结果:\n")
  cat("总基因数:", total_genes, "\n")
  cat("受批次显著影响的基因数:", significant_genes, "\n")
  cat("比例:", percentage, "%\n")
  
  # 检查关键基因是否受批次影响
  key_genes <- c("ADAMTSL1", "FOXG1", "HEY2", "NOTCH1", "FGF10", "EYA1",
                 
                 "NRCAM", "GATA3", "TBX2", "RORB", "LFNG", "HES1",
                 "HES5", "LGR5", "FGF20", "NR2F1", "MEIS2", "SOX2",
                 "JAG2", "ATOH1", "SIX1", "GFI1", "SLC17A8", "POU4F3",
                 "OTOF", "TMC1", "DLL1"
                 
  )
  key_gene_results <- results[results$gene %in% key_genes, ]
  
  cat("\n关键基因批次效应分析:\n")
  print(key_gene_results[, c("gene", "adj_p_value", "significant")])
  
  return(list(
    total_genes = total_genes,
    batch_affected = significant_genes,
    percentage = percentage,
    key_gene_results = key_gene_results
  ))
}

# 生成汇总报告
inner_summary <- summarize_batch_effects(inner_batch_results, "Inner轨迹（图2e）")

# 可视化结果
library(ggplot2)
library(patchwork)

# 创建批次效应可视化
plot_batch_summary <- function(inner_summary) {
  summary_df <- data.frame(
    Analysis = c("Inner Trajectory"),
    Total_Genes = c(inner_summary$total_genes),
    Batch_Affected = c(inner_summary$batch_affected),
    Percentage = c(inner_summary$percentage)
  )
  
  p1 <- ggplot(summary_df, aes(x = Analysis, y = Percentage, fill = Analysis)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = paste0(Percentage, "%")), vjust = -0.5) +
    labs(title = "批次效应基因比例",
         y = "受批次影响基因比例 (%)") +
    theme_minimal() +
    scale_fill_brewer(palette = "Set2")
  
  p2 <- ggplot(summary_df, aes(x = Analysis, y = Batch_Affected, fill = Analysis)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = Batch_Affected), vjust = -0.5) +
    labs(title = "受批次影响基因数量",
         y = "基因数量") +
    theme_minimal() +
    scale_fill_brewer(palette = "Set2")
  
  return(p1 + p2)
}

batch_plot <- plot_batch_summary(inner_summary)
print(batch_plot)

# 最终统计报告
cat("\n=== 最终批次效应分析报告 ===\n")
cat("Inner轨迹（图2e）:", inner_summary$batch_affected, "/", inner_summary$total_genes, 
    "(", inner_summary$percentage, "%) 基因受批次显著影响\n")

########################################################################################

# 改进的批次效应分析 - 使用更严格的阈值和方法
refined_batch_effect_analysis <- function(seurat_obj, heatmap_genes, batch_var = "sample") {
  
  # 提取表达矩阵和元数据
  expr_matrix <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")
  metadata <- seurat_obj@meta.data
  
  # 过滤热图基因中实际存在的基因
  valid_genes <- intersect(heatmap_genes, rownames(expr_matrix))
  cat("Analyzing", length(valid_genes), "genes for batch effects\n")
  
  # 初始化结果数据框
  batch_effect_results <- data.frame(
    gene = valid_genes,
    p_value = NA,
    adj_p_value = NA,
    effect_size = NA,  # 添加效应大小
    mean_expr = NA,    # 平均表达量
    significant_strict = FALSE,  # 严格标准
    significant_liberal = FALSE, # 宽松标准
    stringsAsFactors = FALSE
  )
  
  # 对每个基因进行批次效应检验
  for (i in seq_along(valid_genes)) {
    gene <- valid_genes[i]
    
    tryCatch({
      # 提取基因表达量
      gene_expr <- as.numeric(expr_matrix[gene, ])
      
      # 构建数据框
      test_df <- data.frame(
        expression = gene_expr,
        batch = metadata[[batch_var]]
      )
      
      # 移除NA值
      test_df <- test_df[complete.cases(test_df), ]
      
      # 计算平均表达量
      batch_effect_results$mean_expr[i] <- mean(gene_expr, na.rm = TRUE)
      
      # 使用Kruskal-Wallis检验
      if (length(unique(test_df$batch)) > 1 && length(test_df$expression) > 0) {
        kw_test <- kruskal.test(expression ~ batch, data = test_df)
        batch_effect_results$p_value[i] <- kw_test$p.value
        
        # 计算效应大小 (epsilon-squared)
        n <- nrow(test_df)
        h_stat <- kw_test$statistic
        epsilon_squared <- (h_stat - (length(unique(test_df$batch)) - 1)) / (n - 1)
        batch_effect_results$effect_size[i] <- epsilon_squared
      }
      
    }, error = function(e) {
      # 静默处理错误
    })
  }
  
  # 多重检验校正 - 使用更严格的标准
  valid_p_values <- !is.na(batch_effect_results$p_value)
  if (sum(valid_p_values) > 0) {
    batch_effect_results$adj_p_value[valid_p_values] <- p.adjust(
      batch_effect_results$p_value[valid_p_values], method = "BH"
    )
    
    # 严格标准: FDR < 0.01 且 效应大小 > 0.1
    batch_effect_results$significant_strict[valid_p_values] <- 
      batch_effect_results$adj_p_value[valid_p_values] < 0.01 & 
      batch_effect_results$effect_size[valid_p_values] > 0.1
    
    # 宽松标准: FDR < 0.05
    batch_effect_results$significant_liberal[valid_p_values] <- 
      batch_effect_results$adj_p_value[valid_p_values] < 0.05
  }
  
  return(batch_effect_results)
}

# 应用改进的分析方法
cat("=== 应用改进的批次效应分析 ===\n")
inner_batch_refined <- refined_batch_effect_analysis(sce_monocle2, inner_heatmap_genes$gene, batch_var = "sample")

# 改进的结果汇总函数
summarize_refined_batch_effects <- function(results, analysis_name) {
  total_genes <- nrow(results)
  significant_strict <- sum(results$significant_strict, na.rm = TRUE)
  significant_liberal <- sum(results$significant_liberal, na.rm = TRUE)
  
  percentage_strict <- round(significant_strict / total_genes * 100, 2)
  percentage_liberal <- round(significant_liberal / total_genes * 100, 2)
  
  cat(analysis_name, "分析结果:\n")
  cat("总基因数:", total_genes, "\n")
  cat("严格标准 (FDR < 0.01 & 效应大小 > 0.1):", significant_strict, "(", percentage_strict, "%)\n")
  cat("宽松标准 (FDR < 0.05):", significant_liberal, "(", percentage_liberal, "%)\n")
  
  # 检查关键基因
  key_genes <- c("ATOH1", "SOX2", "NOTCH1", "LGR5", "TBX2", "IKZF2", 
                 "FGF20", "HES1", "HES5", "LFNG", "NR2F1")
  key_gene_results <- results[results$gene %in% key_genes, ]
  
  cat("\n关键基因批次效应分析:\n")
  if (nrow(key_gene_results) > 0) {
    print(key_gene_results[, c("gene", "adj_p_value", "effect_size", "significant_strict", "significant_liberal")])
  } else {
    cat("未找到指定的关键基因\n")
  }
  
  return(list(
    total_genes = total_genes,
    batch_affected_strict = significant_strict,
    batch_affected_liberal = significant_liberal,
    percentage_strict = percentage_strict,
    percentage_liberal = percentage_liberal,
    key_gene_results = key_gene_results
  ))
}

# 生成改进的汇总报告
inner_refined_summary <- summarize_refined_batch_effects(inner_batch_refined, "Inner轨迹（图2e）")

# 可视化效应大小分布
library(ggplot2)

# 创建效应大小分布图
plot_effect_size_distribution <- function(results, analysis_name) {
  ggplot(results[!is.na(results$effect_size), ], aes(x = effect_size)) +
    geom_histogram(bins = 50, fill = "lightblue", color = "black") +
    geom_vline(xintercept = 0.1, linetype = "dashed", color = "red") +
    labs(title = paste(analysis_name, "- 批次效应大小分布"),
         x = "效应大小 (epsilon-squared)",
         y = "基因数量") +
    theme_minimal()
}

effect_plot <- plot_effect_size_distribution(inner_batch_refined, "Inner轨迹")
print(effect_plot)

# 保存详细结果
write.csv(inner_batch_refined, "inner_trajectory_batch_effect_refined_analysis.csv", row.names = FALSE)

# 提取真正受批次影响的基因（严格标准）
strong_batch_genes <- inner_batch_refined$gene[inner_batch_refined$significant_strict]
cat("\n强烈受批次影响的基因 (严格标准):", length(strong_batch_genes), "\n")
if (length(strong_batch_genes) > 0) {
  cat("前20个基因:", head(strong_batch_genes, 20), "\n")
}

# 最终报告
cat("\n=== 最终改进的批次效应分析报告 ===\n")
cat("使用严格标准 (FDR < 0.01 & 效应大小 > 0.1):\n")
cat("Inner轨迹:", inner_refined_summary$batch_affected_strict, "/", 
    inner_refined_summary$total_genes, "(", inner_refined_summary$percentage_strict, "%) 基因受批次显著影响\n")

cat("\n使用宽松标准 (FDR < 0.05):\n")
cat("Inner轨迹:", inner_refined_summary$batch_affected_liberal, "/", 
    inner_refined_summary$total_genes, "(", inner_refined_summary$percentage_liberal, "%) 基因受批次显著影响\n")
####################################################################################
# 扩展的多阈值批次效应分析
extended_multi_threshold_batch_analysis <- function(seurat_obj, heatmap_genes, batch_var = "sample") {
  
  # 提取表达矩阵和元数据
  expr_matrix <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")
  metadata <- seurat_obj@meta.data
  
  # 过滤热图基因中实际存在的基因
  valid_genes <- intersect(heatmap_genes, rownames(expr_matrix))
  cat("Analyzing", length(valid_genes), "genes for batch effects\n")
  
  # 初始化结果数据框 - 添加更多效应大小阈值
  batch_effect_results <- data.frame(
    gene = valid_genes,
    p_value = NA,
    adj_p_value = NA,
    effect_size = NA,
    mean_expr = NA,
    
    # 多个效应大小阈值
    significant_effect_0_5 = FALSE,  # FDR < 0.05 & effect > 0.5 (很强)
    significant_effect_0_4 = FALSE,  # FDR < 0.05 & effect > 0.4
    significant_effect_0_3 = FALSE,  # FDR < 0.05 & effect > 0.3
    significant_effect_0_25 = FALSE, # FDR < 0.05 & effect > 0.25
    significant_effect_0_2 = FALSE,  # FDR < 0.05 & effect > 0.2
    significant_effect_0_15 = FALSE, # FDR < 0.05 & effect > 0.15
    significant_effect_0_1 = FALSE,  # FDR < 0.05 & effect > 0.1
    significant_effect_0_05 = FALSE, # FDR < 0.05 & effect > 0.05
    significant_fdr_only = FALSE,    # FDR < 0.05 (无效应大小要求)
    significant_pval_only = FALSE,   # pval < 0.05
    
    stringsAsFactors = FALSE
  )
  
  # 对每个基因进行批次效应检验
  for (i in seq_along(valid_genes)) {
    gene <- valid_genes[i]
    
    tryCatch({
      # 提取基因表达量
      gene_expr <- as.numeric(expr_matrix[gene, ])
      
      # 构建数据框
      test_df <- data.frame(
        expression = gene_expr,
        batch = metadata[[batch_var]]
      )
      
      # 移除NA值
      test_df <- test_df[complete.cases(test_df), ]
      
      # 计算平均表达量
      batch_effect_results$mean_expr[i] <- mean(gene_expr, na.rm = TRUE)
      
      # 使用Kruskal-Wallis检验
      if (length(unique(test_df$batch)) > 1 && length(test_df$expression) > 0) {
        kw_test <- kruskal.test(expression ~ batch, data = test_df)
        batch_effect_results$p_value[i] <- kw_test$p.value
        
        # 计算效应大小 (epsilon-squared)
        n <- nrow(test_df)
        h_stat <- kw_test$statistic
        epsilon_squared <- (h_stat - (length(unique(test_df$batch)) - 1)) / (n - 1)
        batch_effect_results$effect_size[i] <- epsilon_squared
        
        # p-value only threshold
        batch_effect_results$significant_pval_only[i] <- kw_test$p.value < 0.05
      }
      
    }, error = function(e) {
      # 静默处理错误
    })
  }
  
  # 多重检验校正
  valid_p_values <- !is.na(batch_effect_results$p_value)
  if (sum(valid_p_values) > 0) {
    batch_effect_results$adj_p_value[valid_p_values] <- p.adjust(
      batch_effect_results$p_value[valid_p_values], method = "BH"
    )
    
    # 多个效应大小阈值标准
    batch_effect_results$significant_effect_0_5[valid_p_values] <- 
      batch_effect_results$adj_p_value[valid_p_values] < 0.05 & 
      batch_effect_results$effect_size[valid_p_values] > 0.5
    
    batch_effect_results$significant_effect_0_4[valid_p_values] <- 
      batch_effect_results$adj_p_value[valid_p_values] < 0.05 & 
      batch_effect_results$effect_size[valid_p_values] > 0.4
    
    batch_effect_results$significant_effect_0_3[valid_p_values] <- 
      batch_effect_results$adj_p_value[valid_p_values] < 0.05 & 
      batch_effect_results$effect_size[valid_p_values] > 0.3
    
    batch_effect_results$significant_effect_0_25[valid_p_values] <- 
      batch_effect_results$adj_p_value[valid_p_values] < 0.05 & 
      batch_effect_results$effect_size[valid_p_values] > 0.25
    
    batch_effect_results$significant_effect_0_2[valid_p_values] <- 
      batch_effect_results$adj_p_value[valid_p_values] < 0.05 & 
      batch_effect_results$effect_size[valid_p_values] > 0.2
    
    batch_effect_results$significant_effect_0_15[valid_p_values] <- 
      batch_effect_results$adj_p_value[valid_p_values] < 0.05 & 
      batch_effect_results$effect_size[valid_p_values] > 0.15
    
    batch_effect_results$significant_effect_0_1[valid_p_values] <- 
      batch_effect_results$adj_p_value[valid_p_values] < 0.05 & 
      batch_effect_results$effect_size[valid_p_values] > 0.1
    
    batch_effect_results$significant_effect_0_05[valid_p_values] <- 
      batch_effect_results$adj_p_value[valid_p_values] < 0.05 & 
      batch_effect_results$effect_size[valid_p_values] > 0.05
    
    batch_effect_results$significant_fdr_only[valid_p_values] <- 
      batch_effect_results$adj_p_value[valid_p_values] < 0.05
  }
  
  return(batch_effect_results)
}

# 应用扩展的多阈值分析
cat("=== 扩展的多阈值批次效应分析 ===\n")
inner_batch_extended <- extended_multi_threshold_batch_analysis(sce_monocle2, inner_heatmap_genes$gene, batch_var = "sample")

# 扩展的多阈值结果汇总函数
summarize_extended_batch_effects <- function(results, analysis_name) {
  total_genes <- nrow(results)
  
  # 计算各个阈值的显著基因数
  counts <- list(
    effect_0_5 = sum(results$significant_effect_0_5, na.rm = TRUE),
    effect_0_4 = sum(results$significant_effect_0_4, na.rm = TRUE),
    effect_0_3 = sum(results$significant_effect_0_3, na.rm = TRUE),
    effect_0_25 = sum(results$significant_effect_0_25, na.rm = TRUE),
    effect_0_2 = sum(results$significant_effect_0_2, na.rm = TRUE),
    effect_0_15 = sum(results$significant_effect_0_15, na.rm = TRUE),
    effect_0_1 = sum(results$significant_effect_0_1, na.rm = TRUE),
    effect_0_05 = sum(results$significant_effect_0_05, na.rm = TRUE),
    fdr_only = sum(results$significant_fdr_only, na.rm = TRUE),
    pval_only = sum(results$significant_pval_only, na.rm = TRUE)
  )
  
  percentages <- lapply(counts, function(x) round(x / total_genes * 100, 2))
  
  cat(analysis_name, "扩展多阈值分析结果:\n")
  cat("总基因数:", total_genes, "\n\n")
  
  cat("效应大小阈值分析 (FDR < 0.05):\n")
  cat("1. 效应大小 > 0.5 (极强):", counts$effect_0_5, "(", percentages$effect_0_5, "%)\n")
  cat("2. 效应大小 > 0.4 (很强):", counts$effect_0_4, "(", percentages$effect_0_4, "%)\n")
  cat("3. 效应大小 > 0.3 (强):", counts$effect_0_3, "(", percentages$effect_0_3, "%)\n")
  cat("4. 效应大小 > 0.25 (中强):", counts$effect_0_25, "(", percentages$effect_0_25, "%)\n")
  cat("5. 效应大小 > 0.2 (中等):", counts$effect_0_2, "(", percentages$effect_0_2, "%)\n")
  cat("6. 效应大小 > 0.15 (中弱):", counts$effect_0_15, "(", percentages$effect_0_15, "%)\n")
  cat("7. 效应大小 > 0.1 (弱):", counts$effect_0_1, "(", percentages$effect_0_1, "%)\n")
  cat("8. 效应大小 > 0.05 (很弱):", counts$effect_0_05, "(", percentages$effect_0_05, "%)\n\n")
  
  cat("无效应大小要求:\n")
  cat("9. 仅FDR < 0.05:", counts$fdr_only, "(", percentages$fdr_only, "%)\n")
  cat("10. 仅pval < 0.05:", counts$pval_only, "(", percentages$pval_only, "%)\n")
  
  # 检查关键基因
  key_genes <- c("ATOH1", "SOX2", "NOTCH1", "LGR5", "TBX2", "IKZF2", 
                 "FGF20", "HES1", "HES5", "LFNG", "NR2F1")
  key_gene_results <- results[results$gene %in% key_genes, ]
  
  cat("\n关键基因在不同效应大小阈值下的批次效应:\n")
  if (nrow(key_gene_results) > 0) {
    key_summary <- key_gene_results[, c("gene", "adj_p_value", "effect_size")]
    
    # 添加各个阈值的结果
    thresholds <- c("0.5", "0.4", "0.3", "0.25", "0.2", "0.15", "0.1", "0.05")
    for(thresh in thresholds) {
      col_name <- paste0("significant_effect_", gsub("\\.", "_", thresh))
      key_summary[[paste0("effect_>", thresh)]] <- key_gene_results[[col_name]]
    }
    
    print(key_summary)
  } else {
    cat("未找到指定的关键基因\n")
  }
  
  return(list(
    total_genes = total_genes,
    counts = counts,
    percentages = percentages,
    key_gene_results = key_gene_results
  ))
}

# 生成扩展的汇总报告
inner_extended_summary <- summarize_extended_batch_effects(inner_batch_extended, "Inner轨迹（图2e）")

# 创建效应大小阈值比较可视化
plot_effect_threshold_comparison <- function(results, analysis_name) {
  threshold_data <- data.frame(
    Effect_Threshold = c(">0.5", ">0.4", ">0.3", ">0.25", ">0.2", ">0.15", ">0.1", ">0.05", "FDR only", "pval only"),
    Genes_Affected = c(
      sum(results$significant_effect_0_5, na.rm = TRUE),
      sum(results$significant_effect_0_4, na.rm = TRUE),
      sum(results$significant_effect_0_3, na.rm = TRUE),
      sum(results$significant_effect_0_25, na.rm = TRUE),
      sum(results$significant_effect_0_2, na.rm = TRUE),
      sum(results$significant_effect_0_15, na.rm = TRUE),
      sum(results$significant_effect_0_1, na.rm = TRUE),
      sum(results$significant_effect_0_05, na.rm = TRUE),
      sum(results$significant_fdr_only, na.rm = TRUE),
      sum(results$significant_pval_only, na.rm = TRUE)
    ),
    Percentage = c(
      round(sum(results$significant_effect_0_5, na.rm = TRUE)/nrow(results)*100, 1),
      round(sum(results$significant_effect_0_4, na.rm = TRUE)/nrow(results)*100, 1),
      round(sum(results$significant_effect_0_3, na.rm = TRUE)/nrow(results)*100, 1),
      round(sum(results$significant_effect_0_25, na.rm = TRUE)/nrow(results)*100, 1),
      round(sum(results$significant_effect_0_2, na.rm = TRUE)/nrow(results)*100, 1),
      round(sum(results$significant_effect_0_15, na.rm = TRUE)/nrow(results)*100, 1),
      round(sum(results$significant_effect_0_1, na.rm = TRUE)/nrow(results)*100, 1),
      round(sum(results$significant_effect_0_05, na.rm = TRUE)/nrow(results)*100, 1),
      round(sum(results$significant_fdr_only, na.rm = TRUE)/nrow(results)*100, 1),
      round(sum(results$significant_pval_only, na.rm = TRUE)/nrow(results)*100, 1)
    )
  )
  
  threshold_data$Effect_Threshold <- factor(threshold_data$Effect_Threshold, 
                                            levels = threshold_data$Effect_Threshold)
  
  ggplot(threshold_data, aes(x = Effect_Threshold, y = Percentage, fill = Effect_Threshold)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = paste0(Percentage, "%")), vjust = -0.5, size = 2.5) +
    labs(title = paste(analysis_name, "- 不同效应大小阈值下的批次效应"),
         y = "受批次影响基因比例 (%)",
         x = "效应大小阈值 (FDR < 0.05)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    scale_fill_brewer(palette = "Spectral")
}

effect_threshold_plot <- plot_effect_threshold_comparison(inner_batch_extended, "Inner轨迹")
print(effect_threshold_plot)

# 保存扩展结果
write.csv(inner_batch_extended, "outer_trajectory_batch_effect_extended_thresholds.csv", row.names = FALSE)

# 最终详细报告
cat("\n" + rep("=", 80) + "\n")
cat("                        扩展多阈值批次效应分析最终报告\n")
cat(rep("=", 80) + "\n\n")

cat("基于您要求的特定阈值:\n")
cat("• FDR < 0.05 & 效应大小 > 0.25:", inner_extended_summary$counts$effect_0_25, "/", 
    inner_extended_summary$total_genes, "(", inner_extended_summary$percentages$effect_0_25, "%)\n")
cat("• FDR < 0.05 & 效应大小 > 0.2:", inner_extended_summary$counts$effect_0_2, "/", 
    inner_extended_summary$total_genes, "(", inner_extended_summary$percentages$effect_0_2, "%)\n")
cat("• FDR < 0.05 & 效应大小 > 0.5:", inner_extended_summary$counts$effect_0_5, "/", 
    inner_extended_summary$total_genes, "(", inner_extended_summary$percentages$effect_0_5, "%)\n\n")

cat("效应大小解释 (epsilon-squared):\n")
cat("• > 0.5: 极强效应 - 批次解释了 >50% 的表达变异\n")
cat("• 0.3-0.5: 强效应 - 批次解释了 30-50% 的表达变异\n")
cat("• 0.1-0.3: 中等效应 - 批次解释了 10-30% 的表达变异\n")
cat("• < 0.1: 弱效应 - 批次解释了 <10% 的表达变异\n\n")

# 检查核心基因在您关注的阈值下的状态
cat("核心基因在关键阈值下的状态:\n")
thresholds_to_check <- c("0.25", "0.2", "0.5")
if (nrow(inner_extended_summary$key_gene_results) > 0) {
  for (thresh in thresholds_to_check) {
    cat("\n在 FDR<0.05 & 效应大小 >", thresh, "标准下:\n")
    col_name <- paste0("significant_effect_", gsub("\\.", "_", thresh))
    
    affected_genes <- inner_extended_summary$key_gene_results$gene[inner_extended_summary$key_gene_results[[col_name]]]
    unaffected_genes <- inner_extended_summary$key_gene_results$gene[!inner_extended_summary$key_gene_results[[col_name]]]
    
    if (length(affected_genes) > 0) {
      cat("  受影响: ", paste(affected_genes, collapse = ", "), "\n")
    }
    if (length(unaffected_genes) > 0) {
      cat("  不受影响: ", paste(unaffected_genes, collapse = ", "), "\n")
    }
  }
}

# 创建审稿人回复用的简洁表格
create_final_reviewer_table <- function(extended_summary) {
  final_table <- data.frame(
    Threshold = c(
      "FDR<0.05 & Effect>0.5",
      "FDR<0.05 & Effect>0.4", 
      "FDR<0.05 & Effect>0.3",
      "FDR<0.05 & Effect>0.25",
      "FDR<0.05 & Effect>0.2",
      "FDR<0.05 & Effect>0.1",
      "FDR<0.05 only"
    ),
    Genes_Affected = c(
      paste0(extended_summary$counts$effect_0_5, "/", extended_summary$total_genes),
      paste0(extended_summary$counts$effect_0_4, "/", extended_summary$total_genes),
      paste0(extended_summary$counts$effect_0_3, "/", extended_summary$total_genes),
      paste0(extended_summary$counts$effect_0_25, "/", extended_summary$total_genes),
      paste0(extended_summary$counts$effect_0_2, "/", extended_summary$total_genes),
      paste0(extended_summary$counts$effect_0_1, "/", extended_summary$total_genes),
      paste0(extended_summary$counts$fdr_only, "/", extended_summary$total_genes)
    ),
    Percentage = c(
      paste0(extended_summary$percentages$effect_0_5, "%"),
      paste0(extended_summary$percentages$effect_0_4, "%"),
      paste0(extended_summary$percentages$effect_0_3, "%"),
      paste0(extended_summary$percentages$effect_0_25, "%"),
      paste0(extended_summary$percentages$effect_0_2, "%"),
      paste0(extended_summary$percentages$effect_0_1, "%"),
      paste0(extended_summary$percentages$fdr_only, "%")
    ),
    Effect_Strength = c("Very Strong", "Strong", "Moderate-Strong", "Moderate", 
                        "Moderate-Weak", "Weak", "Statistical only")
  )
  
  return(final_table)
}

reviewer_final_table <- create_final_reviewer_table(inner_extended_summary)
print(reviewer_final_table)
write.csv(reviewer_final_table, "batch_effect_final_summary_for_reviewers_outer.csv", row.names = FALSE)

######################################################################################


