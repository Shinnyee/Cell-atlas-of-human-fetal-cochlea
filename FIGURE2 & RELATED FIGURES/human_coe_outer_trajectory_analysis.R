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