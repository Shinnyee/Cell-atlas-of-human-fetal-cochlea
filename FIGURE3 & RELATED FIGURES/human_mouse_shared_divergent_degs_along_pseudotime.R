library(monocle)
library(Seurat)
library(ggplot2)
library(sctransform)
library(tidyverse)
options(stringsAsFactors = FALSE)
library(dittoSeq)
library(Seurat)
library(dplyr)
library(sctransform)

load("C:/Users/Dell/Desktop/WORKPLACE/R/HC TRAJOTORY/human_hc_monocle2.Rdata")
cds_human=my_cds_subset
colnames(pData(cds_human))
table(pData(cds_human)$State,pData(cds_human)$Cluster)
library(ggsci)
plot_cell_trajectory(cds_human, markers_linear=F, show_branch_points=F, color_by = "Pseudotime")  + 
  scale_color_gradient(low="grey", high="sienna2")

plot_cell_trajectory(cds_human, markers = c("SOX2", "ATOH1","OTOF","TBX2",
                                      "SLC17A8","SLC26A5","IKZF2","FGFR3",
                                      "TMC1","MAGI3","ANKFN1","DACH1"), 
                     markers_linear=F, use_color_gradient=T, 
                     show_branch_points=F,cell_size = 3)

human_hc_GW_ALL <- readRDS("F:/Chai Lab/PROJECT_INTEGRATION ACROSS SPECIES/Cochlea_Evo_new/HC TRAJOTORY/human_hc_GW_ALL.rds")
sce=human_hc_GW_ALL
markers.to.plot<-c(
  "SOX2", "ATOH1","OTOF","TBX2",
  "SLC17A8","SLC26A5","IKZF2","FGFR3",
  "TMC1","MAGI3","ANKFN1","DACH1"
  
)
Idents(sce) <- "celltype"
VlnPlot(sce,features = markers.to.plot)
p6=plot_cell_trajectory(cds_human, color_by = "celltype",cell_size = 2) +
  scale_color_manual(breaks = c("iOHC", "iIHC","OHC", "IHC"), 
                     values=c("#c8a485", "#9ac7ac", "#f7b741","#55efc8"))
p7=plot_cell_trajectory(cds_human, color_by = "Pseudotime",cell_size = 2)  
p6+p7
############################################################################################
library(dplyr)
load("C:/Users/Dell/Desktop/WORKPLACE/R/HC TRAJOTORY/mouse data/mouse_hc_monocle2.Rdata")
cds_mouse=my_cds_subset
colnames(pData(cds_mouse))
table(pData(cds_mouse)$State,pData(cds_mouse)$Cluster)
library(ggsci)
plot_cell_trajectory(cds_mouse, markers = c("Sox2", "Atoh1","Otof","Tbx2",
                                      "Slc17a8","Slc26a5","Ikzf2","Fgfr3",
                                      "Tmc1","Magi3","Ankfn1","Dach1"), 
                     markers_linear=F, use_color_gradient=T, 
                     show_branch_points=F,cell_size = 1.5)#"Magi3","Ankfn1"
markers.to.plot<-c(
  "Sox2", "Atoh1","Otof","Tbx2",
  "Slc17a8","Slc26a5","Ikzf2","Fgfr3",
  "Tmc1","Magi3","Ankfn1","Dach1"
  
)
mouse_hc_days_ALL <- readRDS("C:/Users/Dell/Desktop/WORKPLACE/R/HC TRAJOTORY/mouse data/mouse_hc_days_ALL.rds")
mouse_hc_days_ALL <- readRDS("F:/Chai Lab/PROJECT_INTEGRATION ACROSS SPECIES/Cochlea_Evo_new/HC TRAJOTORY/mouse data/mouse_hc_days_ALL.rds")
sce_mouse=mouse_hc_days_ALL
Idents(sce_mouse) <- "celltype"
VlnPlot(sce_mouse,features = markers.to.plot,assay = "SCT")
FeaturePlot(sce,features = c("IKZF2","ATOH1","TBX2","DACH1"),label = TRUE,pt.size = 2)
FeaturePlot(sce_mouse,features = c("Ikzf2","Atoh1","Tbx2","Dach1"),label = TRUE,pt.size = 2)

# for human
pData(cds_human)$Cluster=pData(cds_human)$celltype
table(pData(cds_human)$Cluster)

Sys.time()
diff_test_res_human <- differentialGeneTest(cds_human,
                                            fullModelFormulaStr = "~Cluster")
Sys.time()
# Select genes that are significant
sig_genes_human <- subset(diff_test_res_human, pval < 1e-3)
sig_genes_human=sig_genes_human[order(sig_genes_human$pval),]
head(sig_genes_human[,c("gene_short_name", "pval", "qval")] ) 


# for mouse
pData(cds_mouse)$Cluster=pData(cds_mouse)$celltype
table(pData(cds_mouse)$Cluster)

Sys.time()
diff_test_res_mouse <- differentialGeneTest(cds_mouse,
                                            fullModelFormulaStr = "~Cluster")
Sys.time()

# Select genes that are significant
sig_genes_mouse <- subset(diff_test_res_mouse, pval < 1e-35)
sig_genes_mouse=sig_genes_mouse[order(sig_genes_mouse$pval),]
head(sig_genes_mouse[,c("gene_short_name", "pval", "qval")] )

# mouse gene id converted into human gene id
genesV2 <- read.table("mouse_to_human_genes.txt", sep="\t", header=T)
sig_genes_mouse$gene_short_name <- genesV2[match(sig_genes_mouse$gene_short_name, genesV2[,1]),2]
sig_genes_mouse <- subset(sig_genes_mouse, gene_short_name!='NA')
sig_genes_mouse <- dplyr::select(sig_genes_mouse, gene_short_name, everything())
sig_genes_mouse <- sig_genes_mouse[, !(colnames(sig_genes_mouse) %in% 'gene')]
dim(sig_genes_mouse)
# extract conserved genes across species
all.genes <- data.frame(genes = unique(c(sig_genes_human$gene_short_name, 
                                         
                                         sig_genes_mouse$gene_short_name)))

all.genes$Human <- as.character(match(all.genes$genes, sig_genes_human$gene_short_name))

all.genes$Mouse <- as.character(match(all.genes$genes, sig_genes_mouse$gene_short_name))
all.genes$Human[which(is.na(all.genes$Human))] <- FALSE
all.genes$Human[which(all.genes$Human != FALSE)] <- TRUE

all.genes$Mouse[which(is.na(all.genes$Mouse))] <- FALSE
all.genes$Mouse[which(all.genes$Mouse != FALSE)] <- TRUE
all.genes$Human <- as.logical(all.genes$Human)

all.genes$Mouse <- as.logical(all.genes$Mouse)
library(eulerr)
plot(euler(
  all.genes[ ,2:3]),
  quantities = list(cex = 3),
  labels = TRUE,
  
  fills = c("royalblue1", "maroon4", "sienna2")
)
# conserved genes
all.genes <- all.genes[which(all.genes[,2] == TRUE), ]
all.genes <- all.genes[which(all.genes[,3] == TRUE), ]
write.csv(all.genes,file = "shared_degs_monocle2.csv")
shared_genes <- read.csv("shared_degs_monocle2.csv",row.names = 1)
gc()
################################################################################################
library(ggsci)
# put 1125 conserved genes to plot in a pseudotime form.
Time_genes <- all.genes %>% pull(genes) %>% as.character()

plot_pseudotime_heatmap(cds_human[Time_genes,], num_clusters=3, show_rownames=T, return_heatmap=T)
plot_genes_branched_heatmap(cds_human[Time_genes,],
                            branch_point = 1,
                            num_clusters = 4,
                            cores = 4,
                            use_gene_short_name = T,
                            show_rownames = T,branch_colors = c("#979797", "#F05662", "#7990C8"))

# human gene id converted into mouse gene id
genesV2 <- read.table("mouse_to_human_genes.txt", sep="\t", header=T)
all.genes2=all.genes
all.genes2$genes <- genesV2[match(all.genes2$genes, genesV2[,2]),1]
all.genes2 <- subset(all.genes2, genes!='NA')
all.genes2 <- dplyr::select(all.genes2, genes, everything())
all.genes2 <- all.genes2[, !(colnames(all.genes2) %in% 'gene')]
dim(all.genes2)
Time_genes2 <- all.genes2 %>% pull(genes) %>% as.character()

plot_pseudotime_heatmap(cds_mouse[Time_genes2,], num_clusters=3, show_rownames=T, return_heatmap=T)
plot_genes_branched_heatmap(cds_mouse[Time_genes2,],
                            branch_point = 1,
                            num_clusters = 4,
                            cores = 4,
                            use_gene_short_name = T,
                            show_rownames = T,branch_colors = c("#979797", "#F05662", "#7990C8"))
# extract human-specific genes across species
all.genes <- data.frame(genes = unique(c(sig_genes_human$gene_short_name, 
                                         
                                         sig_genes_mouse$gene_short_name)))

all.genes$Human <- as.character(match(all.genes$genes, sig_genes_human$gene_short_name))

all.genes$Mouse <- as.character(match(all.genes$genes, sig_genes_mouse$gene_short_name))
all.genes$Human[which(is.na(all.genes$Human))] <- FALSE
all.genes$Human[which(all.genes$Human != FALSE)] <- TRUE

all.genes$Mouse[which(is.na(all.genes$Mouse))] <- FALSE
all.genes$Mouse[which(all.genes$Mouse != FALSE)] <- TRUE
all.genes$Human <- as.logical(all.genes$Human)

all.genes$Mouse <- as.logical(all.genes$Mouse)
library(eulerr)
plot(euler(
  all.genes[ ,2:3]),
  quantities = list(cex = 3),
  labels = TRUE,
  
  fills = c("royalblue1", "maroon4", "sienna2")
)
# human-specific genes
all.genes <- all.genes[which(all.genes[,2] == TRUE), ]
all.genes <- all.genes[which(all.genes[,3] == FALSE), ]
write.csv(all.genes,file = "human_specific_degs_monocle2.csv")
human_genes <- read.csv("human_specific_degs_monocle2.csv",row.names = 1)
gc()
library(ggsci)
# put 2181 conserved genes to plot in a pseudotime form.
Time_genes <- all.genes %>% pull(genes) %>% as.character()

plot_pseudotime_heatmap(cds_human[Time_genes,], num_clusters=3, show_rownames=T, return_heatmap=T)
plot_genes_branched_heatmap(cds_human[Time_genes,],
                            branch_point = 1,
                            num_clusters = 4,
                            cores = 4,
                            use_gene_short_name = T,
                            show_rownames = T,branch_colors = c("#979797", "#F05662", "#7990C8"))


# extract mouse-specific genes across species
all.genes <- data.frame(genes = unique(c(sig_genes_human$gene_short_name, 
                                         
                                         sig_genes_mouse$gene_short_name)))

all.genes$Human <- as.character(match(all.genes$genes, sig_genes_human$gene_short_name))

all.genes$Mouse <- as.character(match(all.genes$genes, sig_genes_mouse$gene_short_name))
all.genes$Human[which(is.na(all.genes$Human))] <- FALSE
all.genes$Human[which(all.genes$Human != FALSE)] <- TRUE

all.genes$Mouse[which(is.na(all.genes$Mouse))] <- FALSE
all.genes$Mouse[which(all.genes$Mouse != FALSE)] <- TRUE
all.genes$Human <- as.logical(all.genes$Human)

all.genes$Mouse <- as.logical(all.genes$Mouse)
library(eulerr)
plot(euler(
  all.genes[ ,2:3]),
  quantities = list(cex = 3),
  labels = TRUE,
  
  fills = c("royalblue1", "maroon4", "sienna2")
)
# mouse-specific genes
all.genes <- all.genes[which(all.genes[,2] == FALSE), ]
all.genes <- all.genes[which(all.genes[,3] == TRUE), ]
write.csv(all.genes,file = "mouse_specific_degs_monocle2.csv")
mouse_genes <- read.csv("mouse_specific_degs_monocle2.csv",row.names = 1)
gc()
library(ggsci)


# human gene id converted into mouse gene id
genesV2 <- read.table("mouse_to_human_genes.txt", sep="\t", header=T)
all.genes2=all.genes
all.genes2$genes <- genesV2[match(all.genes2$genes, genesV2[,2]),1]
all.genes2 <- subset(all.genes2, genes!='NA')
all.genes2 <- dplyr::select(all.genes2, genes, everything())
all.genes2 <- all.genes2[, !(colnames(all.genes2) %in% 'gene')]
dim(all.genes2)
Time_genes2 <- all.genes2 %>% pull(genes) %>% as.character()
dev.off()
plot_pseudotime_heatmap(cds_mouse[Time_genes2,], num_clusters=3, show_rownames=T, return_heatmap=T)
plot_genes_branched_heatmap(cds_mouse[Time_genes2,],
                            branch_point = 1,
                            num_clusters = 4,
                            cores = 4,
                            use_gene_short_name = T,
                            show_rownames = T,branch_colors = c("#979797", "#F05662", "#7990C8"))
human_hc_GW_ALL <- readRDS("F:/Chai Lab/PROJECT_INTEGRATION ACROSS SPECIES/Cochlea_Evo_new/HC TRAJOTORY/human_hc_GW_ALL.rds")
mouse_hc_days_ALL <- readRDS("F:/Chai Lab/PROJECT_INTEGRATION ACROSS SPECIES/Cochlea_Evo_new/HC TRAJOTORY/mouse data/mouse_hc_days_ALL.rds")
sce_mouse=mouse_hc_days_ALL
sce_human=human_hc_GW_ALL

Idents(sce_mouse) <- "celltype"
Idents(sce_human) <- "celltype"
VlnPlot(sce_human,features = c("POU4F1","STAB2","KCNB1"),assay = "SCT")
VlnPlot(sce_mouse,features = c("Ankrd22","Slc8a2","Kcnq4"),assay = "SCT")

##########################
#######see known early or late expressed genes in pseudo-plot
genes <- c("SOX2", "ATOH1","OTOF","TBX2",
           "SLC17A8","SLC26A5","IKZF2","FGFR3",
           "TMC1")
genes <- data.frame(genes)

a <- merge(x=shared_genes,y=genes, by="genes", all=FALSE)
### look for tfs in shared list
tf_all <- read.csv("tf_gene_list.csv")
shared_tf <- merge(x=shared_genes,y=tf_all, by="genes", all=FALSE)
write.csv(shared_tf,file = "shared_tf_human_mouse.csv")
#for human specific
tf_all <- read.csv("tf_gene_list.csv")
human_tf <- merge(x=human_genes,y=tf_all, by="genes", all=FALSE)
write.csv(human_tf,file = "human_specific_tf.csv")
#for mouse specific
tf_all <- read.csv("tf_gene_list.csv")
mouse_tf <- merge(x=mouse_genes,y=tf_all, by="genes", all=FALSE)
write.csv(mouse_tf,file = "mouse_specific_tf.csv")
human_enriched_genes <- read.csv("hc_human_enriched_genes.csv")
mouse_enriched_genes <- read.csv("hc_mouse_enriched_genes.csv")

human_tf_hc_enriched <- merge(x=human_enriched_genes,y=human_tf, by="genes", all=FALSE)
mouse_tf_hc_enriched <- merge(x=mouse_enriched_genes,y=mouse_tf, by="genes", all=FALSE)
write.csv(human_tf_hc_enriched,file = "human_tf_hc_enriched.csv")
write.csv(mouse_tf_hc_enriched,file = "mouse_tf_hc_enriched.csv")
