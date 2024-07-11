library(Seurat)
library(pathfindR)
library(CellChat)
library(NMF)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(cowplot)
library(dittoSeq)
library(RColorBrewer)
library(Nebulosa)
library(scCustomize)
library(qs)
library(scuttle)
library(ggpubr)
library(reshape2)
library(ggalluvial)




###### colors ########

#set color palletes for plotting 
colorpalette <- c("#f4a582","#d6604d","#7fc97f","#33a02c","#a6cee3","#1f78b4","#984ea3","#ff7f00")
condition_colors <- c('#919cb3', '#233867','#f5af91','#eb5f22')

####### data ########

#read the data
seu <- readRDS("data/seu.RDS")
load('data/bulkRNAseq.RData')

####### young #########

Idents(seu) <- "condition"

#choose cells from two conditions only
young_idents<-WhichCells(object = seu, idents = c("control young", "repopulated young"))
young <- seu[,match(young_idents,colnames(seu))]
young$condition <- droplevels(young$condition)
new_condition <- c("control", "repopulated")
names(new_condition) <- levels(young)

#rename clusters
young <- RenameIdents(young, new_condition)
young$condition <- Idents(young)

##### Figure 2 ########

#dimplot with control and repopulated 
Idents(young) <- "celltypes"
DimPlot(young, split.by = "condition",  cols = alpha(colorpalette, 0.6), label = TRUE, label.size = 4, ncol = 1)

#celltype contribution
dittoBarPlot(young,"celltypes",group.by = "condition", color.panel = colorpalette,var.labels.reorder = c(4,1,5,2,6,3,7,8),
             main = "", xlab = "Condition", ylab = "Percentage of cells") + theme(axis.title.x = element_text(size = 15),
                                                                                  axis.text.x = element_text(size = 12, colour = "black"),
                                                                                  axis.title.y = element_text(size = 15),
                                                                                  axis.text.y = element_text(size = 12),
                                                                                  legend.text = element_text(size = 12))

#plot all chosen genes
#chosen markers from all microglia subtypes
markers <- as.factor(c("Tmem119", "P2ry12", "Jun", "Fos", "Bmp2k", "Notch2", 
                       'H2-D1','H2-Oa','H2-DMa',
                       "Esco1", "Cdk1", "H2afx",
                       "Csf1", "Mif", "Ifit1", "Irf7"))


DotPlot(young, features =markers) + RotatedAxis()


#Add module scores
proliferating_markers <- c("Cdkn2a", "Esco1", "Esco2",  "Cdk1", "H2afx")
inflammation_markers <- c("Mcm5","Ifitm3", "Irf7", "Isg15", "Usp18", "Mif")
premg_markers <- c("Csf1", "Mcm5", "Cst7", "Mif", "Ccl12", "Ccl3", "Ccl4", 
                   "Ifit1", "Ifit3", "Ifit3b", "Ifitm3", "Irf7", "Isg15", "Usp18")
phagocytosis <- c('Gpnmb', 'Lgals3','Fabp5', 'Apoe','Spp1', 'Cstb','Cstdc1')
cytokine <- c('Tnf', 'Ccl2','Ccl3','Ccl4','Ccl7','Ccl12')

young <- AddModuleScore(
  object = young,
  features = list(phagocytosis),
  name = 'Phagocytosis'
)

premg <- RidgePlot(young, features = 'Premg1', cols = colorpalette)  + 
  labs(title = "Premature microglia", subtitle = 'Csf1, Mcm5, Cst7, Mif, Ccl12, Ccl3, \
       Ccl4, Ifit1, Ifit3, Ifit3b, Ifitm3, Irf7, \
       Isg15, Usp18') + 
  theme(plot.subtitle = element_text(face = "italic"), legend.position = 'none',
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),)

proli <- RidgePlot(young, features = 'Prolifmg1', cols = colorpalette)  + 
  labs(title = "Proliferation", subtitle = 'Ki67, Cdkn2a, Esco1, \
       Esco2, Cdk1, H2afx') + 
  theme(plot.subtitle = element_text(face = "italic"), legend.position = 'none',
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

infl <-  RidgePlot(young, features = 'Inflammation1', cols = colorpalette)  + 
  labs(title = "Inflammation", subtitle = 'Mcm5, Ifitm3, Irf7, Isg15, \
       Usp18, Mif') + 
  theme(plot.subtitle = element_text(face = "italic"), legend.position = 'none',
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

cyt <- RidgePlot(young, features = 'Cytokines1', cols = colorpalette)  + 
  labs(title = "Cytokines", subtitle = 'Tnf, Ccl2, Ccl3, Ccl4, \
       Ccl7, Ccl12') + 
  theme(plot.subtitle = element_text(face = "italic"), legend.position = 'none',
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

phagocyt <- RidgePlot(young, features = 'Phagocytosis1', cols = colorpalette)  + 
  labs(title = "Phagocytosis", subtitle = 'Gpnmb, Lgals3, Fabp5, Apoe, \
       Spp1, Cstb, Cstdc1') + 
  theme(plot.subtitle = element_text(face = "italic"), legend.position = 'none',
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

plot_grid(premg, proli, infl, cyt, phagocyt, ncol = 5)


#Feature density plot of Mki67
Plot_Density_Custom(young, features = "Mki67", custom_palette = brewer.pal(5, "Purples")) + ggtitle("Ki67") + 
  theme(plot.title = element_text(size=20, face="italic", hjust = 0.5))

#Reactome pathway of pre-MG
#search for the genes which identifies activated cell types
celltypes <- FindAllMarkers(young, logfc.threshold = 0.25,  min.pct = 0.25, test.use = "negbinom")

preMG <- celltypes[which(celltypes$cluster == "pre-MG"),]

preMG_df <- as.data.frame(rownames(preMG))
preMG_df <- cbind(preMG_df, preMG$avg_log2FC, preMG$p_val_adj)
names(preMG_df) <- c('Gene.symbol','logFC','adj.P.Val')
preMG_df$Gene.symbol <- toupper(preMG_df$Gene.symbol)
preMG_Reactome <- run_pathfindR(preMG_df, gene_sets = 'Reactome')

preMG_r <- preMG_Reactome

preMG_10 <- preMG_r[which(preMG_r$ID == 'R-HSA-1169410'),]
preMG_10 <- rbind(preMG_10, preMG_r[which(preMG_r$ID == 'R-HSA-909733'),],preMG_r[which(preMG_r$ID == 'R-HSA-975138'),], 
                  preMG_r[which(preMG_r$ID == 'R-HSA-9020702'),],preMG_r[which(preMG_r$ID == 'R-HSA-9758274'),], 
                  preMG_r[which(preMG_r$ID == 'R-HSA-1266695'),],preMG_r[which(preMG_r$ID == 'R-HSA-912694'),],
                  preMG_r[which(preMG_r$ID == 'R-HSA-5218859'),],preMG_r[which(preMG_r$ID == 'R-HSA-445989'),],
                  preMG_r[which(preMG_r$ID == 'R-HSA-877300'),],preMG_r[which(preMG_r$ID == 'R-HSA-5357801'),],
                  preMG_r[which(preMG_r$ID == 'R-HSA-6785807'),],preMG_r[which(preMG_r$ID == 'R-HSA-877312'),])


enrichment_chart(
  result_df = preMG_10,
  top_terms = 13,
  num_bubbles = 5
) +scale_color_gradientn(colors=viridis::viridis(n=15)) + theme(axis.text.y =element_text(size=12))


##### Figure 3 ######
#dimplot with control and repopulated 
Idents(seu) <- "celltypes"
DimPlot(seu, split.by = "condition", ncol = 2 , cols = alpha(colorpalette, 0.5), 
        label = TRUE, label.size = 3.5)

#celltype contribution
dittoBarPlot(seu,"celltypes",group.by = "condition", color.panel = colorpalette, x.reorder = c(2,1,4,3),
             var.labels.reorder = c(4,1,5,2,6,3,7,8), 
             main = "",xlab = "Condition", ylab = "Percentage of cells") + theme(axis.title.x = element_text(size = 15),
                                                                                 axis.text.x = element_text(size = 12, colour = "black"),
                                                                                 axis.title.y = element_text(size = 15),
                                                                                 axis.text.y = element_text(size = 12),
                                                                                 legend.text = element_text(size = 12))
#MG3 reactome pathway
#Reactome pathway of pre-MG
#search for the genes which identifies activated cell types
celltypes_all <- FindAllMarkers(seu, logfc.threshold = 0.25,  min.pct = 0.25, test.use = "negbinom")

MG3 <- as.data.frame(rownames(celltypes_all))
MG3 <- cbind(MG3, celltypes_all$avg_log2FC, celltypes_all$p_val_adj)
names(MG3) <- c('Gene.symbol','logFC','adj.P.Val')
MG3$Gene.symbol <- toupper(MG3$Gene.symbol)
MG3_Reactome <- run_pathfindR(MG3, gene_sets = 'Reactome')

mg3_10 <- MG3_Reactome[which(MG3_Reactome$ID == 'R-HSA-72613'),]
mg3_10<- rbind(mg3_10, MG3_Reactome[which(MG3_Reactome$ID == 'R-HSA-927802'),],
               MG3_Reactome[which(MG3_Reactome$ID == 'R-HSA-9711097'),],
               MG3_Reactome[which(MG3_Reactome$ID == 'R-HSA-168255'),],
               MG3_Reactome[which(MG3_Reactome$ID == 'R-HSA-72163'),],
               MG3_Reactome[which(MG3_Reactome$ID == 'R-HSA-69610'),],
               MG3_Reactome[which(MG3_Reactome$ID == 'R-HSA-169911'),],
               MG3_Reactome[which(MG3_Reactome$ID == 'R-HSA-450531'),],
               MG3_Reactome[which(MG3_Reactome$ID == 'R-HSA-1169091'),],
               MG3_Reactome[which(MG3_Reactome$ID == 'R-HSA-1234174'),])

#chromatin remodeling markers from MG3 GO-BP expression in heatmap
chromatin_rem <- c('Actl6a', 'Hdac2', 'Myc', 'Smarcb1', 'Ruvbl1', 'Mcrs1', 'Ruvbl2', 'Nudt5', 'Uchl5', 
                   'Pole3', 'Gatad2a', 'Pih1d1', 'Smarcad1')


seu.sce <- as.SingleCellExperiment(seu)
condition_mean <- aggregateAcrossCells(as(seu.sce, "SingleCellExperiment"),  
                                       ids = seu.sce$condition, 
                                       statistics = "mean",
                                       use.assay.type = "counts")


dittoHeatmap(condition_mean,
             assay = "counts", 
             genes = chromatin_rem,
             annot.by = "condition",
             scale = "none",
             heatmap.colors = rev(brewer.pal(11,'RdBu')),
             annot.colors = condition_colors)



#extract counts
data.input <- seu@assays$RNA@data
## re-do it separately for repopulated young condition
meta <- seu@meta.data
cell.use <- rownames(meta)[meta$condition == "repopulated old"]

data.input <- data.input[,cell.use]
meta <- meta[cell.use,]

#create cellchat object
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltypes")
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) 
future::plan("multisession", workers = 4) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2


##### Figure 4 #####
#Volcano plot from bluk RNAseq data analysis 
EnhancedVolcano(top.table,
                lab = top.table$external_gene_name,
                x = 'logFC',
                y = 'P.Value',
                FCcutoff = 1)  

#enrichment analysis using pathfindR
enrichment_chart(Kegg, top_terms = 15, 
                 num_bubbles = 5) + 
  scale_color_gradientn(colors=viridis::viridis(n=15)) + 
  theme(axis.text.y =element_text(size=12))


#Heatmap of 100 expressed genes with the highest and lowest logFC
pheatmap(mat = log10(dge$counts+1),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = T,
         annotation_row = dge$logFC,
         annotation_colors = list(logFC = rev(brewer.pal(5, "Spectral"))),
         labels_col =c("ctrl", "ctrl", "rep", "rep", "rep"),
         color = rev(brewer.pal(9,"RdBu")))


#apoptosis markers from Kegg pathway activated in bulk RNAseq enrichment analysis showed in scRNAseq data using heatmap
apoptosis <- c('Apip', 'Bak1', 'Bcl2', 'Bcl2l1', 'Casp8', 'Cdkn2a', 'Dapk3', 'Dnm1l', 'Gsn', 'Lmna', 'Lmnb1', 
               'Ppp3r1', 'Psmc1', 'Psmd10', 'Psmd3', 'Psmd9', 'Psme3', 'Tfdp1')


dittoHeatmap(condition_mean,
             assay = "counts", 
             genes = apoptosis,
             annot.by = "condition",
             scale = "none",
             heatmap.colors = rev(brewer.pal(11,'RdBu')),
             annot.colors = condition_colors)


#Feature Density Plot of apoptosome markers
apoptosome <- c('Bik','Bad','Bid','Bax','Bcl2','Bcl2l1','Cycs','Apaf1','Casp9','Fas','Casp8')

Plot_Density_Joint_Only(seu, features = apoptosome, custom_palette = brewer.pal(5,"Purples")) + 
  ggtitle("Bik, Bad, Bid, Bax, Bcl2, Bcl2l1, Cycs, Apaf1, Casp9, Fas, Casp8") + 
  theme(plot.title = element_text(size=14, face="italic"))


#Ridge plot of neurodegenerative processes markers in microglia
aging <- c('Spp1', 'Igf1',  'Lamp1', 'Apoe', 'Cd63', 'Cst7')
multi <- c('Apoe', 'Lpl','Trem2','Cd68','Cd81','Cebpa','Hmgb1','Hspd1','Snx3')
alzheimer <- c('Apoe', 'Lpl', 'Trem2', 'Spp1', 'Cst7', 'Ctsd')

seu <- AddModuleScore(
  object = seu,
  features = list(multi),
  name = 'Multiple sclerosis'
)

Idents(seu) <- 'condition'


a <-  RidgePlot(seu, features =  'Aging brain1', cols = condition_colors) + 
  labs(title = "Aging brain", subtitle = 'Spp1, Igf1, Lamp1, Apoe, Cd63, Cst7') + 
  theme(plot.subtitle = element_text(face = "italic"), legend.position = 'none', 
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

b <- RidgePlot(seu, features =  'Alzheimer1', cols = condition_colors) + 
  labs(title = "Alzheimer", subtitle = 'Apoe, Lpl, Trem2, Spp1, Cst7, Ctsd') + 
  theme(plot.subtitle = element_text(face = "italic"), legend.position = 'none',
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

c <- RidgePlot(seu, features =  'Multiple sclerosis1', cols = condition_colors) + 
  labs(title = "Multiple sclerosis", subtitle = 'Apoe, Lpl, Trem2, Cd68, Cd81, Cebpa, Hmgb1, Hspd1,Snx3') + 
  theme(plot.subtitle = element_text(face = "italic"), legend.position = 'none',
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

ggarrange(a, b, c, ncol = 3)



