library(Seurat)
library(ggplot2)
library(pathfindR)
source("/mnt/wspolny/patrycja.rosa/BLZ_analysis/src/R/SeuratRenameGenes.R")

setwd('/mnt/wspolny/patrycja.rosa/BLZ_analysis')
seu <- readRDS('data/seumerged_plots.rds')


#filter out unnecesary clusters
clusters_to_take <- c(1,2,3,4,5,6,7,8,9,11,14)
seu<- seu[,which(seu$seurat_clusters %in% clusters_to_take)]
seu$seurat_clusters <- droplevels(seu$seurat_clusters)

#annotate clusters using manual annotations 
new_cluster_ids <- c("MG1", "Act-MG1", "MG1", "MG2", "Act-MG3", "Act-MG2", "prolif-MG", 
                     "MG3", "Act-MG3", "MG1", "pre-MG")
names(new_cluster_ids) <- levels(seu)

#rename cluster numbers with annotations
seu <- RenameIdents(seu, new_cluster_ids)
seu@meta.data$celltypes <- seu@active.ident


Idents(young) <- "celltypes"
#dge expressed genes
MG1_young <- FindMarkers(object = young, ident.1 = "repopulated", ident.2 = "control", 
                         group.by = "condition", subset.ident = c("MG_1", "Act-MG_1"), test.use = 'negbinom')
MG2_young <-  FindMarkers(object = young, ident.1 = "repopulated", ident.2 = "control", 
                          group.by = "condition", subset.ident = c("MG_2", "Act-MG_2"), test.use = 'negbinom')
MG3_young <-  FindMarkers(object = young, ident.1 = "repopulated", ident.2 = "control", 
                          group.by = "condition", subset.ident = c("MG_3", "Act-MG_3"), test.use = 'negbinom')

#pathfindR
MG1 <- as.data.frame(rownames(MG1_young))
MG1 <- cbind(MG1, MG1_young$avg_log2FC, MG1_young$p_val_adj)
names(MG1) <- c('Gene.symbol','logFC','adj.P.Val')
MG1$Gene.symbol <- toupper(MG1$Gene.symbol)
MG1_KEGG <- run_pathfindR(MG1, gene_sets = 'KEGG')
MG1_Reactome <- run_pathfindR(MG1, gene_sets = 'Reactome')
#enrichment_chart(MG1_Reactome, top_terms = 50)
MG1_Reactome$celltype <- 'MG1'
MG1_Reactome$condition <- 'young'

#pathfindR
MG2 <- as.data.frame(rownames(MG2_young))
MG2 <- cbind(MG2, MG2_young$avg_log2FC, MG2_young$p_val_adj)
names(MG2) <- c('Gene.symbol','logFC','adj.P.Val')
MG2$Gene.symbol <- toupper(MG2$Gene.symbol)
MG2_Reactome <- run_pathfindR(MG2, gene_sets = 'Reactome')
#enrichment_chart(MG2_Reactome, top_terms = 50)
MG2_Reactome$celltype <- 'MG2'
MG2_Reactome$condition <- 'young'

#pathfindR
MG3 <- as.data.frame(rownames(MG3_young))
MG3 <- cbind(MG3, MG3_young$avg_log2FC, MG3_young$p_val_adj)
names(MG3) <- c('Gene.symbol','logFC','adj.P.Val')
MG3$Gene.symbol <- toupper(MG3$Gene.symbol)
MG3_Reactome <- run_pathfindR(MG3, gene_sets = 'Reactome')
#enrichment_chart(MG3_Reactome, top_terms = 50)
MG3_Reactome$celltype <- 'MG3'
MG3_Reactome$condition <- 'young'

#search for the genes which identifies activated cell types
celltypes <- FindAllMarkers(seu, logfc.threshold = 0.25,  min.pct = 0.25, test.use = "negbinom")

Act_MG1 <- celltypes[which(celltypes$cluster == "Act-MG_1"),]
Act_MG2 <- celltypes[which(celltypes$cluster == "Act-MG_2"),]
Act_MG3 <- celltypes[which(celltypes$cluster == "Act-MG_3"),]


Act_MG1_df <- as.data.frame(rownames(Act_MG1))
Act_MG1_df <- cbind(Act_MG1_df, Act_MG1$avg_log2FC, Act_MG1$p_val_adj)
names(Act_MG1_df) <- c('Gene.symbol','logFC','adj.P.Val')
Act_MG1_df$Gene.symbol <- toupper(Act_MG1_df$Gene.symbol)
Act_MG1_KEGG <- run_pathfindR(Act_MG1_df, gene_sets = 'KEGG')
Act_MG1_GO <- run_pathfindR(Act_MG1_df, gene_sets = 'GO-BP')
Act_MG1_Reactome <- run_pathfindR(Act_MG1_df, gene_sets = 'Reactome')
#enrichment_chart(Act_MG1_Reactome, top_terms = 50)
Act_MG1_Reactome$celltype <- 'Act_MG1'
Act_MG1_Reactome$condition <- 'young'

Act_MG2_df <- as.data.frame(rownames(Act_MG2))
Act_MG2_df <- cbind(Act_MG2_df, Act_MG2$avg_log2FC, Act_MG2$p_val_adj)
names(Act_MG2_df) <- c('Gene.symbol','logFC','adj.P.Val')
Act_MG2_df$Gene.symbol <- toupper(Act_MG2_df$Gene.symbol)
Act_MG2_Reactome <- run_pathfindR(Act_MG2_df, gene_sets = 'Reactome')
#enrichment_chart(Act_MG2_Reactome, top_terms = 50)
Act_MG2_Reactome$celltype <- 'Act_MG2'
Act_MG2_Reactome$condition <- 'young'


Act_MG3_df <- as.data.frame(rownames(Act_MG3))
Act_MG3_df <- cbind(Act_MG3_df, Act_MG3$avg_log2FC, Act_MG3$p_val_adj)
names(Act_MG3_df) <- c('Gene.symbol','logFC','adj.P.Val')
Act_MG3_df$Gene.symbol <- toupper(Act_MG3_df$Gene.symbol)
Act_MG3_KEGG <- run_pathfindR(Act_MG3_df, gene_sets = 'KEGG')
Act_MG3_GO <- run_pathfindR(Act_MG3_df, gene_sets = 'GO-BP')
Act_MG3_Reactome <- run_pathfindR(Act_MG3_df, gene_sets = 'Reactome')
#enrichment_chart(Act_MG3_Reactome, top_terms = 50)
Act_MG3_Reactome$celltype <- 'Act_MG3'
Act_MG3_Reactome$condition <- 'young'

#MG3
MG3_genes <- celltypes[which(celltypes$cluster == "MG_3"),]
MG3_df <- as.data.frame(rownames(MG3_genes))
MG3_df <- cbind(MG3_df, MG3_genes$avg_log2FC, MG3_genes$p_val_adj)
names(MG3_df) <- c('Gene.symbol','logFC','adj.P.Val')
MG3_df$Gene.symbol <- toupper(MG3_df$Gene.symbol)
MG3_KEGG <- run_pathfindR(MG3_df, gene_sets = 'KEGG')
MG3_GO <- run_pathfindR(MG3_df, gene_sets = 'GO-BP')
MG3all_Reactome <- run_pathfindR(MG3_df, gene_sets = 'Reactome')
#enrichment_chart(Act_MG1_Reactome, top_terms = 50)


preMG <- celltypes[which(celltypes$cluster == "pre-MG"),]

preMG_df <- as.data.frame(rownames(preMG))
preMG_df <- cbind(preMG_df, preMG$avg_log2FC, preMG$p_val_adj)
names(preMG_df) <- c('Gene.symbol','logFC','adj.P.Val')
preMG_df$Gene.symbol <- toupper(preMG_df$Gene.symbol)
preMG_Reactome <- run_pathfindR(preMG_df, gene_sets = 'Reactome')
#enrichment_chart(preMG_Reactome, top_terms = 50)
preMG_Reactome$celltype <- 'preMG'
preMG_Reactome$condition <- 'young'



prolifMG <- celltypes[which(celltypes$cluster == "prolif-MG"),]

prolifMG_df <- as.data.frame(rownames(prolifMG))
prolifMG_df <- cbind(prolifMG_df, prolifMG$avg_log2FC, prolifMG$p_val_adj)
names(prolifMG_df) <- c('Gene.symbol','logFC','adj.P.Val')
prolifMG_df$Gene.symbol <- toupper(prolifMG_df$Gene.symbol)
prolifMG_Reactome <- run_pathfindR(prolifMG_df, gene_sets = 'Reactome')
#enrichment_chart(prolifMG_Reactome, top_terms = 50)
prolifMG_Reactome$celltype <- 'prolifMG'
prolifMG_Reactome$condition <- 'young'


Idents(seu) <- "celltypes"

#dge expressed genes
MG1_old <- FindMarkers(object = seu, ident.1 = "repopulated old", ident.2 = "control old", 
                         group.by = "condition", subset.ident = c("MG_1", "Act-MG_1"), test.use = 'negbinom')
MG2_old <-  FindMarkers(object = seu, ident.1 = "repopulated old", ident.2 = "control old", 
                          group.by = "condition", subset.ident = c("MG_2", "Act-MG_2"), test.use = 'negbinom')
MG3_old <-  FindMarkers(object = seu, ident.1 = "repopulated old", ident.2 = "control old", 
                          group.by = "condition", subset.ident = c("MG_3", "Act-MG_3"), test.use = 'negbinom')


MG1_new <- as.data.frame(rownames(MG1_old))
MG1_new <- cbind(MG1_new, MG1_old$avg_log2FC, MG1_old$p_val_adj)
names(MG1_new) <- c('Gene.symbol','logFC','adj.P.Val')
MG1_new$Gene.symbol <- toupper(MG1_new$Gene.symbol)
MG1_new_Reactome <- run_pathfindR(MG1_new, gene_sets = 'Reactome')
#enrichment_chart(MG1_new_Reactome, top_terms = 50)
MG1_new_Reactome$celltype <- 'MG1'
MG1_new_Reactome$condition <- 'old'



MG2_new <- as.data.frame(rownames(MG2_old))
MG2_new <- cbind(MG2_new, MG2_old$avg_log2FC, MG2_old$p_val_adj)
names(MG2_new) <- c('Gene.symbol','logFC','adj.P.Val')
MG2_new$Gene.symbol <- toupper(MG2_new$Gene.symbol)
MG2_new_Reactome <- run_pathfindR(MG2_new, gene_sets = 'Reactome')
#enrichment_chart(MG2_new_Reactome, top_terms = 50)
MG2_new_Reactome$celltype <- 'MG2'
MG2_new_Reactome$condition <- 'old'


MG3_new <- as.data.frame(rownames(MG3_old))
MG3_new <- cbind(MG3_new, MG3_old$avg_log2FC, MG3_old$p_val_adj)
names(MG3_new) <- c('Gene.symbol','logFC','adj.P.Val')
MG3_new$Gene.symbol <- toupper(MG3_new$Gene.symbol)
MG3_new_Reactome <- run_pathfindR(MG3_new, gene_sets = 'Reactome')
#enrichment_chart(MG3_new_Reactome, top_terms = 50)
MG3_new_Reactome$celltype <- 'MG3'
MG3_new_Reactome$condition <- 'old'


#dge expressed genes
PreMG_all <- FindMarkers(object = seu, ident.1 = "repopulated old", ident.2 = "repopulated young", 
                       group.by = "condition", subset.ident = "pre-MG", test.use = 'negbinom')

ProlifMG_all <- FindMarkers(object = seu, ident.1 = "repopulated old", ident.2 = "repopulated young", 
                         group.by = "condition", subset.ident = "prolif-MG", test.use = 'negbinom')


PreMG_new <- as.data.frame(rownames(PreMG_all))
PreMG_new <- cbind(PreMG_new, PreMG_all$avg_log2FC, PreMG_all$p_val_adj)
names(PreMG_new) <- c('Gene.symbol','logFC','adj.P.Val')
PreMG_new$Gene.symbol <- toupper(PreMG_new$Gene.symbol)
PreMG_new_Reactome <- run_pathfindR(PreMG_new, gene_sets = 'Reactome')
#enrichment_chart(PreMG_new_Reactome, top_terms = 50)
PreMG_new_Reactome$celltype <- 'PreMG'
PreMG_new_Reactome$condition <- 'repopulated'

ProlifMG_new <- as.data.frame(rownames(ProlifMG_all))
ProlifMG_new <- cbind(ProlifMG_new, ProlifMG_all$avg_log2FC, ProlifMG_all$p_val_adj)
names(ProlifMG_new) <- c('Gene.symbol','logFC','adj.P.Val')
ProlifMG_new$Gene.symbol <- toupper(ProlifMG_new$Gene.symbol)
ProlifMG_new_Reactome <- run_pathfindR(ProlifMG_new, gene_sets = 'Reactome')
#enrichment_chart(ProlifMG_new_Reactome, top_terms = 50)
ProlifMG_new_Reactome$celltype <- 'ProlifMG'
ProlifMG_new_Reactome$condition <- 'repopulated'

# all_reactome <- MG1_Reactome
# all_reactome <- rbind(all_reactome, MG2_Reactome, MG3_Reactome, Act_MG1_Reactome, Act_MG2_Reactome, 
#                       Act_MG3_Reactome, preMG_Reactome, prolifMG_Reactome, MG1_new_Reactome, MG2_new_Reactome,
#                       MG3_new_Reactome, ProlifMG_new_Reactome)
# write.table(all_reactome, file = "Reactome_paths.tsv", sep = ",",row.names = FALSE)


save.image(file = "data/pathfindR_enrichment.RData")


#repo 


repo <- FindMarkers(object = seu, ident.1 = "repopulated old", ident.2 = "repopulated young")
