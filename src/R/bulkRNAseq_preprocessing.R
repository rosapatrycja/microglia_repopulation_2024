library(edgeR)
library(biomaRt)
library(ggplot2)
library(RColorBrewer)
library(limma)
library(statmod)
library(VennDiagram)
library(pheatmap)
library(EnhancedVolcano)
library(tmod)
library(cowplot)
library(tidyverse)


source("/mnt/wspolny/patrycja.rosa/BLZ_RNAseq/analysis/src/counts_heatmap.R")
source("/mnt/wspolny/patrycja.rosa/BLZ_RNAseq/analysis/src/maxlogfc_genes.R")

# for plotting, take the genes with abs(log2 fold change) > 2 and adjusted p-value < 0.05
logfc = 1
alpha = 0.05 

setwd("/mnt/wspolny/patrycja.rosa/BLZ_RNAseq")

#read counts file 
cnt <- read.table("results/counts/all.tsv", header = TRUE, quote = "\"'")
rownames(cnt) <- cnt$gene

counts <- cnt[,-1]

#add some metadata
conditions <- c("control", "control", "control", "repopulated", "repopulated","repopulated")
conditions <- cbind(colnames(counts),conditions)
colnames(conditions) <- c("sample_id", "condition")
rownames(conditions) <- colnames(counts)


#genes
genes <- cnt$gene
genes <- unique(genes)





#using biomart for getting gene names
mart<-useMart(biomart = "ensembl",
              dataset = "mmusculus_gene_ensembl")



genes_mart <- getBM(
  filters="ensembl_gene_id",
  uniqueRows = TRUE,
  attributes=c("ensembl_gene_id", "external_gene_name","description","entrezgene_id", "chromosome_name"),
  values=genes,
  mart=mart)


rownames(genes_mart) = make.names(genes_mart$ensembl_gene_id, unique=TRUE)
genes_mart <- genes_mart[!(genes_mart$external_gene_name == ""),]
genes_mart <- genes_mart[!(is.na(genes_mart$entrezgene_id)), ]
genes_mart$description <- gsub("predicted.*", "", genes_mart$description, perl = TRUE)
genes_mart$description <- gsub("RIKEN.*", "", genes_mart$description, perl = TRUE)
genes_mart <- genes_mart[!(genes_mart$description == ""), ]
genes_mart$external_gene_name <- gsub("Mir.*", "", genes_mart$external_gene_name, perl = TRUE)
genes_mart$external_gene_name <- gsub("Zfp.*", "", genes_mart$external_gene_name, perl = TRUE)
genes_mart$external_gene_name <- gsub("Vmn.*", "", genes_mart$external_gene_name, perl = TRUE)
genes_mart$external_gene_name <- gsub("Slc.*", "", genes_mart$external_gene_name, perl = TRUE)
genes_mart$external_gene_name <- gsub("Rik.*", "", genes_mart$external_gene_name, perl = TRUE)
genes_mart <- genes_mart[!(genes_mart$external_gene_name == ""), ]

#counts
counts <- counts[which(rownames(counts) %in% rownames(genes_mart)),]

#genes

genes_mart <- genes_mart[order(which(rownames(genes_mart) %in% rownames(counts))),]

#make dge object for further analysis
dge_object <- DGEList(counts = counts, samples = conditions, genes = genes_mart)

#change rownames in counts and genes
rownames(dge_object$counts) <- dge_object$genes$external_gene_name


dge_object$samples$lib.size <- colSums(dge_object$counts)


#normalization
data_norm<-calcNormFactors(dge_object)

#counts per million
#filter out genes with low variance
cutoff <- 1
drop <- which(apply(edgeR::cpm(data_norm), 1, max) < cutoff)
data_norm <- data_norm[-drop,]

#remove outlier sample
data_norm <- data_norm[,-which(data_norm$samples$sample_id == "S2")]

#labels for plots
labels <- c(control = "Control", repopulated = "Repopulated")

#PCA
gene_pca <- prcomp(t(data_norm$counts))
pca_out <- as.data.frame(gene_pca$x)
pca_out$condition <- data_norm$samples$condition
pca_out$sample <- data_norm$samples$sample_id
percentage <- round(gene_pca$sdev / sum(gene_pca$sdev) * 100, 2)
percentage <- paste0( colnames(pca_out), " (", paste0( as.character(percentage), "%", ")") )
pc<-ggplot(pca_out,aes(x=PC1,y=PC2, color= condition))
pc<-pc+geom_point(size=4)+scale_color_brewer(palette= 'Set1', labels=labels)+ xlab(percentage[1]) + ylab(percentage[2])
pc
pc2<-ggplot(pca_out,aes(x=PC1,y=PC2, color= sample))
pc2<-pc2+geom_point(size=4)+scale_color_brewer(palette= 'Set1')+ xlab(percentage[1]) + ylab(percentage[2])
pc2

#model matrix
mm <- model.matrix(~0 + data_norm$samples$condition)
colnames(mm) <- c('control', 'repopulated')

#transform to cpm, fit linear model and smooth the curve 
y <- voom(data_norm, mm, plot = T)


#Fitting linear models in limma
fit <- lmFit(y, mm)
head(coef(fit))

#make contrast
contr <- makeContrasts(repopulated - control, levels = colnames(coef(fit)))
contr

#estimate contrast for each gene
tmp <- contrasts.fit(fit, contr)

tmp <- eBayes(tmp)


#diff expressed genes
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)


#how many DE are here 
length(which(top.table$adj.P.Val < 0.05))

tt1 <- top.table[which(top.table$adj.P.Val < 0.05),]

tt_up <- tt1[which(tt1$logFC > 0.5),]

tt_down <- tt1[which(tt1$logFC < -0.5),]


#volcano
EnhancedVolcano(top.table,
                lab = top.table$external_gene_name,
                x = 'logFC',
                y = 'P.Value',
                FCcutoff = 1)              

#show gene expression
data_cpm <- data_norm
cpms <- edgeR::cpm(data_norm)
data_cpm$counts <- cpms

#pheatmap
dge_data <- data_norm[which(rownames(data_norm) %in% tt1$external_gene_name),]
counts_heatmap(dge_data$counts, cluster = TRUE, rownames_value = F)

data_up <- data_norm[which(rownames(data_norm) %in% tt_up$external_gene_name),]
counts_heatmap(data_up$counts, cluster = TRUE, rownames_value = F)

data_down<- data_norm[which(rownames(data_norm) %in% tt_down$external_gene_name),]
counts_heatmap(data_down$counts, cluster = TRUE, rownames_value = F)



dge <- maxlogfc_genes(tt1,50,data_cpm)
up <- maxlogfc_genes(tt_up,50,data_norm)
down <- maxlogfc_genes(tt_down, 50, data_norm)


pheatmap(mat = log10(dge$counts+1),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = T,
         annotation_row = dge$logFC,
         annotation_colors = list(logFC = rev(brewer.pal(5, "Spectral"))),
         labels_col =c("ctrl", "ctrl", "rep", "rep", "rep"),
         color = rev(brewer.pal(9,"RdBu")))





#gene ontologies

#create dataframe with values needed to run pathfindR
diffexpgenes <- as.data.frame(rownames(top.table))
diffexpgenes <- cbind(diffexpgenes, top.table$external_gene_name, top.table$logFC, top.table$adj.P.Val)
names(diffexpgenes) <- c('Entrezid','Gene.symbol', 'logFC', 'adj.P.Val')

diffexpgenes$Gene.symbol <- toupper(diffexpgenes$Gene.symbol)
diffexpgenes <- diffexpgenes[,-1]
#run enrichment analysis 
Kegg <- run_pathfindR(diffexpgenes, pin_name_path = "KEGG")
GO <- run_pathfindR(diffexpgenes, 
                         gene_sets = 'GO-BP')
Reactome <- run_pathfindR(diffexpgenes, gene_sets = 'Reactome')



enrichment_chart(Kegg, top_terms = 15, 
                 num_bubbles = 5) + 
  scale_color_gradientn(colors=viridis::viridis(n=15)) + 
  theme(axis.text.y =element_text(size=12))




d1 <- as.matrix(data_norm[,])
df_2 <- as.data.frame(t(d1))
df_2$group <- data_norm$samples$condition
df_2



df_2 %>%
  pivot_longer(cols = -group, names_to = "genes") %>%
  filter(genes == "Bid") %>%
  ggplot(aes(x = group, y = value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = genes), width = 0.2) +
  theme_minimal()
