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


source("../analysis/src/counts_heatmap.R")
source("../analysis/src/maxlogfc_genes.R")

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

#using BiomaRt for getting gene names
mart<-useMart(biomart = "ensembl",
              dataset = "mmusculus_gene_ensembl")


#select attributes from BiomaRt
genes_mart <- getBM(
  filters="ensembl_gene_id",
  uniqueRows = TRUE,
  attributes=c("ensembl_gene_id", "external_gene_name","description","entrezgene_id", "chromosome_name"),
  values=genes,
  mart=mart)


#Arrange matrix and filter out empty entrezid records, predicted genes and many others, which are unnecessary 
#to extractonly known gene names  
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


#filter out genes with low variance
cutoff <- 1
drop <- which(apply(edgeR::cpm(data_norm), 1, max) < cutoff)
data_norm <- data_norm[-drop,]

#remove outlier sample
data_norm <- data_norm[,-which(data_norm$samples$sample_id == "S2")]

#labels for plots
labels <- c(control = "Control", repopulated = "Repopulated")

#PCA dimension reduction
gene_pca <- prcomp(t(data_norm$counts))
pca_out <- as.data.frame(gene_pca$x)
pca_out$condition <- data_norm$samples$condition
pca_out$sample <- data_norm$samples$sample_id
percentage <- round(gene_pca$sdev / sum(gene_pca$sdev) * 100, 2)
percentage <- paste0( colnames(pca_out), " (", paste0( as.character(percentage), "%", ")") )


#scale counts to counts per million
data_cpm <- data_norm
cpms <- edgeR::cpm(data_norm)
data_cpm$counts <- cpms

#model matrix
mm <- model.matrix(~0 + data_norm$samples$condition)
colnames(mm) <- c('control', 'repopulated')

#transform to cpm, fit linear model and smooth the curve 
y <- voom(data_norm, mm, plot = T)


#fitting linear models in limma
fit <- lmFit(y, mm)
head(coef(fit))

#make contrast
contr <- makeContrasts(repopulated - control, levels = colnames(coef(fit)))
contr

#estimate contrast for each gene
tmp <- contrasts.fit(fit, contr)

tmp <- eBayes(tmp)


#diff expressed genes summed in table 
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)


#how many DE are here 
length(which(top.table$adj.P.Val < 0.05))

tt1 <- top.table[which(top.table$adj.P.Val < 0.05),]
tt_up <- tt1[which(tt1$logFC > 0.5),]
tt_down <- tt1[which(tt1$logFC < -0.5),]



#filter genes with the lowest, hisghets and combined low and high values of logFC
dge <- maxlogfc_genes(tt1,50,data_cpm)
up <- maxlogfc_genes(tt_up,50,data_norm)
down <- maxlogfc_genes(tt_down, 50, data_norm)



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



