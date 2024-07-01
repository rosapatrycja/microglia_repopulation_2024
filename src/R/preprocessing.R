# Repopulation of microglia after BLZ depletion - Females

# Demultiplex hashed samples
# Experimental design:
# Rep1 CTRL Y1 (Ab #1) BLZ+7 Y1 (Ab #2) CTRL O1 (Ab #3) BLZ+7 O1 (Ab #4)
# Rep2 CTRL Y2 (Ab #1) BLZ+7 Y2 (Ab #2) CTRL O2 (Ab #3) BLZ+7 O2 (Ab #4)
# Rep3 CTRL Y3 (Ab #1) BLZ+7 Y3 (Ab #2) CTRL O3 (Ab #3) BLZ+7 O3 (Ab #4)
# Rep4 CTRL Y4 (Ab #1) BLZ+7 Y4 (Ab #2) CTRL O4 (Ab #3) BLZ+7 O4 (Ab #4)

# Sample CTRL O4 was mistakenly loaded together with Rep3 samples. Thus, resulting library – Rep3,
# consist of two samples tagged with same hash antibody – CTRL O3 and CTRL O4, that will be
# indistinguishable and thus will form one pooled replicate.

# Ab#1 = HTO1; Ab#2 = HTO2; Ab#3 = HTO3; Ab#4 = HTO4

getwd()
setwd("~/BLZ_MG_repopulation/sc.mm.mg.blz/counts/")

library(parallel)
library(future)
library(Seurat)
library(Matrix)
options(future.globals.maxSize=64*1024^3)
plan(multicore)
library(ggplot2)
library(cowplot)
library(magrittr)
library(dplyr)
library(purrr)
library(multtest)
library(metap)
library(limma)
library(biomaRt)
library(ape)
library(ggtree)
# Analysis parameters
mc <- 32 # number of given cores
pca_dims <- 30
umap_n_neighbors <- 30
clustering_resolution <- 0.5
n_features <- 2000


source("../analysis/functions_demultiplex.R")

# Read in a list of cell cycle markers
# Use Seurat built-in list
ccGenes <- cc.genes.updated.2019
# Change gene names to mm (Mus musculus) gene names
ccGenes <- mclapply(ccGenes, function(x)  {paste0(substring(x, 1, 1), tolower(substring(x, 2, nchar(x))))}, mc.cores = mc)

# Read raw data (gene/cell count matrix from cellranger, filtered: use only detected cellular barcodes)
samples <- dir(path = ".",
               pattern = "filtered_feature_bc_matrix$",
               full.names = T,
               recursive = T,
               include.dirs = T)
samplesRawData <- readRawData(samples = samples)
names(samplesRawData) <- substr(samples, 3, nchar(samples) - 27)

# Split input data into RNA and HTO (hashtag oligo) counts
rnaCounts <- mclapply(samplesRawData, function(x) x[grepl("ENSMUSG", rownames(x)), ], mc.cores = mc)
htosCounts <- mclapply(samplesRawData, function(x) x[grepl("HTO", rownames(x)), ], mc.cores = mc)

# Get annotation (mapping of gene names to Ensembl ID) from samples genes identifiers
annot <- getAnnot(samples)

# Map ccGenes to ENSMUSG (Ensembl id)
# Divide cell cycle genes list into markers of G2/M phase and markers of S phase
sGenes <- annot[annot$Gene.name %in% ccGenes$s.genes, "Gene.stable.ID"]
sGenes <- sGenes[!is.na(sGenes)]
g2mGenes <- annot[annot$Gene.name %in% ccGenes$g2m.genes, "Gene.stable.ID"]
g2mGenes <- g2mGenes[!is.na(g2mGenes)]

#tables with marker genes 
# Define selected (previously reported) microglia and macrophages markers
microglia_markers <- annot[match(c("Tmem119", "P2ry12", "Sall1", "Pros1", "Crybb1", "Cx3cr1", "P2ry13" ), annot$Gene.name), 
                           "Gene.stable.ID"]
macrophages_markers <- annot[match(c("Itga4", "Tgfbi", "Cxcl2", "Ccr2", "Il10", "Fgr", "Apoc1","Itgax","Cd163", "F13a1", "Cd14", "Il4ra"),
                                   annot$Gene.name),"Gene.stable.ID"]

bam_markers <- annot[match(c("Siglec1"), annot$Gene.name), "Gene.stable.ID"]

monocytes_markers <- annot[match(c("Ccr2", "Itgax", "Fcgr1", "March1", "Itga4", "Ptprc", "Cd274", "Cd33", "Itgam", "Ly6c1", "Ly6c2", "Il4r"), 
                                 annot$Gene.name), "Gene.stable.ID"]

mg1_markers <- annot[match(c("Tmem119", "Crybb1", "Cst3", "P2ry12", "Pros1"), annot$Gene.name),"Gene.stable.ID"]

mg2_markers <- annot[match(c("Jun", "Junb", "Jund", "Fos", "Egr1", "Klf6", "Aft3"), annot$Gene.name), "Gene.stable.ID"]

mg3_markers <- annot[match(c("Bmp2k", "Bhlhe41", "Ncoa3", "Tram1", "Notch2"), annot$Gene.name), "Gene.stable.ID"]

mg456_markers <- annot[match(c("Plp1", "Pltp", "Mbp", "Cd63", "Cd9"), annot$Gene.name), "Gene.stable.ID"]

premg_markers <- annot[match(c("Tmem119", "P2ry12", "Crybb1", "Csf1", "Mcm5", "Ifit3", "Cst7", "Mif", "Ccl12", "Ccl3", "Ccl4", "Ifit1", "Ifit3", "Ifit3b", "Ifitm3", "Irf7", "Isg15", "Usp18"), 
                             annot$Gene.name), "Gene.stable.ID"]


# Set up Seurat objects
seuObjects <- mclapply(seq_along(rnaCounts), function(i) {
  CreateSeuratObject(counts = rnaCounts[[i]],
                     project = paste0("BLZ_MG_repopulation", names(rnaCounts)[i]))
}, mc.cores = mc)
#FIXME tutaj jest warning: brak rownames 
names(seuObjects) <- paste0("BLZ_MG_repopulation", names(rnaCounts))

# Analyze percentage of mitochondrial genes in cells and no. of genes in cells
mitoFeatures <- annot[grep(pattern = "^mt-", annot$Gene.name), ][, 1]
percentMito <- mclapply(seq_along(seuObjects), function(i) {
  mitoReads <- Matrix::colSums(x = GetAssayData(object = seuObjects[[i]],
                                                slot = 'counts')[rownames(seuObjects[[i]]) %in% mitoFeatures, ])
  totalReads <- Matrix::colSums(x = GetAssayData(object = seuObjects[[i]],
                                                 slot = 'counts'))
  pctMito <- mitoReads / totalReads
  pctMito
}, mc.cores = mc)

# Create metrics used in QC
# percent.mito
# log10GenesPerUMI
seuObjects <- mclapply(seq_along(seuObjects), function(i) {
  seuObjects[[i]][['percent.mito']] <- percentMito[[i]]
  seuObjects[[i]]
}, mc.cores = mc)
names(seuObjects) <- names(samplesRawData)

seuObjects <- mclapply(seq_along(seuObjects), function(i) {
  seuObjects[[i]][['log10GenesPerUMI']] <- log10(seuObjects[[i]]$nFeature_RNA) / log10(seuObjects[[i]]$nCount_RNA)
  seuObjects[[i]]
}, mc.cores = mc)
names(seuObjects) <- names(samplesRawData)

# Add HTO data as a new assay independent from RNA
for(i in 1:length(seuObjects)) {
  seuObjects[[i]][["HTO"]] <- CreateAssayObject(counts = htosCounts[[i]])
}

# Normalize HTO data and demultiplex replicates into samples
for(i in 1:length(seuObjects)) {
  # Ref: https://satijalab.org/seurat/archive/v3.1/hashing_vignette.html
  # Normalize HTO data, here we use centered log-ratio (CLR) transformation
  seuObjects[[i]] <- NormalizeData(seuObjects[[i]], assay = "HTO", normalization.method = "CLR")
  
  # If you have a very large dataset we suggest using k_function = 'clara'. This is a k-medoid
  # clustering function for large applications You can also play with additional parameters (see
  # documentation for HTODemux()) to adjust the threshold for classification Here we are using the
  # default settings
  seuObjects[[i]] <- MULTIseqDemux(seuObjects[[i]], assay = "HTO", autoThresh=TRUE)
}

# Global classification results
mclapply(seuObjects, function(x) {table(x$MULTI_ID)}, mc.cores = mc)



# Remove negative and doublet cells from the object
seuObjectsAllCells <- seuObjects
seuObjects <- mclapply(seuObjects, function(x) subset(x, idents = "Negative", invert = TRUE), mc.cores = mc)
seuObjects <- mclapply(seuObjects, function(x) subset(x, idents = "Doublet", invert = TRUE), mc.cores = mc)
names(seuObjects) <- names(seuObjects)
names(seuObjects) <- paste0(names(seuObjects), "_", "singlets")

#setwd into that directory, where all computations will be saved
setwd("~/BLZ_analysis")
saveRDS(seuObjects, file = "results/rds/seuObject.rds")
saveRDS(seuObjectsAllCells, file = "results/rds/seuObjectsAllCells.rds")

#add grouping factor
seuObjects[[1]]$sample <- paste('rep1',seuObjects[[1]]$MULTI_ID, sep="_")
seuObjects[[2]]$sample <- paste('rep2',seuObjects[[2]]$MULTI_ID, sep="_")
seuObjects[[3]]$sample <- paste('rep3',seuObjects[[3]]$MULTI_ID, sep="_")
seuObjects[[4]]$sample <- paste('rep4',seuObjects[[4]]$MULTI_ID, sep="_")

seuObjects[[1]]$replicate <- 'rep1'
seuObjects[[2]]$replicate <- 'rep2'
seuObjects[[3]]$replicate <- 'rep3'
seuObjects[[4]]$replicate <- 'rep4'

# Add sample shortID
seuObjects[[1]]$shortID <- substr(as.character(seuObjects[[1]]$orig.ident), 15,
                                         nchar(as.character(seuObjects[[1]]$orig.ident)))
seuObjects[[2]]$shortID <- substr(as.character(seuObjects[[2]]$orig.ident), 15,
                                         nchar(as.character(seuObjects[[2]]$orig.ident)))
seuObjects[[3]]$shortID <- substr(as.character(seuObjects[[3]]$orig.ident), 15,
                                         nchar(as.character(seuObjects[[3]]$orig.ident)))
seuObjects[[4]]$shortID <- substr(as.character(seuObjects[[4]]$orig.ident), 15,
                                         nchar(as.character(seuObjects[[4]]$orig.ident)))

# Add selected (previously reported) genes to var.features
seuObjects[[1]]@assays$RNA@var.features <- unique(c(seuObjects[[1]]@assays$RNA@var.features,
sGenes, g2mGenes))
seuObjects[[2]]@assays$RNA@var.features <- unique(c(seuObjects[[2]]@assays$RNA@var.features,
sGenes, g2mGenes))
seuObjects[[3]]@assays$RNA@var.features <- unique(c(seuObjects[[3]]@assays$RNA@var.features,
sGenes, g2mGenes))
seuObjects[[4]]@assays$RNA@var.features <- unique(c(seuObjects[[4]]@assays$RNA@var.features,
sGenes, g2mGenes))



#filter out clusters with high number of mitochondrial genes
#normalize data
#Find variable features 

seuObjects <-  mclapply(X = seuObjects, FUN = function(x) {
  x <- subset(x, subset = percent.mito < 0.075)
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = n_features)
}, mc.cores = mc)


#save filtered out object on new variable
seuObjects_hash <- seuObjects


#combine to the one object
seu.merged <- merge(x = seuObjects[[1]], 
                    y = seuObjects[2:length(seuObjects)])

#Cell Cycle scoring
seu.merged<- CellCycleScoring(object = seu.merged, s.features = sGenes, 
                                           g2m.features = g2mGenes, set.ident = TRUE)

# Calculate difference between G2M and S phase scores to separate non-cycling and cycling cells
# Approch described in Seurat's Cell-Cycle Scoring and Regression vignette 
seu.merged$CC_Difference <- seu.merged$S.Score - seu.merged$G2M.Score

    
seu.merged <- FindVariableFeatures(seu.merged, selection.method = "vst", nfeatures = n_features)
seu.merged@assays$RNA@var.features <- unique(c(seu.merged@assays$RNA@var.features, 
                                                             sGenes, g2mGenes))

#scale data to plot umap
all.genes <- rownames(seu.merged)
seu.merged <- ScaleData(object = seu.merged,features = all.genes)


# Run the standard workflow for visualization and clustering
seu.merged <- RunPCA(seu.merged, npcs = pca_dims, verbose = FALSE)
seu.merged <- RunUMAP(seu.merged, reduction = "pca", dims = 1:30)
seu.merged <- FindNeighbors(seu.merged, reduction = "pca", dims = 1:umap_n_neighbors)
seu.merged <- FindClusters(seu.merged, resolution = clustering_resolution)

#remove that one replicate which was mislead
seu.merged<- seu.merged[,which(seu.merged$sample != "rep4_HTO2")]

#change numbers of clusters from 0 to 1 
seu.merged$seurat_clusters <- as.factor(as.numeric(as.character(seu.merged$seurat_clusters)) + 1)
Idents(seu.merged) <- seu.merged$seurat_clusters

#rename MULTI_ID into condition
condition <- seu.merged$MULTI_ID
condition <- sapply(condition, function(x) gsub('HTO1', 'control young',x))
condition <- sapply(condition, function(x) gsub('HTO2', 'control old',x))
condition <- sapply(condition, function(x) gsub('HTO3', 'repopulated young',x))
condition <- sapply(condition, function(x) gsub('HTO4', 'repopulated old',x))

#order the data by derived conditions
condition<- factor(condition, levels = c("control young", "control old", "repopulated young", "repopulated old"))

#merge it to the data
seu.merged@meta.data$condition <- condition

saveRDS(seu.merged, file = "results/rds/seumerged.rds")

#find markers for every cluster compared to all remaining cells, report only the positive
# ones
seu.markers <- FindAllMarkers(object = seu.merged, logfc.threshold = 0.25,  min.pct = 0.25, only.pos = TRUE)
top20genes <- seu.markers %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC)

colnames(top20genes) <- c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "cluster", "Gene.stable.ID")

top20_merged <- merge(top20genes, annot, by = "Gene.stable.ID")
top20_merged <- top20_merged[,c(8,1,7,2,6,3,4,5)]
top20_merged <- top20_merged[order(top20_merged$cluster, top20_merged$avg_log2FC),]


write.table(top20_merged, file="results/tables/top20genestable.tsv", quote = F, sep = "\t", row.names = FALSE, col.names = TRUE)

#filter out unnecesary clusters
clusters_to_take <- c(1,2,3,4,5,6,7,8,9,11,14)
seu.merged <- seu.merged[,which(seu.merged$seurat_clusters %in% clusters_to_take)]
seu.merged$seurat_clusters <- droplevels(seu.merged$seurat_clusters)

#annotate clusters using manual annotations 
new_cluster_ids <- c("MG1", "Act-MG1", "MG1", "MG2", "Act-MG3", "Act-MG2", "prolif-MG", 
                     "MG3", "Act-MG3", "MG1", "pre-MG")
names(new_cluster_ids) <- levels(seu.merged)

#rename cluster numbers with annotations
seu.merged <- RenameIdents(seu.merged, new_cluster_ids)
seu.merged@meta.data$celltypes <- seu.merged@active.ident

