library(Seurat)
library(Matrix)
library(ggplot2)
library(viridis)

setwd("/home/rstudio/BLZ_analysis")


seurat_obj <- readRDS('data/seu_06072023.RDS')
out_data_dir <- "/home/rstudio/BLZ_analysis/data"


# Idents(seurat_obj)
# new_cluster_ids <- c("MG1", "Act-MG1", "MG2", "Act-MG3", "Act-MG2", "prolif-MG", 
#                      "MG3",  "pre-MG")
# 
# names(new_cluster_ids) <- levels(seurat_obj)
# 
# #rename cluster numbers with annotations
# seu <- RenameIdents(seurat_obj, new_cluster_ids)
# seu@meta.data$celltypes <- seu@active.ident
# 
# 
# seurat_obj$celltypes <- factor(seurat_obj$celltypes, 
#                         levels = c("MG1", "Act-MG1", "MG2", "Act-MG2", "MG3", "Act-MG3", "pre-MG", "prolif-MG"))
# 


# save metadata table:
seurat_obj$barcode <- colnames(seurat_obj)
seurat_obj$UMAP_1 <- seurat_obj@reductions$umap@cell.embeddings[,1]
seurat_obj$UMAP_2 <- seurat_obj@reductions$umap@cell.embeddings[,2]




write.csv(seurat_obj@meta.data, file='data/metadata.csv', quote=F, row.names=F)

# write expression counts matrix

counts_matrix <- GetAssayData(seurat_obj, assay='RNA', slot='counts')
writeMM(counts_matrix, file='data/counts.mtx')

# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(seurat_obj@reductions$pca@cell.embeddings, file='data/pca.csv', quote=F, row.names=F)


# write gene names
write.table(
  data.frame('gene'=rownames(counts_matrix)),file='data/gene_names.csv',
  quote=F,row.names=F,col.names=F
)



