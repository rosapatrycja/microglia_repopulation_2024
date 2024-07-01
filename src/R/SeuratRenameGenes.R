library(Seurat)

SeuratRenameGenes <- function(obj, names_from, names_to, assays = "RNA") {

  #' Replace gene names in different slots of a Seurat object
  #'
  #' Returns the modified Seurat object
  #' @param obj Seurat object
  #' @param names_from Vector of gene names to be replaced
  #' @param names_to Vector of gene name replacements
  #' @param assays Vector of assays; "RNA" by default

  names_from <- as.character(names_from)
  names_to <- as.character(names_to)
  stopifnot(length(names_from) == length(names_to))

  # replace gene names in a character vector
  replace_genes <- function(v)
  {
    sel <- match(v, names_from)
    v <- ifelse(is.na(sel), v, names_to[sel])
    return(v)
  }

  for (assay in assays)
  {
    as <- obj@assays[[assay]]

    if (.hasSlot(as, "counts"))
      rownames(as@counts) <- replace_genes(rownames(as@counts))
    if (.hasSlot(as, "data"))
      rownames(as@data) <- replace_genes(rownames(as@data))
    if (.hasSlot(as, "scale.data"))
      rownames(as@scale.data) <- replace_genes(rownames(as@scale.data))
    if (.hasSlot(as, "var.features"))
      as@var.features <- replace_genes(as@var.features)

    obj@assays[[assay]] <- as
  }

  if (!is.null(obj@reductions$pca))
    row.names(obj@reductions$pca@feature.loadings) <-
      replace_genes(row.names(obj@reductions$pca@feature.loadings))

  return(obj)
}


#
#  Test run
#

# gene_names <- read.table("data/gene_names.tsv", sep = "\t", header = T)
# seu <- readRDS("/home/rosa/BLZ_analysis/data/seumerged.rds")
# seu <- SeuratRenameGenes(seu, gene_names$gene_id, gene_names$gene_name)
