maxlogfc_genes <- function(mat,N,cpm_matrix) {
  #sort genes by decreasing absolute value of logfc
  sortedgenes <- mat[order(abs(mat$logFC),decreasing = TRUE),]
  #N genes with  max logfoldchange
  sortedgenes_N <- sortedgenes[1:N,]
  #order by increasing log fold change value
  sortedgenes_N <- sortedgenes_N[order(sortedgenes_N$logFC, decreasing = TRUE), ]
  rownames(sortedgenes_N) <- sortedgenes_N$external_gene_name
  #sorted log fold change values
  sorted_logfc <- sortedgenes_N[,7,drop = FALSE]
  #counts matrix
  sortedcounts <- cpm_matrix[which(rownames(cpm_matrix) %in% sortedgenes_N$external_gene_name),]
  #sort by log fold change value
  sortedcounts <- sortedcounts[sortedgenes_N$external_gene_name,]
  return(list(counts=sortedcounts$counts,logFC=sorted_logfc))
}





