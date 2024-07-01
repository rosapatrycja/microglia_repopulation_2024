counts_heatmap <- function(counts_matrix, cluster, rownames_value) { 
  
  pheatmap(mat = log10(1+counts_matrix), 
           cluster_rows = cluster,
           cluster_cols = FALSE, 
           clustering_distance_rows = "euclidean",
           clustering_method = "complete", 
           show_rownames = rownames_value, 
           labels_col =c("ctrl", "ctrl", "rep", "rep", "rep"))
  
}
