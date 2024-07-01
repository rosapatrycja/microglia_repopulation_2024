significant_genes <- function(qlf_data,alpha = Inf,logfc = -Inf) {
  #choose genes with abs(logfoldchange) >2 and adjusted p-value<0.05
  thresh_data <- qlf_data[which(abs(qlf_data$table$logFC)>logfc & qlf_data$table$adjPValue<alpha),]
  
  
}