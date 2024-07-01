fdr_adjust <- function(qlf_data) {
  
  #multiple testing correction method
  FDR.data <- p.adjust(qlf_data$table$PValue, method = "BH")
  #add p-adjusted value into summary table 
  qlf_data$table$adjPValue <- FDR.data
  qlf_data
}