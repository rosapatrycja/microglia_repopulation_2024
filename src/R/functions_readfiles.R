# Read raw 10x counts data
readRawData <- function(samples) {
  samplesRawData <- mclapply(samples, function(s) {
  matrixDir = s
  barcodePath <- paste0(matrixDir, "/", "barcodes.tsv.gz")
  featuresPath <- paste0(matrixDir, "/", "features.tsv.gz")
  matrixPath <- paste0(matrixDir, "/", "matrix.mtx.gz")
  mat <- readMM(file = matrixPath)
  featureNames = read.delim(featuresPath,
                            header = FALSE,
                            stringsAsFactors = FALSE)
  barcodeNames = read.delim(barcodePath,
                            header = FALSE,
                            stringsAsFactors = FALSE)
  colnames(mat) = barcodeNames$V1
  rownames(mat) = featureNames$V1 # V1 - ENSG, V2 - gene_names
  mat
  }, mc.cores = mc)
  samplesRawData
}

# Fetch Ensembl ID<-->gene name mapping from input data (all samples)
getAnnot <- function(samples) {
  annot <- data.frame()
  for(s in samples) {
    matrixDir = s
    barcodePath <- paste0(matrixDir, "/", "barcodes.tsv.gz")
    featuresPath <- paste0(matrixDir, "/", "features.tsv.gz")
    matrixPath <- paste0(matrixDir, "/", "matrix.mtx.gz")
    mat <- readMM(file = matrixPath)
    featureNames = read.delim(featuresPath,
                              header = FALSE,
                              stringsAsFactors = FALSE)
    barcodeNames = read.delim(barcodePath,
                              header = FALSE,
                              stringsAsFactors = FALSE)
    colnames(mat) = barcodeNames$V1
    rownames(mat) = featureNames$V1 # V1 - Ensembl ID, V2 - gene_names
  
    annot <- rbind(annot, data.frame(Gene.stable.ID = featureNames$V1, Gene.name = featureNames$V2))
  }
 
  annot <- annot %>% distinct
  annot$Gene.name <- annot$Gene.name %>% as.character
  annot$Gene.stable.ID <- annot$Gene.stable.ID %>% as.character
  annot
}