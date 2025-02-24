#' Function that runs thresholding of genes by sample pairs
#'
#' @param seuratObject Bayesian Meta-Analysis object with basic bayesian meta-analysis run
#' @param uniqueSampleID Column in seurat metadata object of the unique sample ID
#' @param pairingColumn Column in seurat metadata object that samples are paired on
#' @param completePairsOnly BOOL where TRUE subsets seuratObject to samples from only the complete pairs
#' @param percentThreshold Percentage of cells that should express the gene across each pair to be included
#' @return Vector of genes that pass thresholds
#'
thresholdGenesAcrossPairs = function(seuratObject, uniqueSampleID = "sample",
                                     pairingColumn, percentThreshold = 0.01,
                                     completePairsOnly = FALSE) {
  if (completePairsOnly){
    seuratObject = completePairSubset(seuratObject, uniqueSampleID, pairingColumn)
  }
  if(!confirmPairs(seuratObject, uniqueSampleID, pairingColumn)) {
    stop("Improper sample pairs. Check samples and pairing structure.")
  }
  seuratObject@meta.data$pairing_column_fxn = seuratObject@meta.data[[pairingColumn]]
  pairList = unique(seuratObject@meta.data$pairing_column_fxn)
  pairGenes <- list()
  for (pairOne in pairList) {
    pairOneSrt <- subset(seuratObject, subset = pairing_column_fxn == pairOne)
    pairOneSrt <- as.matrix(pairOneSrt@assays$RNA$counts)
    pairOneGenes <- pairOneSrt[rowSums(pairOneSrt > 0) > (percentThreshold * ncol(pairOneSrt)), , drop = FALSE]
    pairGenes[[pairOne]] <- rownames(pairOneGenes)
  }
  return(Reduce(intersect, pairGenes))
}
