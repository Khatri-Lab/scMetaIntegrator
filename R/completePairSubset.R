#' Function that runs detects complete pairs
#'
#' @param seuratObject Seurat v5 Object
#' @param uniqueSampleID Column in seurat metadata object of the unique sample ID
#' @param pairingColumn Column in seurat metadata object that samples are paired on
#' @return seuratObject of data with complete pairs
#' @export
completePairSubset = function(seuratObject, uniqueSampleID = "sample", pairingColumn) {
  seuratObject@meta.data$pairing_column_fxn = seuratObject@meta.data[[pairingColumn]]
  seuratObject@meta.data$sample_fxn = seuratObject@meta.data[[uniqueSampleID]]
  meta = seuratObject@meta.data[, c("sample_fxn", "pairing_column_fxn")]
  rownames(meta) = NULL
  meta = unique(meta)
  pairingCounts = table(meta$pairing_column_fxn)
  completePairs <- names(pairingCounts[pairingCounts == 2])
  message(paste0("Removing the following samples: ", names(pairingCounts[pairingCounts != 2])))
  if (length(completePairs) == 0) {
    stop("No complete pairs.")
  }
  seuratObject = subset(seuratObject, subset = pairing_column_fxn %in% completePairs)
  seuratObject@meta.data$pairing_column_fxn = NULL
  seuratObject@meta.data$sample_fxn = NULL
  return(seuratObject)
}
