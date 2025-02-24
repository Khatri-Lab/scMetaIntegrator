#' Function that runs detects complete pairs
#'
#' @param seuratObject Seurat v5 Object
#' @param uniqueSampleID Column in seurat metadata object of the unique sample ID
#' @param pairingColumn Column in seurat metadata object that samples are paired on
#' @return boolean TRUE if two samples per pair present, FALSE if not
#'
confirmPairs = function(seuratObject, uniqueSampleID = "sample", pairingColumn) {
  seuratObject@meta.data$pairing_column_fxn = seuratObject@meta.data[[pairingColumn]]
  seuratObject@meta.data$sample_fxn = seuratObject@meta.data[[uniqueSampleID]]
  meta = seuratObject@meta.data[, c("sample_fxn", "pairing_column_fxn")]
  rownames(meta) = NULL
  meta = unique(meta)
  return(all(table(meta$pairing_column_fxn) == 2))
}
