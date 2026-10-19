
#' Function that runs the full scMI pipeline and outputs simplified DE results
#'
#' @param seuratObject Seurat v5 Object
#' @param uniqueSampleID Column in seurat metadata object of the unique sample ID
#' @param pairingColumn Column in seurat metadata object that samples are paired on
#' @param comparisonGroupColumn Column in seurat metadata object to use for DE comparisons
#' @param groupOne Identity in comparison group to define DEGs
#' @param groupTwo A second identity for comparison, ie control group
#' @param thresholdGenes Can include a vector of genes to include for analysis. If left empty, all genes will be considered
#' @param completePairsOnly BOOL where TRUE subsets seuratObject to samples from only the complete pairs
#' @param assay assay to use
#' @param maxCores Number of cores passed to MetaIntegrator::runMetaAnalysis. NULL uses the MetaIntegrator default
#' @param runLeaveOneOut BOOL passed to MetaIntegrator::runMetaAnalysis. FALSE (default) skips the extra leave-one-out meta-analyses.
#' @return data.frame of differential expression output
#' @export
runscMI = function(seuratObject, uniqueSampleID = "sample", pairingColumn,
                   comparisonGroupColumn, groupOne = "Condition", groupTwo = "Control",
                   thresholdGenes = NULL, completePairsOnly = FALSE, assay = "RNA",
                   maxCores = NULL, runLeaveOneOut = FALSE) {
  deg_out = generateFullMA(seuratObject, uniqueSampleID, pairingColumn, comparisonGroupColumn,
                           groupOne = groupOne, groupTwo = groupTwo,
                           thresholdGenes = thresholdGenes, completePairsOnly = completePairsOnly,
                           assay = assay, maxCores = maxCores, keepExpr = FALSE,
                           runLeaveOneOut = runLeaveOneOut)$metaAnalysis$pooledResults
  
  return(deg_out[, c("gene", "effectSize", "effectSizeFDR", "effectSizePval", "effectSizeBF", "numPairs")])
}

#' Function that runs the full scMI pipeline and outputs complete metaObject
#'
#' Memory-efficient implementation. The sparse expression layer and metadata are extracted
#' once, the Seurat object is released, and only one pair at a time is converted to a dense
#' matrix. The defaults return the same full metaObject as earlier versions (per-pair
#' expression matrices and leave-one-out results included), so downstream MetaIntegrator
#' functions such as \code{filterGenes(isLeaveOneOut = TRUE)}, \code{forestPlot},
#' \code{calculateScore} and \code{pooledROCPlot} work unchanged.
#'
#' @param seuratObject Seurat v5 Object
#' @param uniqueSampleID Column in seurat metadata object of the unique sample ID
#' @param pairingColumn Column in seurat metadata object that samples are paired on
#' @param comparisonGroupColumn Column in seurat metadata object to use for DE comparisons
#' @param groupOne Identity in comparison group to define DEGs
#' @param groupTwo A second identity for comparison, ie control group
#' @param thresholdGenes Can include a vector of genes to include for analysis. If left empty, all genes will be considered
#' @param completePairsOnly BOOL where TRUE subsets seuratObject to samples from only the complete pairs
#' @param assay assay to use
#' @param maxCores Number of cores passed to MetaIntegrator::runMetaAnalysis. NULL uses the MetaIntegrator default
#' @param keepExpr BOOL where FALSE drops the per-pair expression matrices from the returned object.
#'   Keep TRUE if the object will be used for plots or gene scores
#' @param runLeaveOneOut BOOL passed to MetaIntegrator::runMetaAnalysis. TRUE runs the k
#'   leave-one-out meta-analyses needed by MetaIntegrator::filterGenes. FALSE (default) skips them.
#'   The pooled results are identical either way
#' @return metaObject of the processed analyses
#' @export
generateFullMA = function(seuratObject, uniqueSampleID = "sample", pairingColumn,
                          comparisonGroupColumn, groupOne = "Condition",
                          groupTwo = "Control", thresholdGenes = NULL, completePairsOnly = FALSE,
                          assay = "RNA", maxCores = NULL, keepExpr = TRUE,
                          runLeaveOneOut = FALSE) {
  seuratObject@meta.data$groupVar = seuratObject@meta.data[[comparisonGroupColumn]]
  seuratObject = subset(seuratObject, subset = groupVar %in% c(groupOne, groupTwo))
  message("Generating full MA object")
  if (completePairsOnly){
    seuratObject = completePairSubset(seuratObject, uniqueSampleID, pairingColumn)
  }
  if(!confirmPairs(seuratObject, uniqueSampleID, pairingColumn)) {
    stop("Improper sample pairs. Check samples and pairing structure.")
  }
  
  # Extract metadata and the sparse expression layer once, then release the Seurat object
  meta_data = seuratObject@meta.data
  meta_data$pairingVar = meta_data[[pairingColumn]]
  expr_sparse = seuratObject@assays[[assay]]$data
  if (!is.null(thresholdGenes)) {
    expr_sparse = expr_sparse[rownames(expr_sparse) %in% thresholdGenes, , drop = FALSE]
  }
  rm(seuratObject)
  invisible(gc())
  
  # Build each pair's MetaIntegrator dataset object directly into the metaObject
  MI_template = MetaIntegrator::tinyMetaObject$originalData$PBMC.Study.1
  pairList = unique(meta_data$pairingVar)
  MI_obj_full = list(originalData = list())
  for (pairOne in pairList) {
    pair_cells = rownames(meta_data)[meta_data$pairingVar == pairOne]
    
    MI_obj1 = MI_template
    MI_obj1$pheno = meta_data[pair_cells, , drop = FALSE]
    MI_obj1$expr = suppressWarnings(as.matrix(expr_sparse[, pair_cells, drop = FALSE]))
    
    MI_obj1$class = ifelse(MI_obj1$pheno[[comparisonGroupColumn]] == groupTwo, 0, NA)
    MI_obj1$class = ifelse(MI_obj1$pheno[[comparisonGroupColumn]] == groupOne, 1, MI_obj1$class)
    names(MI_obj1$class) = rownames(MI_obj1$pheno)
    
    MI_obj1$keys = rownames(MI_obj1$expr)
    names(MI_obj1$keys) = rownames(MI_obj1$expr)
    MI_obj1$formattedName = as.character(pairOne)
    
    if (!MetaIntegrator::checkDataObject(MI_obj1, "Dataset")) {
      stop("Error with datasetObject for pair: ", pairOne)
    }
    MI_obj_full$originalData[[as.character(pairOne)]] = MI_obj1
    rm(MI_obj1)
  }
  rm(expr_sparse)
  invisible(gc())
  
  if (!MetaIntegrator::checkDataObject(MI_obj_full, "Meta", "Pre-Analysis")) {
    stop("Error with metaObject!")
  }
  
  message("Running meta-analysis")
  if (.Platform$OS.type == "windows") {
    MI_obj_full = MetaIntegrator::runMetaAnalysis(MI_obj_full, runLeaveOneOutAnalysis = runLeaveOneOut,
                                                  maxCores = 1)
  } else if (is.null(maxCores)) {
    MI_obj_full = MetaIntegrator::runMetaAnalysis(MI_obj_full, runLeaveOneOutAnalysis = runLeaveOneOut)
  } else {
    MI_obj_full = MetaIntegrator::runMetaAnalysis(MI_obj_full, runLeaveOneOutAnalysis = runLeaveOneOut,
                                                  maxCores = maxCores)
  }
  
  MI_obj_full$metaAnalysis$pooledResults$gene = rownames(MI_obj_full$metaAnalysis$pooledResults)
  MI_obj_full$metaAnalysis$pooledResults$numPairs = MI_obj_full$metaAnalysis$pooledResults$numStudies
  MI_obj_full$metaAnalysis$pooledResults$method = "scMI"
  MI_obj_full$metaAnalysis$pooledResults$effectSizeBF = stats::p.adjust(MI_obj_full$metaAnalysis$pooledResults$effectSizePval,
                                                                 method = "bonferroni")
  
  if (!keepExpr) {
    for (pair_name in names(MI_obj_full$originalData)) {
      MI_obj_full$originalData[[pair_name]]$expr = NULL
    }
  }
  
  return(MI_obj_full)
}