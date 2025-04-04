
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
#' @return data.frame of differential expression output
#' @export
runscMI = function(seuratObject, uniqueSampleID = "sample", pairingColumn,
                    comparisonGroupColumn, groupOne = "Condition", groupTwo = "Control",
                    thresholdGenes = NULL, completePairsOnly = FALSE, assay = "RNA") {
  deg_out = generateFullMA(seuratObject, uniqueSampleID, pairingColumn, comparisonGroupColumn,
                             groupOne = groupOne, groupTwo = groupTwo,
                             thresholdGenes = thresholdGenes, completePairsOnly = completePairsOnly, assay = assay)$metaAnalysis$pooledResults


  return(deg_out[,c("gene",'effectSize','effectSizeFDR','effectSizePval', "effectSizeBF", "numPairs")])
}

#' Function that runs the full scMI pipeline and outputs complete metaObject
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
#' @return metaObject of the processed analyses
#' @export
generateFullMA = function(seuratObject, uniqueSampleID = "sample", pairingColumn,
                            comparisonGroupColumn, groupOne = "Condition",
                            groupTwo = "Control", thresholdGenes = NULL, completePairsOnly = FALSE,
                          assay = "RNA") {
  seuratObject@meta.data$groupVar = seuratObject@meta.data[[comparisonGroupColumn]]
  seuratObject = subset(seuratObject, subset = groupVar %in% c(groupOne, groupTwo))
  message("Generating full MA object")
  if (completePairsOnly){
    seuratObject = completePairSubset(seuratObject, uniqueSampleID, pairingColumn)
  }
  if(!confirmPairs(seuratObject, uniqueSampleID, pairingColumn)) {
    stop("Improper sample pairs. Check samples and pairing structure.")
  }
  seuratObject@meta.data$pairingVar = seuratObject@meta.data[[pairingColumn]]

  # Get unique pairs from the subset
  pairList = unique(seuratObject$pairingVar)
  seu_list = list()

  # Subset the data for each pairs
  for (pairOne in pairList) {
    seu_list[[pairOne]] = subset(seuratObject, subset = pairingVar == pairOne)
  }

  names(seu_list) = pairList

  # Convert each subset into MI object
  seu_list_MA = lapply(seu_list, seuToMI, thresholdGenes = thresholdGenes, groupOne = groupOne, groupTwo = groupTwo, comparisonGroupColumn = comparisonGroupColumn, assay = assay)

  message("Running meta-analysis")

  # Format the MI objects and run meta-analysis
  if (.Platform$OS.type == "windows") {
    MI_obj_full = formatMI(seu_list_MA)
    MI_obj_full = runMetaAnalysis(MI_obj_full, maxCores = 1)
  }
  if (.Platform$OS.type != "windows") {
    MI_obj_full = formatMI(seu_list_MA) %>% runMetaAnalysis()
  }
  MI_obj_full$metaAnalysis$pooledResults$gene = rownames(MI_obj_full$metaAnalysis$pooledResults)
  MI_obj_full$metaAnalysis$pooledResults$numPairs = MI_obj_full$metaAnalysis$pooledResults$numStudies
  MI_obj_full$metaAnalysis$pooledResults$method = "scMI"
  MI_obj_full$metaAnalysis$pooledResults$effectSizeBF = p.adjust(MI_obj_full$metaAnalysis$pooledResults$effectSizePval, method = "bonferroni")
  return(MI_obj_full)
}

#' Function that convert seuratObject to MetaIntegrator datasetObject
#'
#' @param seuratObject Seurat v5 Object
#' @param thresholdGenes Can include a vector of genes to include for analysis. If left empty, all genes will be considered
#' @param groupOne Identity in comparison group to define DEGs
#' @param groupTwo A second identity for comparison, ie control group
#' @param comparisonGroupColumn Column in seurat metadata object to use for DE comparisons
#' @param assay assay to use
#' @return seuratObject of data with complete pairs
#'
seuToMI = function(seuratObject, thresholdGenes = NULL, groupOne, groupTwo, comparisonGroupColumn = 'condition', assay = "RNA") {
  dat_id = unique(seuratObject$pairingVar)
  # Initialize the MI object from the original data
  MI_obj1 = tinyMetaObject$originalData$PBMC.Study.1
  MI_obj1$pheno = seuratObject@meta.data

  # Extract expression data, filtering by genes to keep if specified
  if (!is.null(thresholdGenes)) {
    expr = as.matrix(seuratObject@assays[[assay]]$data[rownames(seuratObject@assays[[assay]]$data) %in% thresholdGenes, ]) %>% suppressWarnings()
  } else {
    expr = as.matrix(seuratObject@assays[[assay]]$data) %>% suppressWarnings()
  }

  # Assign expression and class to MI object
  MI_obj1$expr = expr
  MI_obj1$class = ifelse(MI_obj1$pheno[[comparisonGroupColumn]] == groupTwo, 0, NA)
  MI_obj1$class = ifelse(MI_obj1$pheno[[comparisonGroupColumn]] == groupOne, 1, MI_obj1$class)
  names(MI_obj1$class) = rownames(MI_obj1$pheno)

  # Set keys for the expression data
  MI_obj1$keys = rownames(MI_obj1$expr)
  MI_obj1$formattedName = dat_id
  names(MI_obj1$keys) = rownames(MI_obj1$expr)

  # Validate the data object
  checkDataObject(MI_obj1, "Dataset")

  # Prepare and validate the output meta object
  outMetaObj = list()
  outMetaObj$originalData = list(dat1 = MI_obj1)
  if (!checkDataObject(outMetaObj, "Meta", "Pre-Analysis")) {
    stop("Error with datasetObject!")
  }

  return(outMetaObj)
}


#' Function that converts a list of MetaIntegrator datasetObjects into a single MetaIntegrator metaObject
#'
#' @param listOfMetaIntegratorObjs list of MetaIntegrator datasetObjects
#' @return MetaIntegrator metaObject
#'
formatMI = function(listOfMetaIntegratorObjs) {
  completeMetaObject = list(originalData = list())
  for (i in seq_along(listOfMetaIntegratorObjs)) {
    completeMetaObject$originalData[[i]] = listOfMetaIntegratorObjs[[i]]$originalData$dat1
    names(completeMetaObject$originalData)[i] = names(listOfMetaIntegratorObjs)[i]
    completeMetaObject$originalData[[i]]$expr %<>% as.matrix()
  }

  # Check data object validity
  if (!checkDataObject(completeMetaObject, "Meta", "Pre-Analysis")) {
    stop("Error with metaObject!")
  }
  return(completeMetaObject)
}

