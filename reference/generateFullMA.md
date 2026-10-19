# Function that runs the full scMI pipeline and outputs complete metaObject

Memory-efficient implementation. The sparse expression layer and
metadata are extracted once, the Seurat object is released, and only one
pair at a time is converted to a dense matrix. The defaults return the
same full metaObject as earlier versions (per-pair expression matrices
and leave-one-out results included), so downstream MetaIntegrator
functions such as `filterGenes(isLeaveOneOut = TRUE)`, `forestPlot`,
`calculateScore` and `pooledROCPlot` work unchanged.

## Usage

``` r
generateFullMA(
  seuratObject,
  uniqueSampleID = "sample",
  pairingColumn,
  comparisonGroupColumn,
  groupOne = "Condition",
  groupTwo = "Control",
  thresholdGenes = NULL,
  completePairsOnly = FALSE,
  assay = "RNA",
  maxCores = NULL,
  keepExpr = TRUE,
  runLeaveOneOut = FALSE
)
```

## Arguments

- seuratObject:

  Seurat v5 Object

- uniqueSampleID:

  Column in seurat metadata object of the unique sample ID

- pairingColumn:

  Column in seurat metadata object that samples are paired on

- comparisonGroupColumn:

  Column in seurat metadata object to use for DE comparisons

- groupOne:

  Identity in comparison group to define DEGs

- groupTwo:

  A second identity for comparison, ie control group

- thresholdGenes:

  Can include a vector of genes to include for analysis. If left empty,
  all genes will be considered

- completePairsOnly:

  BOOL where TRUE subsets seuratObject to samples from only the complete
  pairs

- assay:

  assay to use

- maxCores:

  Number of cores passed to MetaIntegrator::runMetaAnalysis. NULL uses
  the MetaIntegrator default

- keepExpr:

  BOOL where FALSE drops the per-pair expression matrices from the
  returned object. Keep TRUE if the object will be used for plots or
  gene scores

- runLeaveOneOut:

  BOOL passed to MetaIntegrator::runMetaAnalysis. TRUE runs the k
  leave-one-out meta-analyses needed by MetaIntegrator::filterGenes.
  FALSE (default) skips them. The pooled results are identical either
  way

## Value

metaObject of the processed analyses
