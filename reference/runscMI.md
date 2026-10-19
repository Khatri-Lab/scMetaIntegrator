# Function that runs the full scMI pipeline and outputs simplified DE results

Function that runs the full scMI pipeline and outputs simplified DE
results

## Usage

``` r
runscMI(
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

- runLeaveOneOut:

  BOOL passed to MetaIntegrator::runMetaAnalysis. FALSE (default) skips
  the extra leave-one-out meta-analyses.

## Value

data.frame of differential expression output
